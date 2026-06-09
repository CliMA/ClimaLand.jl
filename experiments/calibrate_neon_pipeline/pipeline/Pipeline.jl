"""
Pipeline.jl — orchestrates one or many runs IN A SINGLE SESSION.

Each step is a pure function that takes the `RunConfig` (and prior step results)
and RETURNS its results — no ENV between steps, no top-level `const`s, so a batch
loop can run many runs in one Julia process without stale-global problems.

For each run it:
  1. resolves restart (seed prior means from a previous run's posterior);
  2. writes a per-run priors.toml (provenance);
  3. runs the enabled steps via `run_step`, threading return values;
  4. appends one row (priors + posteriors + final RMSE) to the master CSV;
  5. on failure, records a `status=failed` row and (by default) continues.

This file is `include`d into Main by the driver (run_calibration uses Distributed
and the step functions define model objects in Main).
"""

using Dates
import TOML

# Include the layers once into Main, in dependency order.
include(joinpath(@__DIR__, "Config.jl"))
include(joinpath(@__DIR__, "ResultsTable.jl"))
include(joinpath(@__DIR__, "Observations.jl"))
include(joinpath(@__DIR__, "Calibration.jl"))
include(joinpath(@__DIR__, "Diagnostics.jl"))
include(joinpath(@__DIR__, "ForwardRun.jl"))

using .Config
using .ResultsTable

# ── Logging wrapper ──────────────────────────────────────────────────────────
function run_step(step_name, f)
    println("\n=== START: $step_name ===")
    t0 = time()
    try
        result = f()
        println("=== DONE: $step_name ($(round(time() - t0; digits = 2)) s) ===")
        return result
    catch err
        println("=== FAILED: $step_name after $(round(time() - t0; digits = 2)) s ===")
        showerror(stdout, err, catch_backtrace())
        println()
        rethrow()
    end
end

# ── Restart seeding ──────────────────────────────────────────────────────────
"""
Return a new RunConfig whose prior MEANS come from the referenced previous run's
posterior, keeping std/bounds from `run`. Parameters absent from the previous
posterior keep their config prior (so adding a parameter on restart works).
"""
function seed_from_restart(run, csv_path, output_root)
    prev_dir = ResultsTable.find_previous_run(csv_path, output_root, run.restart_from)
    prev_means = ResultsTable.read_posterior_means(
        joinpath(prev_dir, "final_parameter_means.txt"))
    new_priors = map(run.priors) do (name, p)
        haskey(prev_means, name) ?
            (name => Config.Prior(prev_means[name], p.std, p.lower, p.upper)) :
            (name => p)
    end
    println("Restart: seeded prior means from $(run.restart_from)  ($prev_dir)")
    return Config.RunConfig(
        run.site, run.start_date, run.stop_date, run.spinup_days,
        run.n_iterations, run.cal_depth, run.settingsdesc, run.dt,
        run.output_root, run.run_identifier, run.restart_from, new_priors)
end

"""
    dedup_run_identifier(run) -> RunConfig

If the leaf output directory for `run` already exists, return a copy of `run`
whose `run_identifier` has a `_2`, `_3`, … suffix appended until the resulting
`output_dir_for` path is free, so a re-run with the same settingsdesc and
run_identifier does not overwrite a previous run. Returns `run` unchanged if the
original path does not yet exist.
"""
function dedup_run_identifier(run)
    Config.output_dir_for(run) |> ispath || return run
    base_id = run.run_identifier
    n = 2
    local new_run
    while true
        new_run = Config.RunConfig(
            run.site, run.start_date, run.stop_date, run.spinup_days,
            run.n_iterations, run.cal_depth, run.settingsdesc, run.dt,
            run.output_root, "$(base_id)_$(n)", run.restart_from, run.priors)
        ispath(Config.output_dir_for(new_run)) || break
        n += 1
    end
    println("Run identifier '$base_id' already used; using '$(new_run.run_identifier)' to avoid overwrite")
    return new_run
end

# ── priors.toml provenance ───────────────────────────────────────────────────
function _write_priors_toml(run, output_dir)
    mkpath(output_dir)
    path = joinpath(output_dir, "priors.toml")
    open(path, "w") do io
        TOML.print(io, Dict("priors" => Config.priors_to_toml_dict(run)))
    end
    return path
end

# Parameters to run the forward model at. When calibration ran, use its
# posterior means; otherwise (forward-only run) fall back to the config
# `[priors]` means so the forward run reflects the configured parameters rather
# than ForwardRun.jl's hard-coded defaults. In both cases pin
# labile_depth_scale=0.0 when labile was not among the run's parameters.
function _forward_params(run, posterior)
    p = if isempty(posterior)
        Dict{String, Float64}(name => pr.mean for (name, pr) in run.priors)
    else
        Dict{String, Float64}(posterior)
    end
    if !(Config.LABILE_PARAM in Config.calibrated_param_names(run))
        p[Config.LABILE_PARAM] = 0.0
    end
    return p
end

# ── One run ──────────────────────────────────────────────────────────────────
"""
    run_one(run, cfg; skip_restart=false) -> NamedTuple

Execute the enabled steps for a single run and append its CSV row.
"""
function run_one(run, cfg::Config.PipelineConfig; skip_restart = false)
    if !skip_restart && run.restart_from !== nothing
        run = seed_from_restart(run, cfg.results_csv, run.output_root)
    end

    # Avoid overwriting a previous run with the same settingsdesc/run_identifier.
    run = dedup_run_identifier(run)

    output_dir = Config.output_dir_for(run)
    base_dir = dirname(output_dir)
    labile_on = Config.LABILE_PARAM in Config.calibrated_param_names(run)

    println("\n" * "#"^78)
    println("# RUN  site=$(run.site)  $(run.start_date)..$(run.stop_date)  id=$(run.run_identifier)")
    println("#      output_dir = $output_dir")
    println("#"^78)

    mkpath(output_dir)
    _write_priors_toml(run, output_dir)

    obs_filepath = joinpath(base_dir, "observations.jld2")

    # 1. Observations (needed by calibrate + eki diagnostics)
    obs_needed = cfg.steps[:calibrate] || cfg.steps[:plot_eki_diagnostics]
    if cfg.steps[:generate_observations]
        obs_filepath = run_step("Generate observations",
            () -> generate_observations(run; base_dir = base_dir))
    elseif obs_needed && !isfile(obs_filepath)
        error("generate_observations disabled but a step needs $obs_filepath")
    end

    # 2. Calibration
    calib = nothing
    if cfg.steps[:calibrate]
        calib = run_step("Run calibration",
            () -> run_calibration(run; obs_filepath = obs_filepath, output_dir = output_dir))
    end
    posterior = calib === nothing ? Dict{String, Float64}() : calib.final_params
    eki_path = calib === nothing ? nothing : calib.eki_path

    # 3. EKI diagnostics
    final_rmse = missing
    if cfg.steps[:plot_eki_diagnostics] && eki_path !== nothing
        diag = run_step("Plot EKI diagnostics",
            () -> plot_eki_diagnostics(run; output_dir = output_dir,
                eki_path = eki_path, obs_filepath = obs_filepath))
        final_rmse = diag.final_rmse
    end

    # 4. Forward run at optimized params (+ figures + scatter RMSE/corr stats)
    scatter_stats = nothing
    if cfg.steps[:run_prior_mean]
        fwd_params = _forward_params(run, posterior)
        fwd_dir = joinpath(output_dir, "prior_mean_optimized")
        fwd = run_step("Forward run at optimized params",
            () -> forward_run(run; output_dir = fwd_dir, params = fwd_params))
        scatter_stats = get(fwd, :scatter_stats, nothing)
    end

    # 5. Record the row
    ResultsTable.append_row!(cfg.results_csv, run;
        status = "ok", output_dir = output_dir,
        posterior = posterior, final_rmse = final_rmse, labile_pin = 0.0,
        scatter_stats = scatter_stats)
    println("Appended results row to $(cfg.results_csv)")

    return (; output_dir, run, posterior, final_rmse, labile_on)
end

# ── Whole pipeline (single session) ──────────────────────────────────────────
"""
    run_pipeline(cfg; stop_on_error=false)

Run every `RunConfig` in `cfg` in this process. On a run failure: log, append a
`status=failed` row, and (default) continue; `stop_on_error=true` re-throws.
Because every step is a pure function (no `const`s, no ENV interface), a batch
runs safely in one session. Returns (run, status) tuples.
"""
function run_pipeline(cfg::Config.PipelineConfig; stop_on_error = false)
    results = Tuple{Config.RunConfig, String}[]
    for (i, run) in enumerate(cfg.runs)
        println("\n", "="^78)
        println("PIPELINE RUN $i / $(length(cfg.runs))")
        println("="^78)
        try
            run_one(run, cfg)
            push!(results, (run, "ok"))
        catch err
            @error "Run failed" site = run.site start = run.start_date exception = err
            try
                ResultsTable.append_row!(cfg.results_csv, run;
                    status = "failed", output_dir = Config.output_dir_for(run))
            catch csv_err
                @error "Also failed to record failure row" exception = csv_err
            end
            push!(results, (run, "failed"))
            stop_on_error && rethrow()
        end
    end

    println("\n", "="^78)
    println("PIPELINE COMPLETE")
    for (run, status) in results
        println("  [$status] $(run.site) $(run.start_date)..$(run.stop_date)")
    end
    println("Master CSV: $(cfg.results_csv)")
    println("="^78)
    return results
end
