"""
Calibration.jl — `run_calibration(run)` wraps ClimaCalibrate UKI for one run.

This is the function form of the historical `run_calibration_wLabile.jl`. The
heavy model code is NOT duplicated: we `@everywhere include` the original
`model_interface_wLabile.jl` and `site_metadata.jl` from the source experiment
folder, then `@everywhere include` `labile_pin_shim.jl` which defines a wrapper
interface `NeonPipelineInterface` serving both labile modes.

Key differences from the original script:
  - priors are variable-length: built from `run.priors` (the params actually
    given), in canonical order, so EKP calibrates exactly that set;
  - labile pinning: when `labile_depth_scale` is NOT calibrated, the shim injects
    `labile_depth_scale = 0.0` into each per-member TOML, recovering the plain
    model (`exp(0·z)=1`);
  - returns `(; output_dir, eki_path, final_params, param_names)` instead of
    writing them to ENV.

`run_calibration` runs in `Main` (it uses `@everywhere`/Distributed), so this
file is `include`d into Main by the driver, not loaded as a module.
"""

using Dates
using Distributed
import Random
import ClimaCalibrate
import ClimaLand
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2
import TOML
using LinearAlgebra

# Path to the original experiment folder whose heavy scripts we reuse.
const SRC_EXPERIMENT_DIR =
    joinpath(pkgdir(ClimaLand), "experiments", "calibrate_neon")
const PIPELINE_DIR = @__DIR__

# ── Model interface + shim, loaded at TOP LEVEL (older world age) ─────────────
# ClimaCalibrate calls observation_map on the MAIN process and forward_model on
# the workers. Both methods (and the NeonPipelineInterface dispatch tag) must be
# defined in a world age OLDER than the calibrate() call — otherwise we get
# "method too new to be called from this world context" / "not implemented".
#
# So we include the model interface + shim HERE, at include time, on the main
# process. The interface methods only READ the worker globals (SITE_ID,
# OUTPUT_DIR, …) when CALLED, not at include time, so it's fine that those
# globals don't exist yet — run_calibration sets them before calibrate() runs.
# The same files are broadcast to the workers inside run_calibration.
include(joinpath(SRC_EXPERIMENT_DIR, "site_metadata.jl"))
include(joinpath(SRC_EXPERIMENT_DIR, "model_interface_wLabile.jl"))
include(joinpath(PIPELINE_DIR, "labile_pin_shim.jl"))

# Declare the config globals the model interface reads, at TOP LEVEL, so their
# BINDINGS exist in an old world age. run_calibration only *assigns* them (a new
# value is fine; a newly-created binding would be "too new" to read mid-call —
# which is exactly the UndefVarError observation_map hit on the main process).
global SITE_ID = ""
global OUTPUT_DIR = ""
global OBS_FILEPATH = ""
global DT = 900.0
global Caldepthnum = "0.00"
global LABILE_ON = false

# ── Prior construction (variable-length) ─────────────────────────────────────
"Build a combined prior distribution from an ordered Vector{Pair{name,Prior}}."
function build_priors(prior_pairs)
    dists = [
        PD.constrained_gaussian(name, p.mean, p.std, p.lower, p.upper)
        for (name, p) in prior_pairs
    ]
    return PD.combine_distributions(dists)
end

# ── Provenance snapshot ──────────────────────────────────────────────────────
# Copy the code ACTUALLY used by this run into <output_dir>/model_scripts:
#   - pipeline/      : the pipeline step functions that drive the run
#   - calibrate_neon/: the original model_interface_wLabile.jl + site_metadata.jl
#                      that the pipeline includes on the workers (genuinely used)
#   - Biogeochemistry: the soil-CO₂ source code
#   - the config TOML that produced this run
function _snapshot_model_scripts(output_dir)
    scripts_dst = joinpath(output_dir, "model_scripts")
    ispath(scripts_dst) && rm(scripts_dst; recursive = true, force = true)
    mkpath(scripts_dst)
    cp(PIPELINE_DIR, joinpath(scripts_dst, "pipeline"); force = true)
    cp(SRC_EXPERIMENT_DIR, joinpath(scripts_dst, "calibrate_neon"); force = true)
    bio_src = joinpath(pkgdir(ClimaLand), "src/standalone/Soil/Biogeochemistry")
    cp(bio_src, joinpath(scripts_dst, "Biogeochemistry"); force = true)
    if !isempty(Config.CONFIG_PATH) && isfile(Config.CONFIG_PATH)
        cp(Config.CONFIG_PATH, joinpath(scripts_dst, basename(Config.CONFIG_PATH));
           force = true)
    end
    return nothing
end

# ── Main entry ───────────────────────────────────────────────────────────────
"""
    run_calibration(run; obs_filepath, output_dir)

Run the UKI calibration. `obs_filepath` is the observations.jld2 from the
observation step; `output_dir` is the run output directory (created here).
Returns `(; output_dir, eki_path, final_params, param_names)`, where
`final_params` is a `Dict{String,Float64}` of posterior means for the
calibrated parameters.
"""
function run_calibration(run; obs_filepath, output_dir)
    param_names = first.(run.priors)
    labile_on = Config.LABILE_PARAM in param_names
    prior = build_priors(run.priors)

    DT_local = run.dt
    caldepth_str = string(run.cal_depth)

    PRIOR_VALUES_FILE = joinpath(output_dir, "prior_values.txt")
    FINAL_PARAMS_FILE = joinpath(output_dir, "final_parameter_means.txt")

    mkpath(output_dir)
    _snapshot_model_scripts(output_dir)

    open(PRIOR_VALUES_FILE, "w") do io
        println(io, "# Calibrated parameters [mean, std, lower, upper]:")
        for (name, p) in run.priors
            println(io, "$name = [$(p.mean), $(p.std), $(p.lower), $(p.upper)]")
        end
        println(io, "labile_calibrated = $labile_on")
        println(io, "Means only:")
        for (name, p) in run.priors
            println(io, "$name = $(p.mean)")
        end
    end

    # ── Load observations ────────────────────────────────────────────────────
    obs_data = JLD2.load(obs_filepath)
    y_obs = obs_data["y_obs"]
    noise_cov = obs_data["noise_cov"]
    obs_dates = obs_data["obs_dates"]
    n_obs = length(obs_dates)
    println("Loaded $n_obs valid observation days for site $(run.site)")

    # ── Create EKP ───────────────────────────────────────────────────────────
    rng = Random.MersenneTwister(1234)
    ekp = EKP.EnsembleKalmanProcess(
        y_obs, noise_cov,
        EKP.Unscented(prior; α_reg = 1.0, update_freq = 1);
        rng = rng,
    )
    N_ens = EKP.get_N_ens(ekp)
    println("Ensemble size: $N_ens (parameters: $param_names)")

    # ── Workers ──────────────────────────────────────────────────────────────
    curr_backend = ClimaCalibrate.get_backend()
    if curr_backend == ClimaCalibrate.CaltechHPCBackend ||
       curr_backend == ClimaCalibrate.GCPBackend
        @info "Check your slurm script: #tasks must equal ensemble size ($N_ens)"
        addprocs(ClimaCalibrate.SlurmManager())
    else
        current_workers = nprocs() == 1 ? 0 : nworkers()
        n_add = max(0, N_ens - current_workers)
        if n_add > 0
            @info "No SLURM detected — adding $n_add local workers for $N_ens members"
            addprocs(n_add)
        end
    end

    # ── Set config globals + load model code on workers ──────────────────────
    # The model interface + shim are ALREADY loaded on the main process at
    # top-level include (so observation_map, which runs on main, is in an old
    # world age). Here we (a) set the config globals everywhere — the interface
    # reads them when CALLED: OUTPUT_DIR/OBS_FILEPATH on main for observation_map,
    # all of them on workers for forward_model — and (b) load the model code on
    # the WORKERS (fresh processes that don't have it yet).
    #
    # NOTE: `@everywhere`/`@eval @everywhere` expand to *toplevel* expressions and
    # cannot be used inside a function body; `Distributed.remotecall_eval` +
    # `Base.include(Main, path)` are plain calls that are valid here.
    all_procs = procs()           # includes the driver process (pid 1)
    worker_procs = workers()      # the ensemble workers only (pids ≥ 2)
    start_str = string(run.start_date)
    stop_str = string(run.stop_date)

    # (a) config globals on every process (main needs OUTPUT_DIR/OBS_FILEPATH).
    globals_expr = quote
        global SITE_ID = $(run.site)
        global OUTPUT_DIR = $output_dir
        global OBS_FILEPATH = $obs_filepath
        global DT = $DT_local
        global Caldepthnum = $caldepth_str
        global LABILE_ON = $labile_on
        global THETA_R = $(run.theta_r)
        ENV["NEON_START_DATE"] = $start_str
        ENV["NEON_STOP_DATE"] = $stop_str
        nothing
    end
    Distributed.remotecall_eval(Main, all_procs, globals_expr)

    # (b) load model interface + shim on the WORKERS (main already has them).
    if !isempty(worker_procs)
        load_expr = quote
            import Distributed
            import ClimaLand
            import TOML
            Base.include(Main, $(joinpath(SRC_EXPERIMENT_DIR, "site_metadata.jl")))
            Base.include(Main, $(joinpath(SRC_EXPERIMENT_DIR, "model_interface_wLabile.jl")))
            Base.include(Main, $(joinpath(PIPELINE_DIR, "labile_pin_shim.jl")))
            nothing
        end
        Distributed.remotecall_eval(Main, worker_procs, load_expr)
    end

    # ── Run calibration ──────────────────────────────────────────────────────
    println("\n=== Starting ClimaCalibrate calibration ===")
    println("  Site: $(run.site)  |  Period: $(run.start_date)..$(run.stop_date)")
    println("  Iterations: $(run.n_iterations)  |  Ensemble: $N_ens  |  Obs days: $n_obs")
    println("  Parameters: $param_names  (labile_on = $labile_on)")
    println("  Output: $output_dir")

    eki = ClimaCalibrate.calibrate(
        ClimaCalibrate.WorkerBackend(),
        ekp,
        NeonPipelineInterface(),
        run.n_iterations,
        prior,
        output_dir,
    )

    # ── Results ──────────────────────────────────────────────────────────────
    final_vals = EKP.get_ϕ_mean_final(prior, eki)
    final_params = Dict{String, Float64}(
        name => val for (name, val) in zip(param_names, final_vals)
    )

    open(FINAL_PARAMS_FILE, "w") do io
        println(io, "=== Calibration Complete ===")
        println(io, "Site: $(run.site)")
        println(io, "Period: $(run.start_date) .. $(run.stop_date)")
        println(io, "Final parameter means:")
        for name in param_names
            println(io, "  $name = $(round(final_params[name], sigdigits = 5))")
        end
    end
    println("Saved final parameter means to $FINAL_PARAMS_FILE")

    eki_path = _find_eki_path(output_dir)

    # Free workers between runs so a batch loop doesn't accumulate procs.
    rmprocs(workers())

    return (; output_dir, eki_path, final_params, param_names)
end

"Locate the eki_file.jld2 from the last completed iteration."
function _find_eki_path(output_dir)
    iter_dirs = filter(
        d -> isdir(joinpath(output_dir, d)) && startswith(d, "iteration_"),
        readdir(output_dir),
    )
    isempty(iter_dirs) && return nothing
    last_iter = maximum(sort(iter_dirs))
    p = joinpath(output_dir, last_iter, "eki_file.jld2")
    return isfile(p) ? p : nothing
end
