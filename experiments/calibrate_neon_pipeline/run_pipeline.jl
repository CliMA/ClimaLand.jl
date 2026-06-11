"""
run_pipeline.jl — the single entry point for the NEON calibration pipeline.

Usage:
    julia --project=.buildkite experiments/calibrate_neon_pipeline/run_pipeline.jl CONFIG.toml [options]

Options:
    --site SITE            only run [[runs]] entries for this site
    --year YEAR            only run entries whose start year == YEAR
    --run-identifier ID    override the run identifier (default: pipeline-start datetime)
    --stop-on-error        halt on the first failing run (default: log & continue)

Examples:
    # one config, all its runs (single or batch):
    julia ... run_pipeline.jl config/cper_allyears.toml
    # just one (site, year) out of a batch config:
    julia ... run_pipeline.jl config/cper_allyears.toml --site NEON-cper --year 2021

Every step is a pure function (no ENV interface, no top-level `const`s), so a
whole batch runs safely in this one Julia session — no subprocess-per-run needed.
"""

using Dates

# ── Argument parsing ─────────────────────────────────────────────────────────
function parse_args(args)
    isempty(args) && error("Usage: run_pipeline.jl CONFIG.toml [--site S] [--year Y] " *
                           "[--run-identifier ID] [--stop-on-error]")
    config_path = args[1]
    opts = Dict{String, Any}(
        "site" => nothing, "year" => nothing,
        "run_identifier" => nothing, "stop_on_error" => false,
    )
    i = 2
    while i <= length(args)
        a = args[i]
        if a == "--site"
            opts["site"] = args[i + 1]; i += 2
        elseif a == "--year"
            opts["year"] = parse(Int, args[i + 1]); i += 2
        elseif a == "--run-identifier"
            opts["run_identifier"] = args[i + 1]; i += 2
        elseif a == "--stop-on-error"
            opts["stop_on_error"] = true; i += 1
        else
            error("Unknown argument: $a")
        end
    end
    return config_path, opts
end

config_path, opts = parse_args(ARGS)

# ── Load the pipeline machinery (into Main) ──────────────────────────────────
include(joinpath(@__DIR__, "pipeline", "Pipeline.jl"))

# ── Build + filter config ────────────────────────────────────────────────────
cfg = Config.load_config(config_path; run_identifier = opts["run_identifier"])

runs = cfg.runs
opts["site"] !== nothing && (runs = filter(r -> r.site == opts["site"], runs))
opts["year"] !== nothing &&
    (runs = filter(r -> year(r.start_date) == opts["year"], runs))
isempty(runs) && error("No runs left after applying --site/--year filters")

cfg_filtered = Config.PipelineConfig(cfg.results_csv, cfg.steps, runs)

println("Loaded config:   $config_path")
println("Run identifier:  $(first(runs).run_identifier)")
println("Runs to execute: $(length(runs))")
for r in runs
    rstr = r.restart_from === nothing ? "" : "  (restart_from=$(r.restart_from))"
    println("  - $(r.site)  $(r.start_date)..$(r.stop_date)  " *
            "params=$(Config.calibrated_param_names(r))$rstr")
end

# ── Go ───────────────────────────────────────────────────────────────────────
run_pipeline(cfg_filtered; stop_on_error = opts["stop_on_error"])
