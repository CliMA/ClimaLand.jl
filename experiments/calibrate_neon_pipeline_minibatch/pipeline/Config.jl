"""
Config.jl — parse a pipeline TOML into a typed `Config` with a resolved list of
`Run`s, and build the priors / output paths each run needs.

MINIBATCH pipeline: a `Run` is one calibration for one SITE over ALL of its
calibration WINDOWS jointly (not one run per year, unlike the AllDepth pipeline
this is based on). The `years = [...]` / `start`/`stop` sugar no longer fans out
into separate runs — every date range in a `[[runs]]` entry becomes one window
in `run.windows`, and EKP draws `minibatch_size` windows per iteration via a
minibatcher (see Observation_flag.jl and Calibration.jl). The driver loops over
`config.runs` (one per `[[runs]]` entry).

Each `Run` carries:
  - the site + the full window list + per-run settings (incl. minibatch knobs),
  - its resolved priors (default `[priors]` merged with any `[runs.priors]`
    override, or seeded from a previous run's posterior via `restart_from`),
  - its `run_identifier` and the computed output directory.

`start_date`/`stop_date` are retained as the SPAN (earliest window start,
latest window stop) for output-dir naming and any code that still reads a single
range; the authoritative calibration windows live in `run.windows`.

The set of *calibrated* parameters for a run is exactly the set of priors
present after merging (see `PARAM_ORDER`). `labile_depth_scale` is the optional
5th parameter: present => calibrated, absent => pinned to k=0 in the model.

See ../PIPELINE_PLAN.md (mirrored from the original experiment folder) for the
full design rationale.
"""
module Config

using Dates
import TOML

export PARAM_ORDER, LABILE_PARAM, PipelineConfig, RunConfig, Prior, load_config,
       output_dir_for, priors_to_toml_dict, calibrated_param_names, obs_depth_code,
       cal_depth_codes_of, parse_cal_depth_codes

# ── Canonical parameter ordering ─────────────────────────────────────────────
# The first four are DAMM soilCO2 params (ClimaParams TOML keys). The fifth,
# labile_depth_scale, is optional and NOT a ClimaParams key — it is read
# directly from the parameter TOML by model_interface_wLabile.jl. The model
# always expects a labile_depth_scale entry in the TOML (k=0 when not
# calibrated), but EKP only calibrates the names that have a prior.
const PARAM_ORDER = [
    "soilCO2_reference_rate",
    "soilCO2_activation_energy",
    "michaelis_constant",
    "O2_michaelis_constant",
    "labile_depth_scale",
]
const LABILE_PARAM = "labile_depth_scale"

# Absolute path of the config file currently loaded (set by load_config); copied
# into each run's model_scripts/ for provenance. Empty until a config is loaded.
global CONFIG_PATH = ""

# Valid calibration depths (m) and their NEON observation depth codes.
const VALID_CAL_DEPTHS = Dict(0.02 => "501", 0.06 => "502")

# NEON soil-CO₂ measurement depth codes (shallow → deep). The AllDepth pipeline
# calibrates on an explicit ordered subset of these; the order defines the
# stacking order of the observation vector and G matrix everywhere downstream.
const VALID_DEPTH_CODES = ["501", "502", "503"]

# ── Prior spec ───────────────────────────────────────────────────────────────
"""A single constrained-gaussian prior: mean, std, and bounds. `Inf`/`-Inf`
allowed for bounds (stored as `Float64`, written to TOML as the strings
\"Inf\"/\"-Inf\" because TOML cannot represent infinity)."""
struct Prior
    mean::Float64
    std::Float64
    lower::Float64
    upper::Float64
end

# ── Run config ───────────────────────────────────────────────────────────────
struct RunConfig
    site::String
    # SPAN of the calibration windows (earliest start, latest stop). Kept for
    # output-dir naming and back-compat; the authoritative windows are `windows`.
    start_date::Date
    stop_date::Date
    # MINIBATCH: the full ordered list of calibration windows for this site. EKP
    # draws `minibatch_size` of these per iteration (see Observation_flag.jl).
    # Each entry is one (spinup + window) simulation on the model side.
    windows::Vector{Tuple{Date, Date}}
    spinup_days::Int
    n_iterations::Int
    # MINIBATCH knobs. `minibatch_size` windows per EKP iteration (default 1);
    # `max_n_samples` caps how many windows are kept (0 => keep all); `rng_seed`
    # seeds both the minibatcher and the EKP RNG for reproducibility.
    minibatch_size::Int
    max_n_samples::Int
    rng_seed::Int
    cal_depth::Float64
    # AllDepth: explicit ordered depth-code list to calibrate on jointly (e.g.
    # ["501","502","503"]). `cal_depth` is retained for output-dir naming and any
    # single-depth diagnostics; the actual obs/model depths come from this list.
    cal_depth_codes::Vector{String}
    settingsdesc::String
    dt::Float64
    output_root::String
    run_identifier::String
    restart_from::Union{String, Nothing}
    # ordered (name => Prior) for the params actually calibrated in this run
    priors::Vector{Pair{String, Prior}}
    # optional per-run soil residual water content override (uniform over depth);
    # `nothing` => use the model's default θ_r. When set, it is also appended to
    # the run_identifier (see load_config) so its output dir is distinct.
    theta_r::Union{Float64, Nothing}
end

# ── Pipeline config ──────────────────────────────────────────────────────────
struct PipelineConfig
    results_csv::String           # absolute path to the master CSV
    steps::Dict{Symbol, Bool}
    runs::Vector{RunConfig}
end

# ── Helpers ──────────────────────────────────────────────────────────────────

"Parse a bound that may be the string \"Inf\"/\"-Inf\" or a number."
function _parse_bound(x)
    x isa AbstractString && return x == "Inf" ? Inf :
                                  x == "-Inf" ? -Inf :
                                  parse(Float64, x)
    return Float64(x)
end

"Render a bound for writing back to a TOML file (Inf -> \"Inf\")."
bound_to_toml(x::Real) = isinf(x) ? (x > 0 ? "Inf" : "-Inf") : Float64(x)

"Build a `Prior` from a TOML table `{mean, std, lower, upper}`."
function _prior_from_table(t::AbstractDict)
    return Prior(
        Float64(t["mean"]),
        Float64(t["std"]),
        _parse_bound(get(t, "lower", -Inf)),
        _parse_bound(get(t, "upper", Inf)),
    )
end

"Merge a default prior table with a per-run override table (field-by-field)."
function _merge_prior(default::AbstractDict, override::AbstractDict)
    merged = Dict{String, Any}()
    for k in ("mean", "std", "lower", "upper")
        if haskey(override, k)
            merged[k] = override[k]
        elseif haskey(default, k)
            merged[k] = default[k]
        end
    end
    return merged
end

"""
Resolve the prior set for one run: start from the default `[priors]` table,
apply this run's `[runs.priors]` overrides field-by-field. A parameter is
included iff it appears in the default table OR in the run override (so a run
can ADD a parameter — e.g. switch on `labile_depth_scale` — that the default
omitted). Returns an ordered Vector{Pair{name, Prior}} in PARAM_ORDER.
"""
function _resolve_priors(
    default_priors::AbstractDict,
    run_priors::AbstractDict,
)
    names = String[]
    for name in PARAM_ORDER
        (haskey(default_priors, name) || haskey(run_priors, name)) || continue
        push!(names, name)
    end
    # Any non-canonical names are a config error — fail loudly.
    for name in union(keys(default_priors), keys(run_priors))
        name in PARAM_ORDER ||
            error("Unknown prior parameter \"$name\" (expected one of $(PARAM_ORDER))")
    end

    resolved = Pair{String, Prior}[]
    for name in names
        d = get(default_priors, name, Dict{String, Any}())
        o = get(run_priors, name, Dict{String, Any}())
        merged = _merge_prior(d, o)
        haskey(merged, "mean") && haskey(merged, "std") ||
            error("Prior \"$name\" is missing mean/std after merge")
        push!(resolved, name => _prior_from_table(merged))
    end
    return resolved
end

"""
    _date_ranges(entry) -> Vector{Tuple{Date,Date}}

Expand a [[runs]] entry into the ordered list of calibration windows. Accepts
`years = [...]` (each year → one Jan-1..Dec-31 window) and/or an explicit
`start`/`stop` pair (one window). If both are given, the explicit pair is
appended after the yearly windows. MINIBATCH: the full list is the site's window
set — it is NOT fanned out into separate runs.
"""
function _date_ranges(entry::AbstractDict)
    windows = Tuple{Date, Date}[]
    if haskey(entry, "years")
        for y in entry["years"]
            push!(windows, (Date(y, 1, 1), Date(y, 12, 31)))
        end
    end
    if haskey(entry, "start") && haskey(entry, "stop")
        push!(windows, (Date(entry["start"]), Date(entry["stop"])))
    end
    isempty(windows) &&
        error("Each [[runs]] needs `years` and/or both `start` and `stop`")
    return windows
end

"Validate a calibration depth, returning the matching NEON obs depth code."
function obs_depth_code(cal_depth::Real)
    for (d, code) in VALID_CAL_DEPTHS
        isapprox(Float64(cal_depth), d; atol = 1e-9) && return code
    end
    error("cal_depth=$cal_depth m not supported; must be one of " *
          "$(sort(collect(keys(VALID_CAL_DEPTHS)))) (mapped to NEON obs depth 501/502)")
end

"""
    parse_cal_depth_codes(x) -> Vector{String}

Normalize a `cal_depth_codes` config value into an ordered, validated list of NEON
depth codes. Accepts a TOML array (`["501","502"]` or `[501, 502]`) or a
comma-separated string (`"501,502,503"`). Preserves order, rejects duplicates and
unknown codes, and requires at least one code.
"""
function parse_cal_depth_codes(x)
    raw = x isa AbstractString ? split(x, ",") : x
    codes = String[strip(string(c)) for c in raw]
    isempty(codes) && error("cal_depth_codes must list at least one depth code")
    for c in codes
        c in VALID_DEPTH_CODES ||
            error("Unknown depth code \"$c\"; must be one of $(VALID_DEPTH_CODES)")
    end
    length(unique(codes)) == length(codes) ||
        error("cal_depth_codes has duplicate codes: $(codes)")
    return codes
end

"The ordered depth-code list a run calibrates on (obs/model stacking order)."
cal_depth_codes_of(run::RunConfig) = run.cal_depth_codes

# ── Output path (single source of truth) ─────────────────────────────────────
"""
    output_dir_for(run::RunConfig)

The leaf run directory. MINIBATCH: one SITE = one output dir = one posterior, so
the range folder is named by the span of years plus a `_minibatch` tag instead
of a single `<start>_<stop>`:

    <root>/<site>/<site>_<firstyear>-<lastyear>_minibatch/SpinUP-<spinup>d/CalDepth-<codes>/<n_iter>-It/<settingsdesc>/output_<id>
"""
function output_dir_for(run::RunConfig)
    # AllDepth: name the CalDepth folder by the stacked depth-code list so
    # multi-depth outputs never collide with the old single-depth layout
    # (e.g. "CalDepth-501-502-503").
    depth_tag = "CalDepth-" * join(run.cal_depth_codes, "-")
    y0 = year(run.start_date)
    y1 = year(run.stop_date)
    span_tag = y0 == y1 ? "$(y0)" : "$(y0)-$(y1)"
    return joinpath(
        run.output_root,
        run.site,
        "$(run.site)_$(span_tag)_minibatch",
        "SpinUP-$(run.spinup_days)d",
        depth_tag,
        "$(run.n_iterations)-It",
        run.settingsdesc,
        "output_$(run.run_identifier)",
    )
end

"The base path (everything above `output_<id>`); used for observations.jld2."
function base_dir_for(run::RunConfig)
    return dirname(output_dir_for(run))
end

# ── Priors -> ClimaParams-style TOML dict ────────────────────────────────────
"""
    priors_to_toml_dict(run; means)

Return an `OrderedDict`-like `Dict` mapping each parameter name to its prior
spec (for writing a `priors.toml`). If `means` is provided (name => value), the
mean is overridden with that value — used when seeding from a restart posterior.
"""
function priors_to_toml_dict(run::RunConfig)
    d = Dict{String, Any}()
    for (name, p) in run.priors
        d[name] = Dict(
            "mean" => p.mean,
            "std" => p.std,
            "lower" => bound_to_toml(p.lower),
            "upper" => bound_to_toml(p.upper),
        )
    end
    return d
end

"Names of the parameters actually calibrated in this run, in PARAM_ORDER."
calibrated_param_names(run::RunConfig) = first.(run.priors)

# ── Loading ──────────────────────────────────────────────────────────────────
"""
    load_config(path; run_identifier=nothing)

Parse a pipeline TOML at `path` into a `PipelineConfig`. `run_identifier`
overrides the config/default identifier (the driver passes the pipeline-start
datetime so a whole batch shares one id).
"""
function load_config(path::AbstractString; run_identifier = nothing)
    isfile(path) || error("Config file not found: $path")
    global CONFIG_PATH = abspath(path)   # remembered for provenance snapshots
    raw = TOML.parsefile(path)

    settings = get(raw, "settings", Dict{String, Any}())
    output_root = settings["output_root"]
    results_csv = get(settings, "results_csv", "calibration_results.csv")
    results_csv = isabspath(results_csv) ? results_csv :
                  joinpath(output_root, results_csv)

    # run identifier: CLI arg > config > pipeline-start datetime default
    id = run_identifier !== nothing ? run_identifier :
         get(settings, "run_identifier", default_run_identifier())

    steps_raw = get(raw, "steps", Dict{String, Any}())
    steps = Dict{Symbol, Bool}(
        :generate_observations =>
            Bool(get(steps_raw, "generate_observations", true)),
        :calibrate => Bool(get(steps_raw, "calibrate", true)),
        :plot_eki_diagnostics =>
            Bool(get(steps_raw, "plot_eki_diagnostics", true)),
        :run_prior_mean => Bool(get(steps_raw, "run_prior_mean", true)),
    )

    default_priors = get(raw, "priors", Dict{String, Any}())

    runs = RunConfig[]
    for entry in get(raw, "runs", Any[])
        site = entry["site"]
        # spinup default 20 (mirrors run_from_env); still overridable per-run/settings.
        spinup = Int(get(entry, "spinup_days", get(settings, "spinup_days", 20)))
        n_iter = Int(get(entry, "n_iterations", settings["n_iterations"]))
        cal_depth = Float64(get(entry, "cal_depth", settings["cal_depth"]))
        obs_depth_code(cal_depth)  # validate early
        # AllDepth: explicit ordered depth-code list (entry overrides settings;
        # defaults to all three). Validated here.
        cal_depth_codes = parse_cal_depth_codes(
            get(entry, "cal_depth_codes",
                get(settings, "cal_depth_codes", VALID_DEPTH_CODES)),
        )
        settingsdesc = String(get(entry, "settingsdesc", settings["settingsdesc"]))
        dt = Float64(get(entry, "dt", get(settings, "dt", 900.0)))
        restart_from = haskey(entry, "restart_from") ?
                       String(entry["restart_from"]) : nothing

        # MINIBATCH knobs (entry overrides settings; sane defaults).
        minibatch_size =
            Int(get(entry, "minibatch_size", get(settings, "minibatch_size", 1)))
        max_n_samples =
            Int(get(entry, "max_n_samples", get(settings, "max_n_samples", 0)))
        rng_seed = Int(get(entry, "rng_seed", get(settings, "rng_seed", 1234)))
        minibatch_size >= 1 || error("minibatch_size must be ≥ 1")

        run_priors = get(entry, "priors", Dict{String, Any}())
        priors = _resolve_priors(default_priors, run_priors)

        # Optional per-run θ_r override (uniform over depth); `nothing` keeps the
        # model default. Set `settingsdesc` per-run to keep its output dir
        # distinct, as with the labile runs.
        theta_r = haskey(entry, "theta_r") ? Float64(entry["theta_r"]) : nothing

        # MINIBATCH: gather ALL windows for this entry into ONE run (no fan-out).
        windows = _date_ranges(entry)
        # cap to max_n_samples (keep earliest) when set (>0)
        if max_n_samples > 0 && length(windows) > max_n_samples
            windows = windows[1:max_n_samples]
        end
        span_start = minimum(first.(windows))
        span_stop = maximum(last.(windows))

        push!(runs, RunConfig(
            site, span_start, span_stop, windows, spinup, n_iter,
            minibatch_size, max_n_samples, rng_seed, cal_depth,
            cal_depth_codes, settingsdesc, dt, output_root, id, restart_from,
            priors, theta_r,
        ))
    end

    isempty(runs) && error("Config has no [[runs]] entries")
    return PipelineConfig(results_csv, steps, runs)
end

"Default run identifier = pipeline-start datetime, e.g. 20260603-141530."
default_run_identifier() = Dates.format(now(), "yyyymmdd-HHMMSS")

# ── Standalone / ENV fallback ────────────────────────────────────────────────
"""
    run_from_env()

Build a single `RunConfig` from environment variables, for running a heavy
script by hand (outside the driver). Mirrors the historical ENV contract:
NEON_SITE_ID, NEON_START_DATE, NEON_STOP_DATE, NEON_SPINUP_DAYS,
NEON_N_ITERATIONS, CALL_DEPTH, NEON_SETTINGSDESC, CALL_OUTPUT_PATH (root),
CALL_PRIORS_FILE (a TOML with a [priors] table; falls back to wide defaults).
"""
function run_from_env()
    site = get(ENV, "NEON_SITE_ID", "NEON-srer")
    start_date = Date(get(ENV, "NEON_START_DATE", "2019-01-01"))
    stop_date = Date(get(ENV, "NEON_STOP_DATE", "2019-12-31"))
    spinup = parse(Int, get(ENV, "NEON_SPINUP_DAYS", "20"))
    n_iter = parse(Int, get(ENV, "NEON_N_ITERATIONS", "10"))
    cal_depth = parse(Float64, get(ENV, "CALL_DEPTH", "0.02"))
    cal_depth_codes = parse_cal_depth_codes(get(ENV, "CALL_DEPTHS", "501,502,503"))
    settingsdesc = get(ENV, "NEON_SETTINGSDESC", "standalone")
    dt = parse(Float64, get(ENV, "CALL_DT", "900.0"))
    output_root = get(
        ENV, "CALL_OUTPUT_ROOT",
        "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_calibration",
    )
    id = get(ENV, "NEON_RUN_IDENTIFIER", default_run_identifier())

    prior_file = get(ENV, "CALL_PRIORS_FILE", "")
    default_priors = isfile(prior_file) ?
                     get(TOML.parsefile(prior_file), "priors", _DEFAULT_PRIORS) :
                     _DEFAULT_PRIORS
    priors = _resolve_priors(default_priors, Dict{String, Any}())

    # standalone: a single-window run (minibatch_size 1, no cap).
    windows = [(start_date, stop_date)]
    minibatch_size = parse(Int, get(ENV, "NEON_MINIBATCH_SIZE", "1"))
    max_n_samples = parse(Int, get(ENV, "NEON_MAX_N_SAMPLES", "0"))
    rng_seed = parse(Int, get(ENV, "NEON_RNG_SEED", "1234"))

    return RunConfig(
        site, start_date, stop_date, windows, spinup, n_iter,
        minibatch_size, max_n_samples, rng_seed, cal_depth,
        cal_depth_codes, settingsdesc, dt, output_root, id, nothing, priors,
        nothing,
    )
end

# Historical defaults from run_calibration_wLabile.jl (4-param; labile off).
const _DEFAULT_PRIORS = Dict{String, Any}(
    "soilCO2_reference_rate" =>
        Dict("mean" => 1.0e-7, "std" => 5.0e-8, "lower" => 1.0e-12, "upper" => 1.0e-6),
    "soilCO2_activation_energy" =>
        Dict("mean" => 110550.0, "std" => 40000.0, "lower" => 10.0, "upper" => 200000.0),
    "michaelis_constant" =>
        Dict("mean" => 0.04, "std" => 0.02, "lower" => 0.0, "upper" => "Inf"),
    "O2_michaelis_constant" =>
        Dict("mean" => 0.04, "std" => 0.02, "lower" => 0.0, "upper" => "Inf"),
)

end # module
