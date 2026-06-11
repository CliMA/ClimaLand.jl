"""
ResultsTable.jl — the joint master CSV of priors + posteriors, one appended row
per run.

Design notes:
  - Fixed SUPERSET of columns (all 5 params, including labile_depth_scale) so the
    schema is stable whether or not a run calibrates labile. Uncalibrated params
    get `missing` prior cells; their posterior cell is the pinned value (0.0 for
    labile) or `missing` if no value is known.
  - Header is written once; rows are appended under a file lock so a batch loop
    (or concurrent pipelines) cannot interleave/corrupt rows.
  - `restart_from` + `run_identifier` make a restart chain traceable; the lookup
    helpers here resolve a previous run's posterior means for seeding.
"""
module ResultsTable

using Dates
using CSV
using DataFrames

# Config is included once into Main by Pipeline.jl (or the driver) BEFORE this
# file, so we reference the already-loaded module rather than re-including it
# (which would create a second, incompatible Config type).
import ..Config
using ..Config: PARAM_ORDER, RunConfig, calibrated_param_names, output_dir_for

export append_row!, read_posterior_means, find_previous_run

# ── Column schema ────────────────────────────────────────────────────────────
const META_COLS = [
    "timestamp", "status", "run_identifier", "restart_from",
    "site", "start", "stop", "cal_depth", "spinup_days", "n_iterations",
    "settingsdesc", "output_dir",
]

const PRIOR_FIELDS = ["mean", "std", "lower", "upper"]

# Forward-run scatter/summary diagnostics. Appended after `final_rmse`; all
# default to `missing` when no forward run was done. Each column is named
# `forward_<field>` where <field> is the matching `scatter_stats` NamedTuple
# field (see _fill_scatter_stats!), so the column↔field mapping is mechanical.
const SCATTER_COLS = [
    # RMSE + correlations (obs vs model CO₂/SWC)
    "forward_rmse_sco2", "forward_rmse_swc",
    "forward_corr_obs_model_sco2", "forward_corr_obs_model_swc",
    "forward_corr_obs_sco2_swc", "forward_corr_model_sco2_swc",
    # mean / min / max of daily-mean series, obs and model:
    # soil CO₂ (ppm)
    "forward_obs_sco2_mean", "forward_obs_sco2_min", "forward_obs_sco2_max",
    "forward_model_sco2_mean", "forward_model_sco2_min", "forward_model_sco2_max",
    # soil water content (m³/m³)
    "forward_obs_swc_mean", "forward_obs_swc_min", "forward_obs_swc_max",
    "forward_model_swc_mean", "forward_model_swc_min", "forward_model_swc_max",
    # soil temperature (K)
    "forward_obs_tsoil_mean", "forward_obs_tsoil_min", "forward_obs_tsoil_max",
    "forward_model_tsoil_mean", "forward_model_tsoil_min", "forward_model_tsoil_max",
]

prior_col(param, field) = "prior_$(param)_$(field)"
post_col(param) = "post_$(param)"

"Full ordered column list (fixed superset)."
function column_names()
    cols = copy(META_COLS)
    for p in PARAM_ORDER, f in PRIOR_FIELDS
        push!(cols, prior_col(p, f))
    end
    for p in PARAM_ORDER
        push!(cols, post_col(p))
    end
    push!(cols, "final_rmse")
    append!(cols, SCATTER_COLS)
    return cols
end

# ── Appending ────────────────────────────────────────────────────────────────
"""
    append_row!(csv_path, run; status, output_dir, posterior, final_rmse, labile_pin)

Append one row. `posterior` is a `Dict{String,Float64}` of calibrated parameter
means (may be empty for a failed run). `labile_pin` is the value to record for
labile when it is NOT calibrated (default 0.0). `final_rmse` may be `missing`.
"""
function append_row!(
    csv_path::AbstractString,
    run::RunConfig;
    status::AbstractString,
    output_dir::AbstractString,
    posterior::AbstractDict = Dict{String, Float64}(),
    final_rmse = missing,
    labile_pin::Float64 = 0.0,
    scatter_stats = nothing,
)
    row = _build_row(run, status, output_dir, posterior, final_rmse, labile_pin)
    _fill_scatter_stats!(row, scatter_stats)
    _locked_append(csv_path, row)
    return nothing
end

# Map a forward_run `scatter_stats` NamedTuple onto the SCATTER_COLS cells. Each
# column "forward_<field>" reads NamedTuple field `<field>` (prefix stripped). A
# `nothing` (no forward run) leaves every scatter cell `missing`.
function _fill_scatter_stats!(row, scatter_stats)
    scatter_stats === nothing && return row
    for col in SCATTER_COLS
        field = Symbol(replace(col, r"^forward_" => ""))
        v = getfield(scatter_stats, field)
        # NaN (empty/constant series) is stored as missing for a clean CSV cell.
        row[col] = (v isa Real && isnan(v)) ? missing : v
    end
    return row
end

function _build_row(run, status, output_dir, posterior, final_rmse, labile_pin)
    calibrated = Set(calibrated_param_names(run))
    prior_lookup = Dict(run.priors)  # name => Prior

    row = Dict{String, Any}(
        "timestamp" => Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"),
        "status" => status,
        "run_identifier" => run.run_identifier,
        "restart_from" => run.restart_from === nothing ? missing : run.restart_from,
        "site" => run.site,
        "start" => string(run.start_date),
        "stop" => string(run.stop_date),
        "cal_depth" => run.cal_depth,
        "spinup_days" => run.spinup_days,
        "n_iterations" => run.n_iterations,
        "settingsdesc" => run.settingsdesc,
        "output_dir" => output_dir,
        "final_rmse" => final_rmse,
    )

    for p in PARAM_ORDER
        if p in calibrated
            pr = prior_lookup[p]
            row[prior_col(p, "mean")] = pr.mean
            row[prior_col(p, "std")] = pr.std
            row[prior_col(p, "lower")] = _toml_safe(pr.lower)
            row[prior_col(p, "upper")] = _toml_safe(pr.upper)
            row[post_col(p)] = get(posterior, p, missing)
        else
            for f in PRIOR_FIELDS
                row[prior_col(p, f)] = missing
            end
            # labile pinned to a known value when off; others unknown
            row[post_col(p)] = p == Config.LABILE_PARAM ? labile_pin : missing
        end
    end
    return row
end

# CSV can't write Inf cleanly into a numeric column shared with missing; store
# the string form for infinite bounds.
_toml_safe(x::Real) = isinf(x) ? (x > 0 ? "Inf" : "-Inf") : x

function _locked_append(csv_path, row::AbstractDict)
    mkpath(dirname(csv_path))
    cols = column_names()
    lock_path = csv_path * ".lock"
    _with_file_lock(lock_path) do
        write_header = !isfile(csv_path)
        # Assemble a 1-row DataFrame in fixed column order.
        df = DataFrame()
        for c in cols
            df[!, c] = [get(row, c, missing)]
        end
        CSV.write(csv_path, df; append = !write_header, writeheader = write_header)
    end
    return nothing
end

"""
A crude but portable cross-process file lock: create the lockfile exclusively,
spin (with backoff) until we win, run `f`, then remove it. Good enough for a
handful of sequential/parallel pipeline processes appending to one CSV.
"""
function _with_file_lock(f, lock_path; timeout_s = 120)
    t0 = time()
    fd = nothing
    while true
        try
            fd = Base.Filesystem.open(
                lock_path,
                Base.JL_O_CREAT | Base.JL_O_EXCL | Base.JL_O_WRONLY,
                0o644,
            )
            break
        catch
            time() - t0 > timeout_s &&
                error("Timed out acquiring CSV lock at $lock_path")
            sleep(0.05 + 0.1 * rand())
        end
    end
    try
        return f()
    finally
        fd === nothing || close(fd)
        isfile(lock_path) && rm(lock_path; force = true)
    end
end

# ── Restart lookup ───────────────────────────────────────────────────────────
"""
    read_posterior_means(final_params_file) -> Dict{String,Float64}

Parse a `final_parameter_means.txt` (lines like `  name = value`) into a Dict.
Only canonical parameter names are kept.
"""
function read_posterior_means(final_params_file::AbstractString)
    isfile(final_params_file) ||
        error("final_parameter_means.txt not found: $final_params_file")
    means = Dict{String, Float64}()
    for line in eachline(final_params_file)
        m = match(r"^\s*(\w+)\s*=\s*([-\d.eE+]+)\s*$", line)
        m === nothing && continue
        name = m.captures[1]
        name in PARAM_ORDER || continue
        means[name] = parse(Float64, m.captures[2])
    end
    isempty(means) &&
        error("No parameter means parsed from $final_params_file")
    return means
end

"""
    find_previous_run(csv_path, output_root, run_identifier;
                      site=nothing, start=nothing, stop=nothing) -> output_dir

Locate the output directory of a previous run by its `run_identifier`. Tries the
master CSV first (fast path), then falls back to scanning `output_root` for
`output_<id>` directories.

A single `run_identifier` is shared by every run in a batch, so on its own it is
NOT unique (many site/year rows). Pass `site` + `start` + `stop` (a run's
site/period) to narrow the match to that one calibration — this is what lets a
forward-only batch seed each site-year from its own posterior. After filtering:
  - one row → use it;
  - several rows all pointing at the SAME output_dir (re-runs that overwrote the
    same final_parameter_means.txt) → use that dir;
  - several rows with DIFFERENT output_dirs → take the newest by timestamp.
The unfiltered (id-only) path keeps the old strict behavior: error on >1 dir.
"""
function find_previous_run(
    csv_path::AbstractString,
    output_root::AbstractString,
    run_identifier::AbstractString;
    site = nothing,
    settingsdesc = nothing,
    start = nothing,
    stop = nothing,
)
    # Fast path: CSV lookup.
    if isfile(csv_path)
        df = CSV.read(csv_path, DataFrame)
        if "run_identifier" in names(df) && "output_dir" in names(df)
            mask = (string.(df.run_identifier) .== run_identifier) .&
                   (string.(df.status) .== "ok")
            # Narrow by site/period when given (compared as strings for robustness
            # against Date vs String parsing of the CSV columns).
            site === nothing || (mask = mask .& (string.(df.site) .== string(site)))
            settingsdesc === nothing ||
                (mask = mask .& (string.(df.settingsdesc) .== string(settingsdesc)))
            start === nothing || (mask = mask .& (string.(df.start) .== string(start)))
            stop === nothing || (mask = mask .& (string.(df.stop) .== string(stop)))
            sub = df[mask, :]

            if nrow(sub) >= 1
                dirs = unique(String.(sub.output_dir))
                length(dirs) == 1 && return dirs[1]
                # Multiple distinct dirs. With site/period given, take the newest
                # by timestamp; without, keep the strict ambiguity error.
                if site !== nothing || settingsdesc !== nothing ||
                   start !== nothing || stop !== nothing
                    if "timestamp" in names(sub)
                        order = sortperm(string.(sub.timestamp))   # ISO ts sorts lexically
                        return String(sub.output_dir[last(order)])
                    end
                    return String(sub.output_dir[end])
                end
                error("Multiple completed runs with identifier \"$run_identifier\" " *
                      "in $csv_path; pass site/start/stop to disambiguate.")
            end
        end
    end

    # Fallback: filesystem scan for output_<id>.
    target = "output_$run_identifier"
    matches = String[]
    for (root, dirs, _) in walkdir(output_root)
        for d in dirs
            d == target && push!(matches, joinpath(root, d))
        end
    end
    isempty(matches) &&
        error("No run found with identifier \"$run_identifier\" under $output_root")
    length(matches) > 1 &&
        error("Multiple dirs named $target under $output_root: $matches")
    return matches[1]
end

end # module
