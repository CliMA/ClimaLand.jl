"""
Pre-generate observation vector for NEON site soil CO₂ calibration.

Reads NEON soil CO₂ concentration observations at the requested depth codes
(501/502/503), computes daily means across plots (001–005) per depth, and stacks
them into a single EKP.Observation (per-depth concatenation). Saves to JLD2 for use
by the calibration driver.

Multi-depth layout (option (a), per-depth concatenation):
    y_obs        = [y_<code1>; y_<code2>; ...]   (each depth keeps its own valid days)
    noise_cov    = Diagonal([noise_<code1>; noise_<code2>; ...])  (block-diagonal)
The per-depth date lists and model layers are saved so observation_map can rebuild
the same stacking order. Depths and the matching model layers come from the shared
lookup (neon_depth_lookup.jl → neon_depth_layer_mapping.csv).

Configuration via environment variables:
    NEON_SITE_ID     — NEON site ID (default: "NEON-srer")
    NEON_START_DATE  — Start date (default: site metadata start date, format: YYYY-mm-dd)
    NEON_SPINUP_DAYS — Number of spinup days (default: 20)
    CALL_DEPTHS      — comma-separated depth codes (default: "501,502,503")

Run once before calibration:
    julia --project=.buildkite experiments/calibrate_neon/generate_observations.jl
"""

using Dates
using Statistics
using LinearAlgebra
import JLD2
import EnsembleKalmanProcesses as EKP
import ClimaLand
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using CSV
using DataFrames

const FT = Float64
climaland_dir = pkgdir(ClimaLand)

# Shared site-metadata + depth-lookup helpers (order matters: lookup uses _neon_site_key).
include(joinpath(climaland_dir, "experiments/calibrate_neon/site_metadata.jl"))
include(joinpath(climaland_dir, "experiments/calibrate_neon/neon_depth_lookup.jl"))
#=
ENV["NEON_SITE_ID"]    = "NEON-cper"
ENV["NEON_START_DATE"] = "2017-01-01"   # Format: YYYY-mm-dd
ENV["NEON_STOP_DATE"]  = "2017-12-31"   # Format: YYYY-mm-dd
ENV["NEON_SPINUP_DAYS"] = "20"        # optional
=#

# ── Configuration ────────────────────────────────────────────────────────────
# Define SITE_ID / SPINUP_DAYS only if not already set in the session, so this
# script is self-sufficient whether run standalone, via set_Station.jl, or
# re-included in a loop (see run_priormean_pipeline_mult.jl). Plain globals
# (not const) so they stay reassignable across runs.
if !@isdefined(SITE_ID)
    SITE_ID = get(ENV, "NEON_SITE_ID", "NEON-srer")
end
if !@isdefined(SPINUP_DAYS)
    SPINUP_DAYS = parse(Int, get(ENV, "NEON_SPINUP_DAYS", "20"))
end
outputpath = get(ENV, "CALL_OUTPUT_PATH", "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_siteruns/$(SITE_ID)/")
site_ID_val = FluxnetSimulations.replace_hyphen(SITE_ID)

# Explicit depth-code list to calibrate on (per-depth concatenation order).
depth_codes = String.(split(get(ENV, "CALL_DEPTHS", "501,502,503"), ","))
depth_lookup = neon_depths_for_site(SITE_ID; codes = depth_codes)
println("Depth codes: ", join(depth_codes, ", "),
    "  (model layers: ", join([string(d.model_layer) for d in depth_lookup], ", "), ")")

time_offset = 0
(site_start_date, site_stop_date) =
    FluxnetSimulations.get_data_dates(SITE_ID, time_offset)
start_date =
    DateTime(get(ENV, "NEON_START_DATE", string(Date(site_start_date))))
stop_date = DateTime(get(ENV, "NEON_STOP_DATE", string(Date(site_stop_date))))

spinup_date = start_date + Day(SPINUP_DAYS)

println("Site: $SITE_ID")
println("Data period: $start_date to $stop_date")
println("Spinup until: $spinup_date")

# ── Load NEON CSV observations ──────────────────────────────────────────────
csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(SITE_ID)
println("Loading NEON data from $csv_path")
obs_df = CSV.read(csv_path, DataFrame)

# Parse timestamps → DateTime once (shared across depths)
obs_df[!, :datetime] =
    DateTime.(string.(Int.(obs_df.timestamp_fmt)), dateformat"yyyymmddHHMM")
obs_df[!, :date] = Date.(obs_df.datetime)

# Row-wise mean across plots (skip missing/NaN)
function rowmean_skipinvalid(row, cols)
    vals = Float64[]
    for c in cols
        v = row[c]
        if !ismissing(v) && !isnan(Float64(v))
            push!(vals, Float64(v))
        end
    end
    return isempty(vals) ? NaN : mean(vals)
end

# Inter-sensor variance per row (for the noise estimate)
function rowvar_skipinvalid(row, cols)
    vals = Float64[]
    for c in cols
        v = row[c]
        if !ismissing(v) && !isnan(Float64(v))
            push!(vals, Float64(v))
        end
    end
    return length(vals) >= 2 ? var(vals) : NaN
end

"""
    daily_series_for_code(code)

Daily mean soil CO₂ (across plots 001–005) and the per-day mean inter-sensor
variance for one depth `code`, trimmed to (spinup, stop] and NaN-day-free.
Returns the filtered/sorted per-depth DataFrame (columns :date, :daily_mean, :daily_var).
"""
function daily_series_for_code(code)
    co2_cols = [Symbol("soilCO2concentrationMean_$(lpad(p,3,'0'))_$code") for p in 1:5]
    mean_sym = Symbol("sco2_mean_$code")
    var_sym = Symbol("sco2_var_$code")
    obs_df[!, mean_sym] = [rowmean_skipinvalid(row, co2_cols) for row in eachrow(obs_df)]
    obs_df[!, var_sym] = [rowvar_skipinvalid(row, co2_cols) for row in eachrow(obs_df)]

    daily = combine(
        groupby(obs_df, :date),
        mean_sym =>
            (x -> begin
                valid = filter(!isnan, x)
                length(valid) >= 24 ? mean(valid) : NaN
            end) => :daily_mean,
        var_sym =>
            (x -> begin
                valid = filter(!isnan, x)
                isempty(valid) ? NaN : mean(valid)
            end) => :daily_var,
    )
    daily = filter(row -> row.date >= Date(spinup_date), daily)
    daily = filter(row -> row.date <= Date(stop_date), daily)
    daily = filter(row -> !isnan(row.daily_mean), daily)
    sort!(daily, :date)
    return daily
end

# ── Build the stacked (per-depth concatenated) observation ───────────────────
y_blocks = Vector{Float64}[]
noise_blocks = Vector{Float64}[]
obs_dates_per_code = Dict{String,Vector{Date}}()
n_obs_per_code = Int[]

for d in depth_lookup
    code = d.code
    daily = daily_series_for_code(code)
    n = nrow(daily)
    n == 0 && error("No valid observation days for site $SITE_ID depth code $code")

    # per-depth noise: mean inter-sensor variance broadcast over that depth's days
    mean_sensor_var = mean(filter(!isnan, daily.daily_var))

    push!(y_blocks, Float64.(daily.daily_mean))
    push!(noise_blocks, fill(mean_sensor_var, n))
    obs_dates_per_code[code] = daily.date
    push!(n_obs_per_code, n)
    println("  code $code (layer $(d.model_layer), z=$(d.z_obs_m) m): ",
        "$n days, noise=$(round(mean_sensor_var, sigdigits=3)) ppm²")
end

y_obs = reduce(vcat, y_blocks)
noise_diag = reduce(vcat, noise_blocks)
n_obs = length(y_obs)
println("Stacked observation vector length: $n_obs ",
    "(per-depth: ", join(n_obs_per_code, "+"), ")")

observation = EKP.Observation(
    Dict(
        "samples" => y_obs,
        "covariances" => Diagonal(noise_diag),
        "names" => "neon_sco2_$(SITE_ID)",
    ),
)

# ── Save ────────────────────────────────────────────────────────────────────
mkpath(outputpath)
obs_filepath = joinpath(outputpath, "observations.jld2")
JLD2.jldsave(
    obs_filepath;
    observation = observation,
    y_obs = y_obs,
    noise_cov = Diagonal(noise_diag),
    # multi-depth metadata — observation_map rebuilds the same stacking order from these
    depth_codes = [d.code for d in depth_lookup],
    model_layers = [d.model_layer for d in depth_lookup],
    z_obs_m = [d.z_obs_m for d in depth_lookup],
    obs_dates_per_code = obs_dates_per_code,
    n_obs_per_code = n_obs_per_code,
    # back-compat: flat concatenated date list in stacking order
    obs_dates = reduce(vcat, [obs_dates_per_code[d.code] for d in depth_lookup]),
    site_id = SITE_ID,
    spinup_days = SPINUP_DAYS,
)
println("Saved to $obs_filepath")
