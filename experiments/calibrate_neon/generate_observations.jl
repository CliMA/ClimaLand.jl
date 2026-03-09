"""
Pre-generate observation vector for NEON site soil CO₂ calibration.

Reads NEON soil CO₂ concentration observations at depth 502 (~6 cm),
computes daily means across plots (001–005), and creates a single
EKP.Observation. Saves to JLD2 for use by the calibration driver.

Configuration via environment variables:
    NEON_SITE_ID     — NEON site ID (default: "NEON-srer")
    NEON_SPINUP_DAYS — Number of spinup days (default: 20)

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
const climaland_dir = pkgdir(ClimaLand)

# ── Configuration ────────────────────────────────────────────────────────────
const SITE_ID = get(ENV, "NEON_SITE_ID", "NEON-srer")
const SPINUP_DAYS = parse(Int, get(ENV, "NEON_SPINUP_DAYS", "20"))

site_ID_val = FluxnetSimulations.replace_hyphen(SITE_ID)

# Get dates from site metadata
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(start_date, stop_date) =
    FluxnetSimulations.get_data_dates(SITE_ID, time_offset)
spinup_date = start_date + Day(SPINUP_DAYS)

println("Site: $SITE_ID")
println("Data period: $start_date to $stop_date")
println("Spinup until: $spinup_date")

# ── Load NEON CSV observations ──────────────────────────────────────────────
csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(SITE_ID)
obs_df = CSV.read(csv_path, DataFrame)

# Extract depth-502 columns for plots 001–005
co2_cols_502 = [
    Symbol("soilCO2concentrationMean_001_502"),
    Symbol("soilCO2concentrationMean_002_502"),
    Symbol("soilCO2concentrationMean_003_502"),
    Symbol("soilCO2concentrationMean_004_502"),
    Symbol("soilCO2concentrationMean_005_502"),
]

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

obs_df[!, :sco2_mean_502] =
    [rowmean_skipinvalid(row, co2_cols_502) for row in eachrow(obs_df)]

# Parse timestamps → DateTime
obs_df[!, :datetime] =
    DateTime.(string.(Int.(obs_df.timestamp_fmt)), dateformat"yyyymmddHHMM")

# Compute inter-sensor variance per row for noise estimate
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

obs_df[!, :sco2_var_502] =
    [rowvar_skipinvalid(row, co2_cols_502) for row in eachrow(obs_df)]

obs_df[!, :date] = Date.(obs_df.datetime)

# Group by date, compute daily means (require ≥24 valid half-hours)
daily_df = combine(
    groupby(obs_df, :date),
    :sco2_mean_502 =>
        (x -> begin
            valid = filter(!isnan, x)
            length(valid) >= 24 ? mean(valid) : NaN
        end) => :daily_mean,
    :sco2_var_502 =>
        (x -> begin
            valid = filter(!isnan, x)
            isempty(valid) ? NaN : mean(valid)
        end) => :daily_var,
)

# Trim to after spinup and remove NaN days
daily_df = filter(row -> row.date >= Date(spinup_date), daily_df)
daily_df = filter(row -> !isnan(row.daily_mean), daily_df)
sort!(daily_df, :date)

y_obs = Float64.(daily_df.daily_mean)
obs_dates = daily_df.date
n_obs = length(y_obs)
println("Valid observation days: $n_obs (after $SPINUP_DAYS-day spinup)")

# ── Noise covariance ────────────────────────────────────────────────────────
mean_sensor_var = mean(filter(!isnan, daily_df.daily_var))
noise_diag = fill(mean_sensor_var, n_obs)
println("Noise variance (inter-sensor): $(round(mean_sensor_var, sigdigits=3)) ppm²")

# ── Build observation ────────────────────────────────────────────────────────
observation = EKP.Observation(
    Dict(
        "samples" => y_obs,
        "covariances" => Diagonal(noise_diag),
        "names" => "neon_sco2_$(SITE_ID)",
    ),
)

println("Observation vector length: $n_obs")

# ── Save ────────────────────────────────────────────────────────────────────
obs_filepath =
    joinpath(climaland_dir, "experiments/calibrate_neon/observations.jld2")
JLD2.jldsave(
    obs_filepath;
    observation = observation,
    obs_dates = obs_dates,
    y_obs = y_obs,
    noise_cov = Diagonal(noise_diag),
    site_id = SITE_ID,
    spinup_days = SPINUP_DAYS,
)
println("Saved to $obs_filepath")
