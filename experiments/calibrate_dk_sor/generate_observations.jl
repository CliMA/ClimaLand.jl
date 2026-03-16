"""
Pre-generate observation vector for DK-Sor calibration.

Reads daily FluxNet observations (NEE, Qle, Qh) for 2004-2013, filters for
days with valid data and daily mean wind speed < 5 m/s, and creates a single
EKP.Observation. Saves to JLD2 for use by the calibration driver.

Run once before calibration:
    julia --project=.buildkite experiments/calibrate_dk_sor/generate_observations.jl
"""

using NCDatasets
using Dates
using Statistics
using LinearAlgebra
import JLD2
import EnsembleKalmanProcesses as EKP
import ClimaLand

const FT = Float64
const climaland_dir = pkgdir(ClimaLand)

flux_nc_path = joinpath(
    climaland_dir,
    "DK_Sor",
    "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc",
)
met_nc_path = joinpath(
    climaland_dir,
    "DK_Sor",
    "DK-Sor_1997-2014_FLUXNET2015_Met.nc",
)

# Calibration years (after 1-year spinup starting 2003)
cal_years = 2004:2013

# ── Read flux observations ───────────────────────────────────────────────────
flux_ds = NCDataset(flux_nc_path, "r")
flux_times_dt = flux_ds["time"][:]

nee_raw = Float64.(coalesce.(flux_ds["NEE_daily"][:], NaN))
qle_raw = Float64.(coalesce.(flux_ds["Qle_daily"][:], NaN))
qh_raw = Float64.(coalesce.(flux_ds["Qh_daily"][:], NaN))

nee_uc_raw = Float64.(coalesce.(flux_ds["NEE_uc_daily"][:], NaN))
qle_uc_raw = Float64.(coalesce.(flux_ds["Qle_uc_daily"][:], NaN))
qh_uc_raw = Float64.(coalesce.(flux_ds["Qh_uc_daily"][:], NaN))
close(flux_ds)

flux_dates = Date.(flux_times_dt)

# ── Compute daily mean wind speed from hourly Met data ────────────────────────
met_ds = NCDataset(met_nc_path, "r")
wind_data = Float64.(coalesce.(met_ds["Wind"][1, 1, :], NaN))
met_times = met_ds["time"][:]
close(met_ds)

met_dates = Date.(met_times)

# Average wind speed per day
daily_wind = Dict{Date, Float64}()
wind_by_day = Dict{Date, Vector{Float64}}()
for (i, d) in enumerate(met_dates)
    v = wind_data[i]
    isnan(v) && continue
    if !haskey(wind_by_day, d)
        wind_by_day[d] = Float64[]
    end
    push!(wind_by_day[d], v)
end
for (d, vals) in wind_by_day
    daily_wind[d] = mean(vals)
end

# ── Build valid-day mask ──────────────────────────────────────────────────────
WIND_THRESHOLD = 5.0  # m/s

cal_mask =
    Date(first(cal_years), 1, 1) .<= flux_dates .<= Date(last(cal_years), 12, 31)

valid_mask = cal_mask .&
    .!isnan.(nee_raw) .& .!isnan.(qle_raw) .& .!isnan.(qh_raw) .&
    (abs.(nee_raw) .< 1e10) .& (abs.(qle_raw) .< 1e10) .& (abs.(qh_raw) .< 1e10)

# Apply wind filter
for i in eachindex(valid_mask)
    if valid_mask[i]
        w = get(daily_wind, flux_dates[i], NaN)
        if isnan(w) || w >= WIND_THRESHOLD
            valid_mask[i] = false
        end
    end
end

obs_dates = flux_dates[valid_mask]
nee_obs = nee_raw[valid_mask]
qle_obs = qle_raw[valid_mask]
qh_obs = qh_raw[valid_mask]
nee_uc = nee_uc_raw[valid_mask]
qle_uc = qle_uc_raw[valid_mask]
qh_uc = qh_uc_raw[valid_mask]

n_obs = length(obs_dates)
println("Valid observation days: $n_obs ($(first(cal_years))-$(last(cal_years)), wind < $(WIND_THRESHOLD) m/s)")

# ── Noise covariance ──────────────────────────────────────────────────────────
nee_var = mean(filter(!isnan, nee_uc .^ 2))
qle_var = mean(filter(!isnan, qle_uc .^ 2))
qh_var = mean(filter(!isnan, qh_uc .^ 2))
println(
    "Noise variances - NEE: $(round(nee_var, sigdigits=3)) (gC/m²/d)², " *
    "Qle: $(round(qle_var, sigdigits=3)) (W/m²)², " *
    "Qh: $(round(qh_var, sigdigits=3)) (W/m²)²",
)

# ── Build single observation ─────────────────────────────────────────────────
y_obs = vcat(nee_obs, qle_obs, qh_obs)
noise_diag = vcat(
    fill(nee_var, n_obs),
    fill(qle_var, n_obs),
    fill(qh_var, n_obs),
)

observation = EKP.Observation(
    Dict(
        "samples" => y_obs,
        "covariances" => Diagonal(noise_diag),
        "names" => "dk_sor_$(first(cal_years))_$(last(cal_years))",
    ),
)

println("Observation vector length: $(length(y_obs)) (3 × $n_obs)")

# ── Save ──────────────────────────────────────────────────────────────────────
obs_filepath =
    joinpath(climaland_dir, "experiments/calibrate_dk_sor/observations.jld2")
JLD2.jldsave(
    obs_filepath;
    observation = observation,
    obs_dates = obs_dates,
    cal_years = collect(cal_years),
    y_obs = y_obs,
    noise_cov = Diagonal(noise_diag),
)
println("Saved to $obs_filepath")
