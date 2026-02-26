"""
Pre-generate observation vector for DK-Sor calibration.

Reads daily FluxNet observations (NEE, Qle, Qh) for 2004-2013, creates one
EKP.Observation per year, and saves to JLD2 for use with ObservationSeries.

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

# Calibration years (after 1-year spinup starting 2003)
cal_years = 2004:2013

# Read all flux observations
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

# Global noise variances (computed across all calibration years)
all_cal_mask =
    (Date(first(cal_years), 1, 1) .<= flux_dates .<= Date(last(cal_years), 12, 31))
nee_var = mean(filter(!isnan, nee_uc_raw[all_cal_mask] .^ 2))
qle_var = mean(filter(!isnan, qle_uc_raw[all_cal_mask] .^ 2))
qh_var = mean(filter(!isnan, qh_uc_raw[all_cal_mask] .^ 2))
println(
    "Noise variances - NEE: $(round(nee_var, sigdigits=3)) (gC/m²/d)², " *
    "Qle: $(round(qle_var, sigdigits=3)) (W/m²)², " *
    "Qh: $(round(qh_var, sigdigits=3)) (W/m²)²",
)

# Build one EKP.Observation per year + save per-year dates
observation_vector = EKP.Observation[]
year_dates_dict = Dict{Int, Vector{Date}}()

for yr in cal_years
    yr_start = Date(yr, 1, 1)
    yr_end = Date(yr, 12, 31)
    in_year = yr_start .<= flux_dates .<= yr_end
    valid_mask =
        in_year .& .!isnan.(nee_raw) .& .!isnan.(qle_raw) .&
        .!isnan.(qh_raw) .& (abs.(nee_raw) .< 1e10) .&
        (abs.(qle_raw) .< 1e10) .& (abs.(qh_raw) .< 1e10)

    dates_yr = flux_dates[valid_mask]
    n_yr = length(dates_yr)

    y_obs_yr = vcat(nee_raw[valid_mask], qle_raw[valid_mask], qh_raw[valid_mask])
    noise_diag_yr = vcat(
        fill(nee_var, n_yr),
        fill(qle_var, n_yr),
        fill(qh_var, n_yr),
    )

    obs = EKP.Observation(
        Dict(
            "samples" => y_obs_yr,
            "covariances" => Diagonal(noise_diag_yr),
            "names" => "dk_sor_$yr",
        ),
    )
    push!(observation_vector, obs)
    year_dates_dict[yr] = dates_yr

    println("  $yr: $n_yr valid days, obs vector length $(length(y_obs_yr))")
end

println("\nTotal observations: $(length(observation_vector)) years ($(first(cal_years))-$(last(cal_years)))")

# Save
obs_filepath =
    joinpath(climaland_dir, "experiments/calibrate_dk_sor/observations.jld2")
JLD2.jldsave(
    obs_filepath;
    observation_vector = observation_vector,
    year_dates = year_dates_dict,
    cal_years = collect(cal_years),
)
println("Saved to $obs_filepath")
