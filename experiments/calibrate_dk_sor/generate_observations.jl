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

# ── Fixed-length observations (365 days per year) ───────────────────────────
# EKP's minibatching requires all Observations to have the same dimension.
# Each year may have different missing data, so we use a fixed 365-day calendar
# (Jan 1 – Dec 31, excluding Feb 29). Missing days get fill value 0 with very
# large noise variance so EKP effectively ignores them.
N_DAYS = 365  # fixed calendar (no Feb 29)
LARGE_VAR = 1e12  # noise variance for missing days — makes them uninformative

# Build a reference calendar (non-leap year dates as (month, day))
ref_calendar = [(month(d), day(d)) for d in Date(2001, 1, 1):Day(1):Date(2001, 12, 31)]
@assert length(ref_calendar) == N_DAYS

# Build per-year lookup: flux_date → index in raw arrays
flux_date_to_idx = Dict(flux_dates[i] => i for i in eachindex(flux_dates))

# Build one EKP.Observation per year + save per-year dates
observation_vector = EKP.Observation[]
year_dates_dict = Dict{Int, Vector{Date}}()

for yr in cal_years
    # Map reference calendar to actual dates in this year
    yr_dates = [Date(yr, m, d) for (m, d) in ref_calendar if
                try Date(yr, m, d); true catch; false end]
    # For leap years, ref_calendar already excludes Feb 29, so length is 365

    nee_yr = zeros(N_DAYS)
    qle_yr = zeros(N_DAYS)
    qh_yr = zeros(N_DAYS)
    nee_noise = fill(LARGE_VAR, N_DAYS)
    qle_noise = fill(LARGE_VAR, N_DAYS)
    qh_noise = fill(LARGE_VAR, N_DAYS)

    n_valid = 0
    for (i, d) in enumerate(yr_dates)
        idx = get(flux_date_to_idx, d, nothing)
        if !isnothing(idx) &&
           !isnan(nee_raw[idx]) && !isnan(qle_raw[idx]) && !isnan(qh_raw[idx]) &&
           abs(nee_raw[idx]) < 1e10 && abs(qle_raw[idx]) < 1e10 && abs(qh_raw[idx]) < 1e10
            nee_yr[i] = nee_raw[idx]
            qle_yr[i] = qle_raw[idx]
            qh_yr[i] = qh_raw[idx]
            nee_noise[i] = nee_var
            qle_noise[i] = qle_var
            qh_noise[i] = qh_var
            n_valid += 1
        end
    end

    y_obs_yr = vcat(nee_yr, qle_yr, qh_yr)
    noise_diag_yr = vcat(nee_noise, qle_noise, qh_noise)

    obs = EKP.Observation(
        Dict(
            "samples" => y_obs_yr,
            "covariances" => Diagonal(noise_diag_yr),
            "names" => "dk_sor_$yr",
        ),
    )
    push!(observation_vector, obs)
    year_dates_dict[yr] = yr_dates

    println("  $yr: $n_valid valid days ($(N_DAYS - n_valid) padded), obs vector length $(length(y_obs_yr))")
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
