"""
Pre-generate minibatched observation windows for DK-Sor calibration.

Builds fixed-size two-year windows (2003-2004, ..., 2011-2012), drops leap days
to keep all windows equal length, and pads missing/invalid observations with
zeros plus very large variance. This mirrors the minibatched ObservationSeries
workflow used in Alexis's calibration setup.

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
const climaland_dir = abspath(joinpath(@__DIR__, "..", ".."))

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

# Two-year fixed windows for minibatching. Keep sample count divisible by 2.
window_pairs = [
    (2003, 2004),
    (2005, 2006),
    (2007, 2008),
    (2009, 2010),
    (2011, 2012),
    (2013, 2014),
]

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

# ── Build validity mask ───────────────────────────────────────────────────────
WIND_THRESHOLD = 5.0  # m/s
valid_mask =
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

# Build a date index for fast lookup
date_to_idx = Dict{Date, Int}(d => i for (i, d) in enumerate(flux_dates))

# Typical variances from valid entries
nee_var_typ = mean(filter(!isnan, nee_uc_raw[valid_mask] .^ 2))
qle_var_typ = mean(filter(!isnan, qle_uc_raw[valid_mask] .^ 2))
qh_var_typ = mean(filter(!isnan, qh_uc_raw[valid_mask] .^ 2))

# Large variance for padded/missing days
BIG_VAR = 1e12

observation_vector = EKP.Observation[]
window_names = String[]
window_dates = Vector{Vector{Date}}()

println("Building per-window observations ($(length(window_pairs)) uniform windows, 2003-2014):")
for (y0, y1) in window_pairs
    dates_full = collect(Date(y0, 1, 1):Day(1):Date(y1, 12, 31))
    # Drop leap day to force identical window length (730 days)
    dates = [d for d in dates_full if !(month(d) == 2 && day(d) == 29)]

    nee = zeros(FT, length(dates))
    qle = zeros(FT, length(dates))
    qh = zeros(FT, length(dates))

    nee_var = fill(FT(BIG_VAR), length(dates))
    qle_var = fill(FT(BIG_VAR), length(dates))
    qh_var = fill(FT(BIG_VAR), length(dates))

    n_valid = 0
    for (k, d) in enumerate(dates)
        idx = get(date_to_idx, d, 0)
        idx == 0 && continue
        if valid_mask[idx]
            n_valid += 1
            nee[k] = nee_raw[idx]
            qle[k] = qle_raw[idx]
            qh[k] = qh_raw[idx]
            nee_var[k] = isnan(nee_uc_raw[idx]) ? nee_var_typ : nee_uc_raw[idx]^2
            qle_var[k] = isnan(qle_uc_raw[idx]) ? qle_var_typ : qle_uc_raw[idx]^2
            qh_var[k] = isnan(qh_uc_raw[idx]) ? qh_var_typ : qh_uc_raw[idx]^2
        end
    end

    y = vcat(nee, qle, qh)
    noise_diag = vcat(nee_var, qle_var, qh_var)

    wname = "$(y0)_$(y1)"
    push!(window_names, wname)
    push!(window_dates, dates)
    println("  Window $(length(window_names)) ($(y0)-01-01–$(y1)-12-31): $(n_valid)/$(length(dates)) valid days")

    push!(
        observation_vector,
        EKP.Observation(
            Dict(
                "samples" => y,
                "covariances" => Diagonal(noise_diag),
                "names" => wname,
            ),
        ),
    )
end

uniform_days = length(window_dates[1])
println("\nTotal windows: $(length(window_pairs)), uniform length: $(3 * uniform_days)")

# ── Save ──────────────────────────────────────────────────────────────────────
obs_filepath =
    joinpath(climaland_dir, "experiments/calibrate_dk_sor/observations.jld2")
JLD2.jldsave(
    obs_filepath;
    observation_vector = observation_vector,
    window_names = window_names,
    window_dates = window_dates,
    window_pairs = window_pairs,
)
println("Saved $obs_filepath  ($(length(observation_vector)) windows, length $(3 * uniform_days) each)")
