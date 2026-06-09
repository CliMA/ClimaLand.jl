"""
Generate minibatched observation windows for DK-Sor CalLMIP Phase 1a calibration.

Builds 16 yearly windows (1997–2012) from the CalLMIP FLUXNET flux obs file.
Missing / gap-filled days are padded with zero value and very large variance,
so the EKP observation noise covariance correctly down-weights them.

Also writes the full-record (1997–2012) aggregated observation vector used
by the CES stage (fixed window, no minibatching).

Run once before calibration:
    julia --project=experiments/callmip_phase1a_v2 \\
          experiments/callmip_phase1a_v2/generate_observations.jl

Output: experiments/callmip_phase1a_v2/observations.jld2
"""

using NCDatasets
using Dates
using Statistics
using LinearAlgebra
import JLD2
import EnsembleKalmanProcesses as EKP

const FT = Float64
const climaland_dir = abspath(joinpath(@__DIR__, "..", ".."))

flux_nc_path = joinpath(
    climaland_dir, "DK_Sor",
    "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc",
)

# ── Yearly calibration windows: 1997–2012 (16 windows) ───────────────────────
# Year 2013 is excluded from calibration (last available obs year but kept as
# internal check by ME-org).  Year 2014 is the held-out temporal validation year
# (no obs available to us).
const CAL_YEARS   = collect(1997:2012)
const MINIBATCH_SIZE = 2

# ── Read flux observations ────────────────────────────────────────────────────
flux_ds      = NCDataset(flux_nc_path, "r")
flux_times   = flux_ds["time"][:]

nee_raw  = Float64.(coalesce.(flux_ds["NEE_daily"][:],    NaN))
qle_raw  = Float64.(coalesce.(flux_ds["Qle_daily"][:],    NaN))
qh_raw   = Float64.(coalesce.(flux_ds["Qh_daily"][:],     NaN))
nee_uc   = Float64.(coalesce.(flux_ds["NEE_uc_daily"][:], NaN))
qle_uc   = Float64.(coalesce.(flux_ds["Qle_uc_daily"][:], NaN))
qh_uc    = Float64.(coalesce.(flux_ds["Qh_uc_daily"][:],  NaN))
close(flux_ds)

flux_dates  = Date.(flux_times)
date_to_idx = Dict{Date, Int}(d => i for (i, d) in enumerate(flux_dates))

# Valid observations mask (not NaN, physically plausible)
valid_mask = (
    .!isnan.(nee_raw) .& .!isnan.(qle_raw) .& .!isnan.(qh_raw) .&
    (abs.(nee_raw) .< 1e6) .& (abs.(qle_raw) .< 1e6) .& (abs.(qh_raw) .< 1e6)
)

# Typical uncertainty variances from valid entries
let v = filter(!isnan, nee_uc[valid_mask] .^ 2); nee_var_typ = isempty(v) ? 1.0   : mean(v); end
let v = filter(!isnan, qle_uc[valid_mask] .^ 2); qle_var_typ = isempty(v) ? 100.0 : mean(v); end
let v = filter(!isnan, qh_uc[valid_mask]  .^ 2); qh_var_typ  = isempty(v) ? 100.0 : mean(v); end

const BIG_VAR = 1e12    # large variance for missing / gap-filled days

# ── Build per-year observation windows ───────────────────────────────────────
observation_vector = EKP.Observation[]
window_years  = Int[]
window_dates  = Vector{Vector{Date}}()

# Aggregated arrays for CES full-record fixed window
obs_dates_all = Date[]
nee_all       = FT[]
qle_all       = FT[]
qh_all        = FT[]
nee_var_all   = FT[]
qle_var_all   = FT[]
qh_var_all    = FT[]

println("Building per-year observation windows ($(length(CAL_YEARS)) years, 1997-2012):")

for yr in CAL_YEARS
    # Daily dates for this year, excluding Feb 29 to keep all windows equal length
    dates_full = collect(Date(yr, 1, 1):Day(1):Date(yr, 12, 31))
    dates = [d for d in dates_full if !(month(d) == 2 && day(d) == 29)]
    n     = length(dates)

    nee  = zeros(FT, n)
    qle  = zeros(FT, n)
    qh   = zeros(FT, n)
    nee_v = fill(FT(BIG_VAR), n)
    qle_v = fill(FT(BIG_VAR), n)
    qh_v  = fill(FT(BIG_VAR), n)

    n_valid = 0
    for (k, d) in enumerate(dates)
        idx = get(date_to_idx, d, 0)
        (idx == 0 || !valid_mask[idx]) && continue
        n_valid += 1
        nee[k]  = nee_raw[idx]
        qle[k]  = qle_raw[idx]
        qh[k]   = qh_raw[idx]
        nee_v[k] = isnan(nee_uc[idx]) ? nee_var_typ : nee_uc[idx]^2
        qle_v[k] = isnan(qle_uc[idx]) ? qle_var_typ : qle_uc[idx]^2
        qh_v[k]  = isnan(qh_uc[idx])  ? qh_var_typ  : qh_uc[idx]^2
    end

    y         = vcat(nee, qle, qh)
    noise_diag = vcat(nee_v, qle_v, qh_v)
    wname     = string(yr)

    push!(window_years, yr)
    push!(window_dates, dates)
    println("  Year $yr: $n_valid / $n valid days")

    append!(obs_dates_all, dates)
    append!(nee_all,  nee);  append!(qle_all, qle);  append!(qh_all, qh)
    append!(nee_var_all, nee_v); append!(qle_var_all, qle_v); append!(qh_var_all, qh_v)

    push!(
        observation_vector,
        EKP.Observation(
            Dict(
                "samples"     => y,
                "covariances" => Diagonal(noise_diag),
                "names"       => wname,
            ),
        ),
    )
end

uniform_days = length(window_dates[1])   # 365 (leap days dropped)
println("\nTotal windows: $(length(CAL_YEARS)),  days/window: $uniform_days")
println("Observation vector length per window: $(3 * uniform_days)")

# ── Full-record aggregated obs/cov for CES GP training ────────────────────────
y_obs_full   = vcat(nee_all, qle_all, qh_all)
noise_full   = Diagonal(vcat(nee_var_all, qle_var_all, qh_var_all))

# ── Save ──────────────────────────────────────────────────────────────────────
obs_filepath = joinpath(@__DIR__, "observations.jld2")
JLD2.jldsave(
    obs_filepath;
    observation_vector = observation_vector,
    window_years       = window_years,
    window_dates       = window_dates,
    y_obs_full         = y_obs_full,
    noise_full         = noise_full,
    obs_dates_full     = obs_dates_all,
)
println("\nSaved: $obs_filepath")
println("  Windows: $(length(observation_vector))")
println("  Full-record obs length: $(length(y_obs_full)) (3 × $(length(obs_dates_all)) days)")
