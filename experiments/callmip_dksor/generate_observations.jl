"""
Generate EKP observations for DK-Sor CalLMIP Phase 1a calibration.

Reads the FLUXNET2015 daily flux file (1997-2013) and builds 11 yearly
observation windows covering 2003-2013 (the period when NEE, Qle, AND Qh
are all available). Each window is a single EKP.Observation with:

  - Observation vector: [NEE_jan..dec, Qle_jan..dec, Qh_jan..dec]
    12 monthly means per variable = 36 entries per year (fixed, no leap issue)
  - NEE units: gC m⁻² d⁻¹  (obs × 12 × 86400 applied in run_calibration.jl)
  - Qle, Qh units: W m⁻²
  - Months with < MIN_VALID_DAYS valid days → obs=0, noise=SIGMA2_MISS

Outputs
-------
  output_calibration/observations.jld2  with keys:
    "calib_years"   — Vector{Int} of calibration years [2003..2013]
    "obs_series"    — Vector{EKP.Observation} (length 11, each 36-entry)

Usage
-----
  julia --project=.buildkite \
        experiments/callmip_dksor/generate_observations.jl
"""

import EnsembleKalmanProcesses as EKP
using NCDatasets
using Dates
using LinearAlgebra
using JLD2
using Statistics

# ── Paths ──────────────────────────────────────────────────────────────────────
const CLIMALAND_DIR = abspath(joinpath(@__DIR__, "..", ".."))
const FLUX_NC_PATH  = joinpath(CLIMALAND_DIR, "DK_Sor",
    "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")
const OUTDIR        = joinpath(CLIMALAND_DIR,
    "experiments/callmip_dksor/output_calibration")

# ── Calibration period ─────────────────────────────────────────────────────────
# Qle and Qh are only available from 2003 onward → 11 years of joint obs.
const CALIB_YEARS = collect(2003:2013)

# 12 monthly means × 3 variables = 36 entries per year (always fixed)
const N_OBS_PER_YEAR = 36

# Large noise for months with insufficient valid data
const SIGMA2_MISS    = 1.0e12
const MIN_VALID_DAYS = 5   # require at least 5 valid days to trust a monthly mean

"""
    build_obs_series(flux_nc_path, calib_years) -> Vector{EKP.Observation}

For each year, compute monthly means of NEE, Qle, Qh from daily FLUXNET data.
Months with fewer than MIN_VALID_DAYS valid observations receive obs=0 and
noise=SIGMA2_MISS. Returns 36-entry obs vectors (always fixed length).
"""
function build_obs_series(flux_nc_path, calib_years)
    obs_series = EKP.Observation[]

    NCDataset(flux_nc_path, "r") do ds
        # NCDatasets decodes time axis directly to DateTime
        dates   = DateTime.(ds["time"][:])
        nee_all = Float64.(coalesce.(ds["NEE_daily"][:], NaN))
        qle_all = Float64.(coalesce.(ds["Qle_daily"][:], NaN))
        qh_all  = Float64.(coalesce.(ds["Qh_daily"][:],  NaN))
        nee_uc  = Float64.(coalesce.(ds["NEE_uc_daily"][:], NaN))
        qle_uc  = Float64.(coalesce.(ds["Qle_uc_daily"][:], NaN))
        qh_uc   = Float64.(coalesce.(ds["Qh_uc_daily"][:],  NaN))

        # Sentinel values (FillValue ≥ 1e19) → NaN
        for arr in (nee_all, qle_all, qh_all, nee_uc, qle_uc, qh_uc)
            arr[arr .>= 1.0e19] .= NaN
        end

        for yr in calib_years
            obs_vec   = Float64[]
            noise_vec = Float64[]

            for (vals, ucs) in [(nee_all, nee_uc),
                                (qle_all, qle_uc),
                                (qh_all,  qh_uc)]
                for mon in 1:12
                    mask  = (Dates.year.(dates) .== yr) .& (Dates.month.(dates) .== mon)
                    v_mon = vals[mask]
                    u_mon = ucs[mask]

                    # Keep only days with finite value AND finite positive uncertainty
                    valid = isfinite.(v_mon) .& isfinite.(u_mon) .& (u_mon .> 0.0)
                    n_valid = sum(valid)

                    if n_valid >= MIN_VALID_DAYS
                        push!(obs_vec,   mean(v_mon[valid]))
                        # Propagate uncertainty: σ²_mean = mean(σ²) / n
                        push!(noise_vec, mean(u_mon[valid].^2) / n_valid)
                    else
                        push!(obs_vec,   0.0)
                        push!(noise_vec, SIGMA2_MISS)
                    end
                end
            end

            @assert length(obs_vec) == N_OBS_PER_YEAR
            push!(obs_series,
                EKP.Observation(obs_vec, Diagonal(noise_vec), "DK_Sor_$(yr)"))
        end
    end

    return obs_series
end

# ── Main ───────────────────────────────────────────────────────────────────────
mkpath(OUTDIR)

@info "Building monthly obs series for years $(CALIB_YEARS[1])–$(CALIB_YEARS[end])…"
obs_series = build_obs_series(FLUX_NC_PATH, CALIB_YEARS)

# ── Validation ─────────────────────────────────────────────────────────────────
@assert length(obs_series) == length(CALIB_YEARS)
for (i, obs) in enumerate(obs_series)
    yr      = CALIB_YEARS[i]
    y_vec   = EKP.get_obs(obs)
    Γ_vec   = diag(EKP.get_obs_noise_cov(obs))
    n_valid = sum(Γ_vec .< SIGMA2_MISS)
    n_miss  = sum(Γ_vec .>= SIGMA2_MISS)
    @assert length(y_vec) == N_OBS_PER_YEAR
    @assert all(Γ_vec .> 0)
    nee_slice = y_vec[1:12]
    @info "  $yr: valid=$n_valid/36 | missing=$n_miss | " *
          "NEE∈[$(round(minimum(nee_slice); digits=2)), $(round(maximum(nee_slice); digits=2))] gC/m²/d"
end

# ── Save ───────────────────────────────────────────────────────────────────────
out_path = joinpath(OUTDIR, "observations.jld2")
jldsave(out_path; calib_years = CALIB_YEARS, obs_series)
@info "Saved → $out_path"
