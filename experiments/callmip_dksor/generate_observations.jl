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

import ClimaLand
import EnsembleKalmanProcesses as EKP
using NCDatasets
using Dates
using LinearAlgebra
using JLD2
using Statistics

# ── Paths ──────────────────────────────────────────────────────────────────────
# The FLUXNET2015 flux obs are resolved from the `callmip_phase1` artifact (same
# accessor used by model_runner.jl), so the pipeline runs from a clean clone with
# no data committed to the repo. See src/Artifacts.jl and the companion
# CliMA/ClimaArtifacts PR for the artifact definition + cluster hosting.
const CLIMALAND_DIR = abspath(joinpath(@__DIR__, "..", ".."))
const FLUX_NC_PATH  = ClimaLand.Artifacts.callmip_phase1_flux_path("DK-Sor"; phase = "1a")
const OUTDIR        = get(ENV, "CALLMIP_OUTDIR", joinpath(CLIMALAND_DIR,
    "experiments/callmip_dksor/output_calibration"))

# ── Calibration period ─────────────────────────────────────────────────────────
# Qle and Qh are only available from 2003 onward → 11 years of joint obs.
# Override with CALLMIP_CALIB_YEARS="2004" (comma-separated) for single-year tests.
const CALIB_YEARS = haskey(ENV, "CALLMIP_CALIB_YEARS") ?
    parse.(Int, split(ENV["CALLMIP_CALIB_YEARS"], ",")) : collect(2003:2013)

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
    # Per year: the calendar dates of the valid obs days, per flux. Threaded into
    # G(θ) so the MODEL monthly mean is taken over the SAME days as the obs mean
    # (sparse/biased-sampling months otherwise inject spurious model-data mismatch).
    valid_dates_series = Dict{String, Vector{Date}}[]

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

        fluxes = [("nee", nee_all, nee_uc),
                  ("lhf", qle_all, qle_uc),
                  ("shf", qh_all,  qh_uc)]

        # ---- First pass: inter-annual (representativeness) variance ----
        # The FLUXNET random uncertainty in the NetCDF (NEE_uc/Qle_uc/Qh_uc, propagated
        # as σ²/n) is unrealistically small — it claims we know each month's mean to
        # ~0.05 gC/m²/d (NEE). That over-confidence makes the CES emulator/MCMC blow up
        # (coverage ~0.18). The honest tolerance also includes REPRESENTATIVENESS: one
        # month is just a noisy sample of "a typical July", and July NEE varies year to
        # year. We add that year-to-year variance as a floor:
        #     σ²_realistic(f,m) = mean_y[σ²_random] + Var_y[ monthly_mean ]
        # See plots/fig_obs_noise_derivation.png. interannual_var[fname][mon] holds the
        # Var_y term (0.0 if a month has <2 valid years → no floor added).
        interannual_var = Dict(f[1] => zeros(12) for f in fluxes)
        for (fname, vals, ucs) in fluxes
            for mon in 1:12
                year_means = Float64[]
                for yr in calib_years
                    mask  = (Dates.year.(dates) .== yr) .& (Dates.month.(dates) .== mon)
                    v_mon = vals[mask]
                    u_mon = ucs[mask]
                    valid = isfinite.(v_mon) .& isfinite.(u_mon) .& (u_mon .> 0.0)
                    if sum(valid) >= MIN_VALID_DAYS
                        push!(year_means, mean(v_mon[valid]))
                    end
                end
                interannual_var[fname][mon] = length(year_means) >= 2 ? var(year_means) : 0.0
            end
        end

        for yr in calib_years
            obs_vec   = Float64[]
            noise_vec = Float64[]
            vdays     = Dict("nee" => Date[], "lhf" => Date[], "shf" => Date[])

            for (fname, vals, ucs) in fluxes
                for mon in 1:12
                    mask  = (Dates.year.(dates) .== yr) .& (Dates.month.(dates) .== mon)
                    mdays = Date.(dates[mask])
                    v_mon = vals[mask]
                    u_mon = ucs[mask]

                    # Keep only days with finite value AND finite positive uncertainty
                    valid = isfinite.(v_mon) .& isfinite.(u_mon) .& (u_mon .> 0.0)
                    n_valid = sum(valid)

                    if n_valid >= MIN_VALID_DAYS
                        push!(obs_vec,   mean(v_mon[valid]))
                        # Realistic noise = FLUXNET random (σ²/n) + inter-annual floor
                        random_var = mean(u_mon[valid].^2) / n_valid
                        push!(noise_vec, random_var + interannual_var[fname][mon])
                        append!(vdays[fname], mdays[valid])   # record valid obs days
                    else
                        push!(obs_vec,   0.0)
                        push!(noise_vec, SIGMA2_MISS)
                        # month dropped → leave no valid dates (model falls back to all days)
                    end
                end
            end

            @assert length(obs_vec) == N_OBS_PER_YEAR
            push!(obs_series,
                EKP.Observation(obs_vec, Diagonal(noise_vec), "DK_Sor_$(yr)"))
            push!(valid_dates_series, vdays)
        end
    end

    return obs_series, valid_dates_series
end

# ── Main ───────────────────────────────────────────────────────────────────────
mkpath(OUTDIR)

@info "Building monthly obs series for years $(CALIB_YEARS[1])–$(CALIB_YEARS[end])…"
obs_series, valid_dates_series = build_obs_series(FLUX_NC_PATH, CALIB_YEARS)

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
jldsave(out_path; calib_years = CALIB_YEARS, obs_series, valid_dates_series)
@info "Saved → $out_path  (with per-flux valid-day dates for model-side matching)"
