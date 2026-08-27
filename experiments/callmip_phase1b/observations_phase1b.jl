"""
Build the multisite EKP observations for CalLMIP Phase 1b, Scenario 1.

For each calibration year and each of the 21 sites, monthly means of daily
NEE/Qle/Qh (site-major layout: [site1 NEE×12 Qle×12 Qh×12, site2 …], 756
entries per year). The recipe per site is exactly Phase 1a
(experiments/callmip_dksor/generate_observations.jl): months with fewer than
MIN_VALID_DAYS valid days get obs = 0 and noise = SIGMA2_MISS; noise for valid
months = FLUXNET random uncertainty (σ²/n) + the inter-annual
representativeness floor Var_y[monthly mean]. Sites with no flux coverage in a
year contribute 36 fully-masked entries (the forward model still simulates
them; the misfit ignores them).

Outputs (CALLMIP1B_OUTDIR, default experiments/callmip_phase1b/output_calibration):
  observations.jld2: calib_years, site_ids, obs_series (Vector{EKP.Observation},
  one per year), valid_dates_series (year → site → flux → Vector{Date}),
  coverage (site × year × flux valid-month counts).

Run: julia --project=.buildkite experiments/callmip_phase1b/observations_phase1b.jl
"""

import ClimaLand
import EnsembleKalmanProcesses as EKP
using NCDatasets
using Dates
using LinearAlgebra
using JLD2
using Statistics

include(joinpath(@__DIR__, "sites.jl"))

const CALIB_YEARS =
    haskey(ENV, "CALLMIP1B_CALIB_YEARS") ?
    parse.(Int, split(ENV["CALLMIP1B_CALIB_YEARS"], ",")) : collect(2003:2014)
const OUTDIR =
    get(ENV, "CALLMIP1B_OUTDIR", joinpath(@__DIR__, "output_calibration"))
const SIGMA2_MISS = 1.0e12
const MIN_VALID_DAYS = 5
const FLUX_NAMES = ("nee", "lhf", "shf")
const N_PER_SITE_YEAR = 36

"""
    site_monthly_blocks(flux_nc_path, calib_years)

One site's per-year (obs36, noise36, valid-dates) using the Phase 1a recipe.
Years outside the file's coverage are fully masked.
"""
function site_monthly_blocks(flux_nc_path, calib_years)
    NCDataset(flux_nc_path, "r") do ds
        dates = DateTime.(ds["time"][:])
        series = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
        for (fname, v, u) in (
            ("nee", "NEE_daily", "NEE_uc_daily"),
            ("lhf", "Qle_daily", "Qle_uc_daily"),
            ("shf", "Qh_daily", "Qh_uc_daily"),
        )
            vals = Float64.(coalesce.(ds[v][:], NaN))
            ucs = Float64.(coalesce.(ds[u][:], NaN))
            for arr in (vals, ucs)
                arr[arr .>= 1.0e19] .= NaN
            end
            series[fname] = (vec(vals), vec(ucs))
        end
        yrs_all = Dates.year.(dates)
        mons_all = Dates.month.(dates)

        # inter-annual representativeness floor per flux+month
        interannual = Dict(f => zeros(12) for f in FLUX_NAMES)
        for f in FLUX_NAMES
            vals, ucs = series[f]
            for mon in 1:12
                ymeans = Float64[]
                for yr in calib_years
                    sel = (yrs_all .== yr) .& (mons_all .== mon)
                    v = vals[sel]
                    u = ucs[sel]
                    ok = isfinite.(v) .& isfinite.(u) .& (u .> 0.0)
                    sum(ok) >= MIN_VALID_DAYS && push!(ymeans, mean(v[ok]))
                end
                interannual[f][mon] = length(ymeans) >= 2 ? var(ymeans) : 0.0
            end
        end

        out = Dict{Int, Any}()
        for yr in calib_years
            obs = Float64[]
            noise = Float64[]
            vdays = Dict(f => Date[] for f in FLUX_NAMES)
            nvalid = Dict(f => 0 for f in FLUX_NAMES)
            for f in FLUX_NAMES
                vals, ucs = series[f]
                for mon in 1:12
                    sel = (yrs_all .== yr) .& (mons_all .== mon)
                    v = vals[sel]
                    u = ucs[sel]
                    mdays = Date.(dates[sel])
                    ok = isfinite.(v) .& isfinite.(u) .& (u .> 0.0)
                    if sum(ok) >= MIN_VALID_DAYS
                        push!(obs, mean(v[ok]))
                        push!(
                            noise,
                            mean(u[ok] .^ 2) / sum(ok) + interannual[f][mon],
                        )
                        append!(vdays[f], mdays[ok])
                        nvalid[f] += 1
                    else
                        push!(obs, 0.0)
                        push!(noise, SIGMA2_MISS)
                    end
                end
            end
            @assert length(obs) == N_PER_SITE_YEAR
            out[yr] = (; obs, noise, vdays, nvalid)
        end
        return out
    end
end

# ── Main ─────────────────────────────────────────────────────────────────────
mkpath(OUTDIR)
site_ids = CALIBRATION_SITE_IDS
@info "Building Phase 1b obs: $(length(site_ids)) sites × years $(CALIB_YEARS[1])–$(CALIB_YEARS[end])"

blocks = Dict(
    id => site_monthly_blocks(
        ClimaLand.Artifacts.callmip_phase1_flux_path(id; phase = "1b"),
        CALIB_YEARS,
    ) for id in site_ids
)

obs_series = EKP.Observation[]
valid_dates_series = Dict{String, Dict{String, Vector{Date}}}[]
coverage = Dict{String, Dict{Int, Dict{String, Int}}}(
    id => Dict{Int, Dict{String, Int}}() for id in site_ids
)
for yr in CALIB_YEARS
    obs_vec = Float64[]
    noise_vec = Float64[]
    vd = Dict{String, Dict{String, Vector{Date}}}()
    for id in site_ids
        b = blocks[id][yr]
        append!(obs_vec, b.obs)
        append!(noise_vec, b.noise)
        vd[id] = b.vdays
        coverage[id][yr] = Dict(b.nvalid)
    end
    @assert length(obs_vec) == length(site_ids) * N_PER_SITE_YEAR
    push!(
        obs_series,
        EKP.Observation(obs_vec, Diagonal(noise_vec), "Phase1b_$(yr)"),
    )
    push!(valid_dates_series, vd)
end

n_unmasked =
    sum(o -> count(<(SIGMA2_MISS), diag(EKP.get_obs_noise_cov(o))), obs_series)
n_total = sum(o -> length(EKP.get_obs(o)), obs_series)
@info "Observation windows built" n_windows = length(obs_series) entries_per_window =
    length(site_ids) * N_PER_SITE_YEAR unmasked = n_unmasked total = n_total

jldsave(
    joinpath(OUTDIR, "observations.jld2");
    calib_years = CALIB_YEARS,
    site_ids,
    obs_series,
    valid_dates_series,
    coverage,
)
@info "Saved $(joinpath(OUTDIR, "observations.jld2"))"

# per-site coverage report (valid months per flux over all years)
println(
    "\nsite, nee_months, lhf_months, shf_months  (of $(12 * length(CALIB_YEARS)))",
)
for id in site_ids
    t = Dict(
        f => sum(get(coverage[id][yr], f, 0) for yr in CALIB_YEARS) for
        f in FLUX_NAMES
    )
    println("$id, $(t["nee"]), $(t["lhf"]), $(t["shf"])")
end
