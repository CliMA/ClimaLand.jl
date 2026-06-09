"""
Analyze posterior ensemble: coverage statistics and uncertainty bands.

Computes per-day p05/p50/p95 uncertainty bands for NEE/LE/H and saves them
for use by write_callmip_netcdf.jl.

Usage:
    julia --project=experiments/callmip_phase1a_v2 \\
          experiments/callmip_phase1a_v2/analyze_posterior_ensemble.jl
"""

using Dates
using Statistics
import JLD2
import ClimaCalibrate

const climaland_dir  = abspath(joinpath(@__DIR__, "..", ".."))
const exp_dir        = @__DIR__
const ens_output_dir = joinpath(exp_dir, "output_posterior_ensemble")
const flux_nc_path   = joinpath(climaland_dir, "DK_Sor",
                                 "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")
const START_DATE = Date(1997, 1, 1)
const STOP_DATE  = Date(2014, 12, 31)
const N_DAYS     = 6574
const mol_CO2_to_gC_day = 12.0 * 86400.0

all_dates = collect(START_DATE:Day(1):STOP_DATE)

# ── Collect ensemble daily diagnostics ────────────────────────────────────────
function load_member(m)
    path = joinpath(
        ClimaCalibrate.path_to_ensemble_member(ens_output_dir, 0, m),
        "callmip_diagnostics.jld2",
    )
    isfile(path) || return nothing
    d = JLD2.load(path)
    return d["dates"], d["surface_data"]
end

n_valid_members = 0
nee_matrix = fill(NaN, N_DAYS, 0)   # (time × members)
qle_matrix = fill(NaN, N_DAYS, 0)
qh_matrix  = fill(NaN, N_DAYS, 0)

for m in 1:50
    result = load_member(m)
    isnothing(result) && continue
    dates, sd = result

    d2i = Dict(dt => i for (i, dt) in enumerate(all_dates))
    nee_vec = fill(NaN, N_DAYS)
    qle_vec = fill(NaN, N_DAYS)
    qh_vec  = fill(NaN, N_DAYS)
    for (k, dt) in enumerate(dates)
        haskey(d2i, dt) || continue
        j = d2i[dt]
        nee_v = get(sd, "nee", Float64[]); length(nee_v) >= k && (nee_vec[j] = nee_v[k] * mol_CO2_to_gC_day)
        qle_v = get(sd, "lhf", Float64[]); length(qle_v) >= k && (qle_vec[j] = qle_v[k])
        qh_v  = get(sd, "shf", Float64[]); length(qh_v)  >= k && (qh_vec[j]  = qh_v[k])
    end
    nee_matrix = hcat(nee_matrix, nee_vec)
    qle_matrix = hcat(qle_matrix, qle_vec)
    qh_matrix  = hcat(qh_matrix,  qh_vec)
    global n_valid_members += 1
end

println("Loaded $n_valid_members valid ensemble members")

function quantile_bands(mat)
    p05 = [nanpercentile(mat[i, :], 5)  for i in 1:N_DAYS]
    p50 = [nanpercentile(mat[i, :], 50) for i in 1:N_DAYS]
    p95 = [nanpercentile(mat[i, :], 95) for i in 1:N_DAYS]
    return p05, p50, p95
end

function nanpercentile(v, p)
    vv = filter(!isnan, v)
    isempty(vv) && return NaN
    quantile(vv, p/100)
end

nee_p05, nee_p50, nee_p95 = quantile_bands(nee_matrix)
qle_p05, qle_p50, qle_p95 = quantile_bands(qle_matrix)
qh_p05,  qh_p50,  qh_p95  = quantile_bands(qh_matrix)

# Coverage: fraction of obs days where obs falls in [p05, p95]
using NCDatasets
ds = NCDatasets.NCDataset(flux_nc_path, "r")
flux_dates = Date.(ds["time"][:])
nee_obs_raw = Float64.(coalesce.(ds["NEE_daily"][:], NaN))
qle_obs_raw = Float64.(coalesce.(ds["Qle_daily"][:], NaN))
close(ds)

d2i_flux = Dict(dt => i for (i, dt) in enumerate(flux_dates))
nee_obs = [haskey(d2i_flux, dt) ? nee_obs_raw[d2i_flux[dt]] : NaN for dt in all_dates]
qle_obs = [haskey(d2i_flux, dt) ? qle_obs_raw[d2i_flux[dt]] : NaN for dt in all_dates]

valid = .!isnan.(nee_obs) .& .!isnan.(nee_p05) .& .!isnan.(nee_p95)
cov_nee = 100 * mean(nee_p05[valid] .<= nee_obs[valid] .<= nee_p95[valid])
println("NEE 90% coverage: $(round(cov_nee; digits=1))%  (target ≥ 70%)")

valid_q = .!isnan.(qle_obs) .& .!isnan.(qle_p05) .& .!isnan.(qle_p95)
cov_qle = 100 * mean(qle_p05[valid_q] .<= qle_obs[valid_q] .<= qle_p95[valid_q])
println("Qle 90% coverage: $(round(cov_qle; digits=1))%")

# Save
out_path = joinpath(exp_dir, "output_posterior_analysis", "posterior_uncertainty_bands.jld2")
isdir(dirname(out_path)) || mkpath(dirname(out_path))
JLD2.jldsave(
    out_path;
    dates = all_dates,
    nee_p05, nee_p50, nee_p95,
    qle_p05, qle_p50, qle_p95,
    qh_p05,  qh_p50,  qh_p95,
    cov_nee, cov_qle,
)
println("Uncertainty bands saved → $out_path")
