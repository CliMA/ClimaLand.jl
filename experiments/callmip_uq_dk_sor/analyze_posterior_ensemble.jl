"""
Analysis and plotting for the posterior forward ensemble.

Reads `posterior_ensemble_diagnostics.jld2` produced by
`run_posterior_ensemble.jl` plus the FLUXNET observations, and:

  1. Computes 5th / 50th / 95th percentile uncertainty bands on
     daily NEE, Qle, and Qh.
  2. Plots ensemble spread vs FLUXNET observations for each flux.
  3. Computes coverage (fraction of observations within the 90 % band).
  4. Saves a structured JLD2 with the band statistics for CalLMIP submission.

Usage
-----
    julia --project=experiments/callmip_uq_dk_sor \\
          experiments/callmip_uq_dk_sor/analyze_posterior_ensemble.jl
"""

using JLD2, Dates, Statistics
using Plots, Printf
import NCDatasets
import ClimaLand

# ── Paths ─────────────────────────────────────────────────────────────────────
const climaland_dir = abspath(joinpath(@__DIR__, "..", ".."))
const cal_dir       = joinpath(climaland_dir, "experiments", "calibrate_dk_sor")  # calibration artifacts
const exp_dir       = joinpath(climaland_dir, "experiments", "callmip_uq_dk_sor")  # UQ outputs
const ensemble_dir  = joinpath(exp_dir, "output_posterior_ensemble")
const out_dir       = joinpath(exp_dir, "output_posterior_analysis")
isdir(out_dir) || mkdir(out_dir)

flux_nc_path = joinpath(
    climaland_dir,
    "DK_Sor",
    "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc",
)

# ── Load ensemble diagnostics ─────────────────────────────────────────────────
diag_file = joinpath(ensemble_dir, "posterior_ensemble_diagnostics.jld2")
isfile(diag_file) || error("Run run_posterior_ensemble.jl first: $diag_file not found")

data = JLD2.load(diag_file)
member_results          = data["member_results"]    # Dict{Int, NamedTuple}
posterior_samples       = data["posterior_samples"]
constrained_ekp_optimal = data["constrained_ekp_optimal"]
param_names             = data["param_names"]
n_samples               = data["n_samples"]

members = sort(collect(keys(member_results)))
@info "Loaded $(length(members)) member runs out of $n_samples requested"

# ── Load FLUXNET observations ──────────────────────────────────────────────────
obs_ds     = NCDatasets.NCDataset(flux_nc_path, "r")
nee_obs_raw = Float64.(NCDatasets.coalesce.(obs_ds["NEE_daily"][:], NaN))
qle_obs_raw = Float64.(NCDatasets.coalesce.(obs_ds["Qle_daily"][:], NaN))
qh_obs_raw  = Float64.(NCDatasets.coalesce.(obs_ds["Qh_daily"][:], NaN))
obs_times   = obs_ds["time"][:]
NCDatasets.close(obs_ds)

obs_dates   = Date.(obs_times)
obs_dict_nee = Dict(zip(obs_dates, nee_obs_raw))
obs_dict_qle = Dict(zip(obs_dates, qle_obs_raw))
obs_dict_qh  = Dict(zip(obs_dates, qh_obs_raw))

# ── Align ensemble on a common date vector ────────────────────────────────────
# Use the date range common to all successful members.
all_dates = sort(unique(reduce(vcat, [collect(member_results[m].dates) for m in members])))
n_dates   = length(all_dates)
@info "Common date range: $(first(all_dates)) – $(last(all_dates)) ($n_dates days)"

function build_matrix(flux_key)
    mat = fill(NaN, n_dates, length(members))
    for (j, m) in enumerate(members)
        r  = member_results[m]
        d  = Dict(zip(r.dates, getfield(r, flux_key)))
        for (i, date) in enumerate(all_dates)
            mat[i, j] = get(d, date, NaN)
        end
    end
    return mat
end

nee_mat = build_matrix(:nee)   # (n_dates × n_members), mol CO₂/m²/s
qle_mat = build_matrix(:qle)   # (n_dates × n_members), W/m²
qh_mat  = build_matrix(:qh)    # (n_dates × n_members), W/m²

# Convert NEE model output to gC/m²/d (same units as FLUXNET) for comparison
nee_mat_gc = nee_mat .* 12.0 .* 86400.0

# ── Percentile bands ─────────────────────────────────────────────────────────
function nanpercentile(v, p)
    vv = filter(!isnan, v)
    isempty(vv) && return NaN
    return quantile(vv, p / 100)
end

function compute_bands(mat)
    p05 = [nanpercentile(mat[i, :], 5)  for i in 1:n_dates]
    p50 = [nanpercentile(mat[i, :], 50) for i in 1:n_dates]
    p95 = [nanpercentile(mat[i, :], 95) for i in 1:n_dates]
    return p05, p50, p95
end

nee_p05, nee_p50, nee_p95 = compute_bands(nee_mat_gc)
qle_p05, qle_p50, qle_p95 = compute_bands(qle_mat)
qh_p05,  qh_p50,  qh_p95  = compute_bands(qh_mat)

# ── Coverage diagnostics ──────────────────────────────────────────────────────
function coverage_90pct(obs_dict, p05, p95)
    n_in = 0; n_total = 0
    for (i, d) in enumerate(all_dates)
        v = get(obs_dict, d, NaN)
        isnan(v) && continue
        n_total += 1
        if !isnan(p05[i]) && !isnan(p95[i]) && p05[i] <= v <= p95[i]
            n_in += 1
        end
    end
    return n_total > 0 ? n_in / n_total : NaN
end

cov_nee = coverage_90pct(obs_dict_nee, nee_p05, nee_p95)
cov_qle = coverage_90pct(obs_dict_qle, qle_p05, qle_p95)
cov_qh  = coverage_90pct(obs_dict_qh,  qh_p05,  qh_p95)

println("\n╔══════════════════════════════════════════╗")
println("║  90 % posterior predictive coverage      ║")
println("╠══════════════════════════════════════════╣")
@printf("║  NEE : %5.1f %%                           ║\n", 100 * cov_nee)
@printf("║  Qle : %5.1f %%                           ║\n", 100 * cov_qle)
println("╚══════════════════════════════════════════╝")
@printf("║  Qh  : %5.1f %%                           ║\n", 100 * cov_qh)

# ── Plots ─────────────────────────────────────────────────────────────────────
function flux_plot(dates, p05, p50, p95, obs_dict, ylabel, flux_name)
    obs_y = [get(obs_dict, d, NaN) for d in dates]
    t     = dates

    pl = plot(t, p50;
        ribbon = (p50 .- p05, p95 .- p50),
        fillalpha = 0.3, lw = 1.5, lc = :steelblue, fc = :steelblue,
        label = "posterior median (90 % band)",
        ylabel = ylabel, xlabel = "Date",
        title  = "DK-Sor: $flux_name posterior uncertainty",
        legend = :topleft, size = (1200, 400),
    )
    scatter!(pl, t, obs_y;
        ms = 2, mc = :black, ma = 0.5, label = "FLUXNET obs",
    )
    return pl
end

pl_nee = flux_plot(all_dates, nee_p05, nee_p50, nee_p95, obs_dict_nee,
                   "NEE (gC m⁻² d⁻¹)", "NEE")
pl_qle = flux_plot(all_dates, qle_p05, qle_p50, qle_p95, obs_dict_qle,
                   "Qle (W m⁻²)", "Latent heat flux")
pl_qh  = flux_plot(all_dates, qh_p05,  qh_p50,  qh_p95,  obs_dict_qh,
                   "Qh (W m⁻²)",  "Sensible heat flux")

savefig(pl_nee, joinpath(out_dir, "posterior_uncertainty_nee.png"))
savefig(pl_qle, joinpath(out_dir, "posterior_uncertainty_qle.png"))
savefig(pl_qh,  joinpath(out_dir, "posterior_uncertainty_qh.png"))

combined = plot(pl_nee, pl_qle, pl_qh; layout = (3, 1), size = (1200, 1000))
savefig(combined, joinpath(out_dir, "posterior_uncertainty_all_fluxes.png"))
@info "Plots saved to $out_dir"

# ── Save structured output for CalLMIP submission ─────────────────────────────
JLD2.jldsave(
    joinpath(out_dir, "posterior_uncertainty_bands.jld2");
    dates             = all_dates,
    nee_p05           = nee_p05,
    nee_p50           = nee_p50,
    nee_p95           = nee_p95,
    qle_p05           = qle_p05,
    qle_p50           = qle_p50,
    qle_p95           = qle_p95,
    qh_p05            = qh_p05,
    qh_p50            = qh_p50,
    qh_p95            = qh_p95,
    coverage_nee_90pct = cov_nee,
    coverage_qle_90pct = cov_qle,
    coverage_qh_90pct  = cov_qh,
    n_members         = length(members),
    param_names       = param_names,
    posterior_samples = posterior_samples,
)
@info "Uncertainty bands saved → $(joinpath(out_dir, "posterior_uncertainty_bands.jld2"))"
@info "Analysis complete."
