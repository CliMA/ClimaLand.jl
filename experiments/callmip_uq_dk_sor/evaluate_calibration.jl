"""
Calibration evaluation and visualisation for CalLMIP Phase 1a — DK-Sor.

Produces a comprehensive set of figures and a skill-score table comparing
the ClimaLand prior and posterior simulations against FLUXNET2015 observations.

Requires (run in order first):
  1. run_callmip_simulations.jl  → output_callmip_sims/
  2. analyze_posterior_ensemble.jl → output_posterior_analysis/
  3. emulate_sample.jl           → output_posterior_uq/ (for parameter plots)

Outputs (all in output_evaluation/):
  prior_vs_posterior_annual.png  — annual-mean time series NEE/Qle/Qh
  seasonal_cycles.png            — mean annual cycle prior/posterior/obs
  parameter_shifts.png           — parameter movement in units of prior σ
  uncertainty_reduction.png      — posterior σ / prior σ per parameter
  skill_score_table.txt          — text table of RMSE/bias/R

Usage (on a compute node):
    julia --project=experiments/callmip_uq_dk_sor \\
          experiments/callmip_uq_dk_sor/evaluate_calibration.jl
"""

println("[1/5] Loading packages..."); flush(stdout)
using JLD2, Dates, Statistics, NCDatasets
using Plots, Printf
println("[2/5] Packages loaded."); flush(stdout)

# ── Paths ─────────────────────────────────────────────────────────────────────
const climaland_dir = abspath(joinpath(@__DIR__, "..", ".."))
const exp_dir       = joinpath(climaland_dir, "experiments", "callmip_uq_dk_sor")
const sims_dir      = joinpath(exp_dir, "output_callmip_sims", "iteration_000")
const analysis_dir  = joinpath(exp_dir, "output_posterior_analysis")
const posterior_dir = joinpath(exp_dir, "output_posterior_uq")
const out_dir       = joinpath(exp_dir, "output_evaluation")
isdir(out_dir) || mkdir(out_dir)

const flux_nc_path = joinpath(
    climaland_dir, "DK_Sor",
    "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc",
)

const CAL_END_YEAR = 2012   # last year of calibration period (inclusive)

# mol CO₂ m⁻² s⁻¹ → gC m⁻² d⁻¹
const mol_to_gC = 12.0 * 86400.0

# ── Helpers ───────────────────────────────────────────────────────────────────
function load_sim(member::Int)
    path = joinpath(sims_dir, "member_$(lpad(member,3,'0'))", "callmip_diagnostics.jld2")
    isfile(path) || error("Missing: $path\nRun run_callmip_simulations.jl first.")
    JLD2.load(path)
end

# callmip_diagnostics.jld2 stores surface vars inside surface_data sub-dict
function get_surface(d, key)
    sd = d["surface_data"]   # Dict{String,Any}
    haskey(sd, key) ? Float64.(sd[key]) : fill(NaN, length(d["dates"]))
end

function rmse(pred, obs)
    mask = .!isnan.(pred) .& .!isnan.(obs)
    sum(mask) == 0 && return NaN
    sqrt(mean((pred[mask] .- obs[mask]).^2))
end

function bias(pred, obs)
    mask = .!isnan.(pred) .& .!isnan.(obs)
    sum(mask) == 0 && return NaN
    mean(pred[mask] .- obs[mask])
end

function pearsonr(pred, obs)
    mask = .!isnan.(pred) .& .!isnan.(obs)
    sum(mask) < 3 && return NaN
    cor(pred[mask], obs[mask])
end

function annual_mean(vals, dates)
    yrs = unique(year.(dates))
    ymean = [let v = filter(!isnan, vals[year.(dates) .== y]); isempty(v) ? NaN : mean(v); end for y in yrs]
    Float64.(yrs), ymean
end

function monthly_mean(vals, dates)
    [let v = filter(!isnan, vals[month.(dates) .== m]); isempty(v) ? NaN : mean(v); end for m in 1:12]
end

# Align two date-indexed vectors onto a shared date axis (by year(d) ≤ / >)
function date_align(src_vals, src_dates, tgt_dates)
    d2v = Dict(zip(src_dates, src_vals))
    [get(d2v, d, NaN) for d in tgt_dates]
end

# ── Load prior and posterior simulations ──────────────────────────────────────
println("[3/5] Loading simulation outputs..."); flush(stdout)

prior_d     = load_sim(1)
posterior_d = load_sim(2)

sim_dates_raw = Date.(prior_d["dates"])   # Vector{Date}
n_days = length(sim_dates_raw)
@info "Simulation period: $(first(sim_dates_raw)) – $(last(sim_dates_raw))  ($n_days days)"

nee_prior_raw = get_surface(prior_d,     "nee") .* mol_to_gC
let raw = get_surface(posterior_d, "nee")
    # Threshold for numerical blowup: >1e-3 mol m⁻² s⁻¹ is ~10³× realistic max.
    # EKI-optimal member has instability in Dec 2005–Mar 2006 and Dec 2009 (56 days).
    blowup_mask = abs.(raw) .> 1e-3
    n_blowup = sum(blowup_mask)
    if n_blowup > 0
        @warn "EKI posterior NEE: $n_blowup blowup day(s) (|NEE| > 1e-3 mol m⁻² s⁻¹) set to NaN — model numerical instability with calibrated parameters"
        raw_clean = copy(raw)
        raw_clean[blowup_mask] .= NaN
        global nee_post_raw = raw_clean .* mol_to_gC
    else
        global nee_post_raw = raw .* mol_to_gC
    end
end
qle_prior_raw = get_surface(prior_d,     "lhf")
qle_post_raw  = get_surface(posterior_d, "lhf")
qh_prior_raw  = get_surface(prior_d,     "shf")
qh_post_raw   = get_surface(posterior_d, "shf")

# ── Load FLUXNET observations ─────────────────────────────────────────────────
isfile(flux_nc_path) || error("FLUXNET file not found: $flux_nc_path")
obs_ds = NCDatasets.NCDataset(flux_nc_path, "r")
nee_obs_nc  = Float64.(NCDatasets.coalesce.(obs_ds["NEE_daily"][:], NaN))
qle_obs_nc  = Float64.(NCDatasets.coalesce.(obs_ds["Qle_daily"][:], NaN))
qh_obs_nc   = Float64.(NCDatasets.coalesce.(obs_ds["Qh_daily"][:], NaN))
obs_times_nc = obs_ds["time"][:]
NCDatasets.close(obs_ds)
obs_dates_nc = Date.(obs_times_nc)

# ── Load posterior uncertainty bands ─────────────────────────────────────────
bands = JLD2.load(joinpath(analysis_dir, "posterior_uncertainty_bands.jld2"))
band_dates = Date.(bands["dates"])   # 4018 days

nee_p05_raw = Float64.(bands["nee_p05"])   # already in gC m⁻² d⁻¹ (converted in analyze_posterior_ensemble.jl)
nee_p50_raw = Float64.(bands["nee_p50"])
nee_p95_raw = Float64.(bands["nee_p95"])
qle_p05_raw = Float64.(bands["qle_p05"])
qle_p50_raw = Float64.(bands["qle_p50"])
qle_p95_raw = Float64.(bands["qle_p95"])
qh_p05_raw  = Float64.(bands["qh_p05"])
qh_p50_raw  = Float64.(bands["qh_p50"])
qh_p95_raw  = Float64.(bands["qh_p95"])

# ── Align everything onto band_dates as the common time axis ─────────────────
common_dates = band_dates   # 4018 days, the best-covered axis

nee_prior = date_align(nee_prior_raw, sim_dates_raw, common_dates)
nee_post  = date_align(nee_post_raw,  sim_dates_raw, common_dates)
qle_prior = date_align(qle_prior_raw, sim_dates_raw, common_dates)
qle_post  = date_align(qle_post_raw,  sim_dates_raw, common_dates)
qh_prior  = date_align(qh_prior_raw,  sim_dates_raw, common_dates)
qh_post   = date_align(qh_post_raw,   sim_dates_raw, common_dates)

nee_obs = date_align(nee_obs_nc, obs_dates_nc, common_dates)
qle_obs = date_align(qle_obs_nc, obs_dates_nc, common_dates)
qh_obs  = date_align(qh_obs_nc,  obs_dates_nc, common_dates)

nee_p05 = nee_p05_raw
nee_p50 = nee_p50_raw
nee_p95 = nee_p95_raw
qle_p05 = qle_p05_raw
qle_p50 = qle_p50_raw
qle_p95 = qle_p95_raw
qh_p05  = qh_p05_raw
qh_p50  = qh_p50_raw
qh_p95  = qh_p95_raw

cal_mask = year.(common_dates) .<= CAL_END_YEAR
val_mask = .!cal_mask
@info "Cal days: $(sum(cal_mask))    Val days: $(sum(val_mask))"

# ── Skill scores ──────────────────────────────────────────────────────────────
println("[4/5] Computing skill scores..."); flush(stdout)

function skill_row(label, pred_all, obs_all)
    @sprintf("  %-28s  Cal: RMSE=%7.3f  bias=%+7.3f  R=%6.3f  |  Val: RMSE=%7.3f  bias=%+7.3f  R=%6.3f",
        label,
        rmse(pred_all[cal_mask], obs_all[cal_mask]),
        bias(pred_all[cal_mask], obs_all[cal_mask]),
        pearsonr(pred_all[cal_mask], obs_all[cal_mask]),
        rmse(pred_all[val_mask], obs_all[val_mask]),
        bias(pred_all[val_mask], obs_all[val_mask]),
        pearsonr(pred_all[val_mask], obs_all[val_mask]))
end

table_lines = [
    "ClimaLand DK-Sor — Calibration Skill Scores",
    "=" ^ 105,
    @sprintf("  %-28s  %-45s  %-45s", "Model run", "Calibration (≤$(CAL_END_YEAR))", "Validation (>$(CAL_END_YEAR))"),
    @sprintf("  %-28s  %-11s  %-8s  %-6s      %-11s  %-8s  %-6s",
        "", "RMSE", "bias", "R", "RMSE", "bias", "R"),
    "-" ^ 105,
    "NEE (gC m⁻² d⁻¹):",
    skill_row("Prior",      nee_prior, nee_obs),
    skill_row("Posterior (EKI)",    nee_post,  nee_obs),
    skill_row("Posterior (ens p50)", nee_p50,  nee_obs),
    "Qle (W m⁻²):",
    skill_row("Prior",      qle_prior, qle_obs),
    skill_row("Posterior (EKI)",    qle_post,  qle_obs),
    skill_row("Posterior (ens p50)", qle_p50,  qle_obs),
    "Qh (W m⁻²):",
    skill_row("Prior",      qh_prior, qh_obs),
    skill_row("Posterior (EKI)",    qh_post,  qh_obs),
    skill_row("Posterior (ens p50)", qh_p50,  qh_obs),
    "-" ^ 105,
    @sprintf("  Posterior 90%% coverage: NEE=%.1f%%  Qle=%.1f%%  Qh=%.1f%%",
        bands["coverage_nee_90pct"] * 100,
        bands["coverage_qle_90pct"] * 100,
        bands["coverage_qh_90pct"] * 100),
]

for l in table_lines; println(l); end
table_path = joinpath(out_dir, "skill_score_table.txt")
open(table_path, "w") do io; foreach(l -> println(io, l), table_lines); end
@info "Saved skill_score_table.txt"

# ── Plotting ──────────────────────────────────────────────────────────────────
println("[5/5] Generating plots..."); flush(stdout)
gr()

# Helper: annual mean of obs aligned on common_dates
yr_obs_nee_x, yr_obs_nee_y = annual_mean(nee_obs, common_dates)
yr_obs_qle_x, yr_obs_qle_y = annual_mean(qle_obs, common_dates)
yr_obs_qh_x,  yr_obs_qh_y  = annual_mean(qh_obs,  common_dates)

cal_end_x = CAL_END_YEAR + 0.5

# ── Plot 1: Annual mean time series ──────────────────────────────────────────
yr_nee_prior_x, yr_nee_prior_y = annual_mean(nee_prior, common_dates)
yr_nee_post_x,  yr_nee_post_y  = annual_mean(nee_post,  common_dates)
yr_nee_p50_x,   yr_nee_p50_y   = annual_mean(nee_p50,   common_dates)
yr_nee_p05_x,   yr_nee_p05_y   = annual_mean(nee_p05,   common_dates)
yr_nee_p95_x,   yr_nee_p95_y   = annual_mean(nee_p95,   common_dates)

yr_qle_prior_x, yr_qle_prior_y = annual_mean(qle_prior, common_dates)
yr_qle_post_x,  yr_qle_post_y  = annual_mean(qle_post,  common_dates)
yr_qle_p50_x,   yr_qle_p50_y   = annual_mean(qle_p50,   common_dates)
yr_qle_p05_x,   yr_qle_p05_y   = annual_mean(qle_p05,   common_dates)
yr_qle_p95_x,   yr_qle_p95_y   = annual_mean(qle_p95,   common_dates)

yr_qh_prior_x, yr_qh_prior_y = annual_mean(qh_prior, common_dates)
yr_qh_post_x,  yr_qh_post_y  = annual_mean(qh_post,  common_dates)
yr_qh_p50_x,   yr_qh_p50_y   = annual_mean(qh_p50,   common_dates)
yr_qh_p05_x,   yr_qh_p05_y   = annual_mean(qh_p05,   common_dates)
yr_qh_p95_x,   yr_qh_p95_y   = annual_mean(qh_p95,   common_dates)

p1 = plot(yr_nee_prior_x, yr_nee_prior_y;
    label="Prior", color=:gray, lw=2, ls=:dash,
    ylabel="NEE (gC m⁻² d⁻¹)", title="Annual Mean NEE",
    legend=:topright, xrotation=45)
plot!(p1, yr_nee_p05_x, yr_nee_p05_y; fillrange=yr_nee_p95_y,
    fillalpha=0.20, color=:dodgerblue, label="Ens 90% band", lw=0)
plot!(p1, yr_nee_p50_x, yr_nee_p50_y; label="Ens median", color=:dodgerblue, lw=1.5, ls=:dot)
plot!(p1, yr_nee_post_x, yr_nee_post_y; label="Posterior (EKI)", color=:blue, lw=2)
scatter!(p1, yr_obs_nee_x, yr_obs_nee_y; label="FLUXNET obs", ms=4, mc=:black, ma=0.8)
vline!(p1, [cal_end_x]; color=:black, lw=1, ls=:dot, label="Cal/Val")

p2 = plot(yr_qle_prior_x, yr_qle_prior_y;
    label="Prior", color=:gray, lw=2, ls=:dash,
    ylabel="Qle (W m⁻²)", title="Annual Mean Latent Heat",
    legend=:topright)
plot!(p2, yr_qle_p05_x, yr_qle_p05_y; fillrange=yr_qle_p95_y,
    fillalpha=0.20, color=:orange, label="Ens 90% band", lw=0)
plot!(p2, yr_qle_p50_x, yr_qle_p50_y; label="Ens median", color=:orange, lw=1.5, ls=:dot)
plot!(p2, yr_qle_post_x, yr_qle_post_y; label="Posterior (EKI)", color=:darkorange, lw=2)
scatter!(p2, yr_obs_qle_x, yr_obs_qle_y; label="FLUXNET obs", ms=4, mc=:black, ma=0.8)
vline!(p2, [cal_end_x]; color=:black, lw=1, ls=:dot, label="Cal/Val")

p3 = plot(yr_qh_prior_x, yr_qh_prior_y;
    label="Prior", color=:gray, lw=2, ls=:dash,
    ylabel="Qh (W m⁻²)", title="Annual Mean Sensible Heat",
    xlabel="Year", legend=:topright)
plot!(p3, yr_qh_p05_x, yr_qh_p05_y; fillrange=yr_qh_p95_y,
    fillalpha=0.20, color=:red, label="Ens 90% band", lw=0)
plot!(p3, yr_qh_p50_x, yr_qh_p50_y; label="Ens median", color=:red, lw=1.5, ls=:dot)
plot!(p3, yr_qh_post_x, yr_qh_post_y; label="Posterior (EKI)", color=:crimson, lw=2)
scatter!(p3, yr_obs_qh_x, yr_obs_qh_y; label="FLUXNET obs", ms=4, mc=:black, ma=0.8)
vline!(p3, [cal_end_x]; color=:black, lw=1, ls=:dot, label="Cal/Val")

pall = plot(p1, p2, p3; layout=(3,1), size=(1200, 900),
    plot_title="ClimaLand DK-Sor: Prior vs Posterior vs FLUXNET (annual means)")
savefig(pall, joinpath(out_dir, "prior_vs_posterior_annual.png"))
@info "Saved prior_vs_posterior_annual.png"

# ── Plot 2: Seasonal cycles ────────────────────────────────────────────────────
mnths     = 1:12
mn_labels = ["J","F","M","A","M","J","J","A","S","O","N","D"]

pc1 = plot(mnths, monthly_mean(nee_prior, common_dates);
    label="Prior", color=:gray, lw=2, ls=:dash,
    ylabel="NEE (gC m⁻² d⁻¹)", title="Seasonal Cycle — NEE",
    xticks=(mnths, mn_labels), legend=:bottomleft)
mm_p05 = monthly_mean(nee_p05, common_dates)
mm_p95 = monthly_mean(nee_p95, common_dates)
plot!(pc1, mnths, mm_p05; fillrange=mm_p95, fillalpha=0.20, color=:dodgerblue, label="Ens 90%", lw=0)
plot!(pc1, mnths, monthly_mean(nee_p50, common_dates); label="Ens median", color=:dodgerblue, lw=1.5, ls=:dot)
plot!(pc1, mnths, monthly_mean(nee_post, common_dates); label="Posterior", color=:blue, lw=2.5)
plot!(pc1, mnths, monthly_mean(nee_obs,  common_dates); label="FLUXNET", color=:black, lw=2, ls=:dashdot)

pc2 = plot(mnths, monthly_mean(qle_prior, common_dates);
    label="Prior", color=:gray, lw=2, ls=:dash,
    ylabel="Qle (W m⁻²)", title="Seasonal Cycle — Latent Heat",
    xticks=(mnths, mn_labels))
qle_p05_m = monthly_mean(qle_p05, common_dates)
qle_p95_m = monthly_mean(qle_p95, common_dates)
plot!(pc2, mnths, qle_p05_m; fillrange=qle_p95_m, fillalpha=0.20, color=:orange, label="Ens 90%", lw=0)
plot!(pc2, mnths, monthly_mean(qle_p50, common_dates); label="Ens median", color=:orange, lw=1.5, ls=:dot)
plot!(pc2, mnths, monthly_mean(qle_post, common_dates); label="Posterior", color=:darkorange, lw=2.5)
plot!(pc2, mnths, monthly_mean(qle_obs,  common_dates); label="FLUXNET", color=:black, lw=2, ls=:dashdot)

pc3 = plot(mnths, monthly_mean(qh_prior, common_dates);
    label="Prior", color=:gray, lw=2, ls=:dash,
    ylabel="Qh (W m⁻²)", title="Seasonal Cycle — Sensible Heat",
    xticks=(mnths, mn_labels), xlabel="Month")
qh_p05_m = monthly_mean(qh_p05, common_dates)
qh_p95_m = monthly_mean(qh_p95, common_dates)
plot!(pc3, mnths, qh_p05_m; fillrange=qh_p95_m, fillalpha=0.20, color=:red, label="Ens 90%", lw=0)
plot!(pc3, mnths, monthly_mean(qh_p50, common_dates); label="Ens median", color=:red, lw=1.5, ls=:dot)
plot!(pc3, mnths, monthly_mean(qh_post, common_dates); label="Posterior", color=:crimson, lw=2.5)
plot!(pc3, mnths, monthly_mean(qh_obs,  common_dates); label="FLUXNET", color=:black, lw=2, ls=:dashdot)

pseas = plot(pc1, pc2, pc3; layout=(1,3), size=(1400, 420),
    plot_title="ClimaLand DK-Sor: Mean Annual Cycle (1997–2013)")
savefig(pseas, joinpath(out_dir, "seasonal_cycles.png"))
@info "Saved seasonal_cycles.png"

# ── Plot 3: Parameter shifts ───────────────────────────────────────────────────
post_jld    = JLD2.load(joinpath(posterior_dir, "posterior_its1to8.jld2"))
param_names = String.(post_jld["param_names"])
post_mean   = vec(mean(post_jld["constrained_posterior"], dims=2))
post_std    = vec(std( post_jld["constrained_posterior"], dims=2))
ekp_optimal = Float64.(post_jld["constrained_ekp_optimal"])

# Prior means and σ — must stay in sync with build_dk_sor_priors() in
# experiments/calibrate_dk_sor/priors.jl.  Order matches param_names above.
# These are the constrained_gaussian(μ, σ, lb, ub) arguments — approximate
# constrained-space mean and standard deviation for each parameter.
n_params   = length(param_names)
prior_means = [0.5, 0.43, 20.0, 0.1, 0.10, 0.001, 0.65, 0.85, 0.97,
               25000.0, 0.01, 0.01, 61000.0, 1.0, 0.1, 2500.0]
prior_sigmas = [0.3, 0.15, 8.0, 0.05, 0.04, 0.0005, 0.12, 0.25, 0.02,
                10000.0, 0.005, 0.005, 15000.0, 0.5, 0.05, 1500.0]
@assert length(prior_means)  == n_params "prior_means length mismatch — update to match priors.jl"
@assert length(prior_sigmas) == n_params "prior_sigmas length mismatch — update to match priors.jl"
shift_ekp  = (ekp_optimal .- prior_means) ./ prior_sigmas
shift_post = (post_mean   .- prior_means) ./ prior_sigmas
xs = 1:n_params
bar_w = 0.35

pshift = bar(xs .- bar_w/2, shift_ekp;
    bar_width=bar_w, label="EKI optimal", color=:steelblue, alpha=0.8,
    ylabel="Δ (units of prior σ)", title="Parameter Shifts from Prior",
    xticks=(xs, param_names), xrotation=40, size=(1400, 520), legend=:topright)
bar!(pshift, xs .+ bar_w/2, shift_post;
    bar_width=bar_w, label="MCMC mean", color=:darkorange, alpha=0.8)
hline!(pshift, [0.0];     color=:black, lw=1,   ls=:dash, label="Prior mean")
hline!(pshift, [-2.0, 2.0]; color=:gray, lw=0.8, ls=:dot,  label="±2σ")
savefig(pshift, joinpath(out_dir, "parameter_shifts.png"))
@info "Saved parameter_shifts.png"

# ── Plot 4: Uncertainty reduction ─────────────────────────────────────────────
sigma_ratio = post_std ./ prior_sigmas

pred_ur = bar(xs, sigma_ratio;
    ylabel="Posterior σ / Prior σ", title="Parameter Uncertainty Reduction (< 1 = data-informed)",
    xticks=(xs, param_names), xrotation=40,
    color=:mediumseagreen, alpha=0.85, label="Post/Prior σ ratio",
    size=(1400, 520))
hline!(pred_ur, [1.0]; color=:black, lw=1.5, ls=:dash, label="Prior level")
savefig(pred_ur, joinpath(out_dir, "uncertainty_reduction.png"))
@info "Saved uncertainty_reduction.png"

# ── Summary ────────────────────────────────────────────────────────────────────
println()
println("╔══════════════════════════════════════════════════════════════╗")
println("║          Evaluation complete — DK-Sor CalLMIP Phase 1a       ║")
println("╠══════════════════════════════════════════════════════════════╣")
@printf("║  Posterior 90%% coverage: NEE=%5.1f%%  Qle=%5.1f%%  Qh=%5.1f%%    ║\n",
    bands["coverage_nee_90pct"] * 100,
    bands["coverage_qle_90pct"] * 100,
    bands["coverage_qh_90pct"] * 100)
println("╠══════════════════════════════════════════════════════════════╣")
println("║  Outputs saved to output_evaluation/:                        ║")
println("║    prior_vs_posterior_annual.png                             ║")
println("║    seasonal_cycles.png                                       ║")
println("║    parameter_shifts.png                                      ║")
println("║    uncertainty_reduction.png                                 ║")
println("║    skill_score_table.txt                                     ║")
println("╚══════════════════════════════════════════════════════════════╝")
