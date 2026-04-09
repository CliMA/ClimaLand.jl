"""
Calibration quality diagnostics for DK-Sor EKI run.

Compares prior ensemble (iter_000) vs posterior ensemble (iter_008) against
observed NEE, Qle, Qh.  Uses:
  - G_ensemble.jld2 per iteration for RMSE curve + obs-day statistics
  - daily_diagnostics.jld2 per member for the full time-series seasonal cycle
  - observations.jld2 for observed values and dates

Outputs:
  calibration_check.png  — 4-panel figure:
    top-3: monthly seasonality obs vs posterior (with prior shaded), per variable
    bottom: RMSE vs EKI iteration for each variable

Usage:
    julia --project=.buildkite experiments/calibrate_dk_sor/plot_calibration_check.jl
"""

using CairoMakie
using Dates
using Statistics: mean, std, quantile, median
import JLD2

# ── Paths ─────────────────────────────────────────────────────────────────────

const REPO     = abspath(joinpath(@__DIR__, "..", ".."))
const OBS_PATH = joinpath(REPO, "experiments/calibrate_dk_sor/observations.jld2")
const OUT_DIR  = joinpath(REPO, "experiments/calibrate_dk_sor")
const ITER_DIR = joinpath(REPO, "experiments/calibrate_dk_sor/output")

# ── Load observations ─────────────────────────────────────────────────────────

println("Loading observations…")
obs   = JLD2.load(OBS_PATH)
y_obs = obs["y_obs"]::Vector{Float64}
obs_dates = obs["obs_dates"]::Vector{Date}
n_obs = length(obs_dates)
@assert length(y_obs) == 3 * n_obs

y_nee_obs = y_obs[1:n_obs]
y_qle_obs = y_obs[n_obs+1:2n_obs]
y_qh_obs  = y_obs[2n_obs+1:end]

println("  $n_obs observation days, $(first(obs_dates)) – $(last(obs_dates))")

# ── Load G_ensemble per iteration ─────────────────────────────────────────────

println("Loading G_ensemble files…")

iters_available = Int[]
for i in 0:20
    p = joinpath(ITER_DIR, "iteration_$(lpad(i, 3, '0'))", "G_ensemble.jld2")
    isfile(p) && push!(iters_available, i)
end
println("  Found G_ensemble for iterations: $iters_available")

function load_G(iter::Int)
    p = joinpath(ITER_DIR, "iteration_$(lpad(iter, 3, '0'))", "G_ensemble.jld2")
    G = JLD2.load(p, "single_stored_object")::Matrix{Float64}
    @assert size(G, 1) == 3 * n_obs "Expected $(3n_obs) rows, got $(size(G,1))"
    return G
end

# RMSE helper (per-variable, per-iteration)
function compute_rmse(G::Matrix{Float64}, y::Vector{Float64})
    # ensemble-mean residuals
    gm = vec(mean(G, dims=2))
    ok = .!isnan.(gm)
    return sqrt(mean((gm[ok] .- y[ok]).^2))
end

# Robust prior RMSE: median across member-wise RMSEs, excluding blown-up members
function compute_rmse_robust(G::Matrix{Float64}, y::Vector{Float64})
    member_rmse = [begin
        g = G[:, m]
        ok = .!isnan.(g)
        sqrt(mean((g[ok] .- y[ok]).^2))
    end for m in 1:size(G,2)]
    # exclude extreme outliers (member RMSE > 10x median)
    med = median(member_rmse)
    valid_m = filter(r -> r <= 10 * med, member_rmse)
    return median(valid_m)
end

rmse_nee = Float64[]; rmse_qle = Float64[]; rmse_qh = Float64[]
for (idx, i) in enumerate(iters_available)
    G   = load_G(i)
    push!(rmse_nee, compute_rmse(G[1:n_obs,         :], y_nee_obs))
    push!(rmse_qle, compute_rmse(G[n_obs+1:2n_obs,   :], y_qle_obs))
    push!(rmse_qh,  compute_rmse(G[2n_obs+1:end,     :], y_qh_obs))
end

# Robust prior RMSE from iter_000
let G0 = load_G(iters_available[1])
    global rmse_nee_robust0 = compute_rmse_robust(G0[1:n_obs,         :], y_nee_obs)
    global rmse_qle_robust0 = compute_rmse_robust(G0[n_obs+1:2n_obs,   :], y_qle_obs)
    global rmse_qh_robust0  = compute_rmse_robust(G0[2n_obs+1:end,     :], y_qh_obs)
end
println("  Robust prior RMSE (median member, outliers filtered):")
println("    NEE=$(round(rmse_nee_robust0,digits=2))  Qle=$(round(rmse_qle_robust0,digits=1))  Qh=$(round(rmse_qh_robust0,digits=1))")
println("  NEE  RMSE iter $(iters_available[1])→$(iters_available[end]): " *
        "$(round(rmse_nee[1], sigdigits=4)) → $(round(rmse_nee[end], sigdigits=4))")
println("  Qle  RMSE iter $(iters_available[1])→$(iters_available[end]): " *
        "$(round(rmse_qle[1], sigdigits=4)) → $(round(rmse_qle[end], sigdigits=4))")
println("  Qh   RMSE iter $(iters_available[1])→$(iters_available[end]): " *
        "$(round(rmse_qh[1], sigdigits=4)) → $(round(rmse_qh[end], sigdigits=4))")

# ── Load daily diagnostics: all members, prior (iter_000) + posterior (iter_008) ─

println("Loading daily diagnostics for prior and posterior ensembles…")

function load_all_members_daily(iter::Int; filter_nee_outliers::Bool=false)
    istr = lpad(iter, 3, '0')
    nees = Vector{Float64}[]; qles = Vector{Float64}[]; qhs = Vector{Float64}[]
    dts_ref = nothing
    n_skipped = 0
    for m in 1:50
        mstr = lpad(m, 3, '0')
        p = joinpath(ITER_DIR, "iteration_$(istr)", "member_$(mstr)", "daily_diagnostics.jld2")
        isfile(p) || break
        d    = JLD2.load(p)
        dts  = d["dates"]::Vector{Date}
        nee  = d["nee"]::Vector{Float64} .* (12.0 * 86400.0)  # mol/m²/s → gC/m²/d
        qle  = d["qle"]::Vector{Float64}
        qh   = d["qh"]::Vector{Float64}
        # Optionally skip blown-up members (|mean NEE| > 100 gC/m²/d)
        if filter_nee_outliers && abs(mean(nee)) > 100.0
            n_skipped += 1
            continue
        end
        if dts_ref === nothing; dts_ref = dts; end
        push!(nees, nee); push!(qles, qle); push!(qhs, qh)
    end
    isempty(nees) && error("No member daily_diagnostics found for iter $iter")
    N = length(nees)
    skip_str = n_skipped > 0 ? " ($n_skipped outlier members excluded)" : ""
    println("  iter $iter: $N members, $(length(dts_ref)) days, " *
            "$(first(dts_ref)) – $(last(dts_ref))$(skip_str)")
    # Convert to (n_days × N_ens) matrices
    nee_mat = hcat(nees...)
    qle_mat = hcat(qles...)
    qh_mat  = hcat(qhs...)
    return dts_ref, nee_mat, qle_mat, qh_mat
end

dts_prior, nee_prior, qle_prior, qh_prior = load_all_members_daily(iters_available[1];   filter_nee_outliers=true)
dts_post,  nee_post,  qle_post,  qh_post  = load_all_members_daily(iters_available[end]; filter_nee_outliers=false)

# ── Monthly climatology helper ─────────────────────────────────────────────────

function monthly_clim(dates::Vector{Date}, mat::Matrix{Float64};
                      cal_start = Date(2004, 1, 1),
                      cal_end   = Date(2013, 12, 31))
    cal_mask = cal_start .<= dates .<= cal_end
    d_cal    = dates[cal_mask]
    m_cal    = mat[cal_mask, :]
    # Ensemble mean per day
    ens_mean = vec(mean(m_cal, dims = 2))
    ens_p5   = vec(mapslices(x -> quantile(x, 0.05), m_cal, dims = 2))
    ens_p95  = vec(mapslices(x -> quantile(x, 0.95), m_cal, dims = 2))
    by_month_mean  = [Float64[] for _ in 1:12]
    by_month_p5    = [Float64[] for _ in 1:12]
    by_month_p95   = [Float64[] for _ in 1:12]
    for (i, dt) in enumerate(d_cal)
        mo = Dates.month(dt)
        isnan(ens_mean[i])  || push!(by_month_mean[mo],  ens_mean[i])
        isnan(ens_p5[i])    || push!(by_month_p5[mo],    ens_p5[i])
        isnan(ens_p95[i])   || push!(by_month_p95[mo],   ens_p95[i])
    end
    mn   = [isempty(v) ? NaN : mean(v) for v in by_month_mean]
    p5   = [isempty(v) ? NaN : mean(v) for v in by_month_p5]
    p95  = [isempty(v) ? NaN : mean(v) for v in by_month_p95]
    return mn, p5, p95
end

function monthly_obs(dates::Vector{Date}, vals::Vector{Float64})
    by_month = [Float64[] for _ in 1:12]
    for (dt, v) in zip(dates, vals)
        isnan(v) || push!(by_month[Dates.month(dt)], v)
    end
    mn = [isempty(v) ? NaN : mean(v)   for v in by_month]
    sd = [length(v) < 2 ? NaN : std(v) for v in by_month]
    return mn, sd
end

println("Computing monthly climatologies…")
MONTHS = 1:12
MO_LAB = ["J","F","M","A","M","J","J","A","S","O","N","D"]

nee_prior_mn,  nee_prior_p5,  nee_prior_p95  = monthly_clim(dts_prior, nee_prior)
qle_prior_mn,  qle_prior_p5,  qle_prior_p95  = monthly_clim(dts_prior, qle_prior)
qh_prior_mn,   qh_prior_p5,   qh_prior_p95   = monthly_clim(dts_prior, qh_prior)

nee_post_mn,   nee_post_p5,   nee_post_p95   = monthly_clim(dts_post,  nee_post)
qle_post_mn,   qle_post_p5,   qle_post_p95   = monthly_clim(dts_post,  qle_post)
qh_post_mn,    qh_post_p5,    qh_post_p95    = monthly_clim(dts_post,  qh_post)

nee_obs_mn,    nee_obs_sd                     = monthly_obs(obs_dates, y_nee_obs)
qle_obs_mn,    qle_obs_sd                     = monthly_obs(obs_dates, y_qle_obs)
qh_obs_mn,     qh_obs_sd                      = monthly_obs(obs_dates, y_qh_obs)

# ── Build figure ──────────────────────────────────────────────────────────────

println("Building figure…")

fig = Figure(size = (1100, 1050), fontsize = 13)

C_obs    = :black
C_prior  = (:steelblue,  0.55)
C_post   = (:firebrick,  0.55)
C_obs_bd = (:black,      0.30)

function add_seasonality_panel!(fig, row, col,
                                obs_mn, obs_sd,
                                prior_mn, prior_p5, prior_p95,
                                post_mn,  post_p5,  post_p95,
                                ylabel, var_label;
                                rmse_prior, rmse_post)
    title_str = "$var_label  |  Prior RMSE = $(round(rmse_prior, digits=2))  →  " *
                "Posterior RMSE = $(round(rmse_post, digits=2))"
    ax = Axis(fig[row, col];
              title   = title_str,
              xlabel  = row == 3 ? "Month" : "",
              ylabel  = ylabel,
              xticks  = (collect(1:12), MO_LAB))

    x = collect(Float64, 1:12)

    # obs ± 1 SD band
    obs_lo = obs_mn .- obs_sd; obs_hi = obs_mn .+ obs_sd
    band!(ax, x, obs_lo, obs_hi; color = C_obs_bd)

    # prior 5th–95th band
    band!(ax, x, prior_p5, prior_p95; color = (:steelblue, 0.25))
    lines!(ax, x, prior_mn; color = :steelblue, linewidth = 2,
           label = "Prior (iter $(iters_available[1]))")

    # posterior 5th–95th band
    band!(ax, x, post_p5, post_p95; color = (:firebrick, 0.25))
    lines!(ax, x, post_mn; color = :firebrick, linewidth = 2,
           label = "Posterior (iter $(iters_available[end]))")

    # obs mean
    lines!(ax, x, obs_mn; color = :black, linewidth = 2.5, label = "Observed")
    scatter!(ax, x, obs_mn; color = :black, markersize = 7)

    axislegend(ax; position = :lt, framevisible = false, labelsize = 11)
    return ax
end

# Row 1: NEE
add_seasonality_panel!(
    fig, 1, 1,
    nee_obs_mn, nee_obs_sd,
    nee_prior_mn, nee_prior_p5, nee_prior_p95,
    nee_post_mn,  nee_post_p5,  nee_post_p95,
    "NEE (gC m⁻² d⁻¹)", "NEE";
    rmse_prior = rmse_nee_robust0, rmse_post = rmse_nee[end])

# Row 2: Qle
add_seasonality_panel!(
    fig, 2, 1,
    qle_obs_mn, qle_obs_sd,
    qle_prior_mn, qle_prior_p5, qle_prior_p95,
    qle_post_mn,  qle_post_p5,  qle_post_p95,
    "Qle (W m⁻²)", "Latent Heat (Qle)";
    rmse_prior = rmse_qle_robust0, rmse_post = rmse_qle[end])

# Row 3: Qh
add_seasonality_panel!(
    fig, 3, 1,
    qh_obs_mn, qh_obs_sd,
    qh_prior_mn, qh_prior_p5, qh_prior_p95,
    qh_post_mn,  qh_post_p5,  qh_post_p95,
    "Qh (W m⁻²)", "Sensible Heat (Qh)";
    rmse_prior = rmse_qh_robust0, rmse_post = rmse_qh[end])

# Row 4: RMSE vs iteration
ax_rmse = Axis(fig[4, 1];
               title   = "Ensemble-mean RMSE vs EKI iteration  (log scale)",
               xlabel  = "EKI iteration",
               ylabel  = "RMSE (log scale)",
               yscale  = log10)

x_iter = Float64.(iters_available)
lines!(ax_rmse, x_iter, rmse_nee; color = :seagreen,  linewidth = 2.5,
       label = "NEE (gC m⁻² d⁻¹)")
scatter!(ax_rmse, x_iter, rmse_nee; color = :seagreen,  markersize = 8)

lines!(ax_rmse, x_iter, rmse_qle; color = :royalblue, linewidth = 2.5,
       label = "Qle (W m⁻²)")
scatter!(ax_rmse, x_iter, rmse_qle; color = :royalblue, markersize = 8)

lines!(ax_rmse, x_iter, rmse_qh;  color = :darkorange, linewidth = 2.5,
       label = "Qh (W m⁻²)")
scatter!(ax_rmse, x_iter, rmse_qh;  color = :darkorange, markersize = 8)

axislegend(ax_rmse; position = :rt, framevisible = false)

n_ens_str = string(size(load_G(iters_available[1]), 2))
Label(fig[0, 1]; text = "DK-Sor EKI Calibration Check  |  " *
      "$(first(obs_dates)) – $(last(obs_dates))  |  N_ens = $(n_ens_str)",
      fontsize = 15, font = :bold)

rowgap!(fig.layout, 8)

out_path = joinpath(OUT_DIR, "calibration_check.png")
save(out_path, fig; px_per_unit = 2)
println("Saved → $out_path")
