"""
Diagnostic plots for NEON site EKI calibration results.

Reads eki_file.jld2 (at OUTPUT_DIR) for all observation and model data:
  - y_obs — observation vector embedded in the EKP ObservationSeries
  - G     — G-matrices from each completed iteration

Three output figures:
  1. plot_timeseries.png         — daily soil CO₂ obs / prior / posterior
  2. plot_rmse_per_iteration.png — ensemble RMSE per iteration
  3. plot_1to1.png               — 1-to-1 scatter with regression lines

Configuration via environment variables:
    NEON_SITE_ID — NEON site ID (default: "NEON-srer")

Usage:
    julia --project=.buildkite experiments/calibrate_neon/plot_eki_diagnostics.jl
"""

using CairoMakie
using Dates
using Statistics
using LinearAlgebra
import JLD2
import ClimaLand

# ── Paths ─────────────────────────────────────────────────────────────────────

const SITE_ID = get(ENV, "NEON_SITE_ID", "NEON-srer")
const climaland_dir = pkgdir(ClimaLand)
const EKI_PATH = joinpath(climaland_dir, "experiments/calibrate_neon/output/iteration_006/eki_file.jld2")
const OBS_PATH = joinpath(climaland_dir, "experiments/calibrate_neon/observations.jld2")
const OUTDIR   = joinpath(climaland_dir, "experiments/calibrate_neon")

# ── Load all data from eki_file.jld2 ─────────────────────────────────────────

"""
Extract y_obs vector and G-matrices from the EKP's JLD2 file.

EKP field layout (positional):
  1 = u,  2 = observation_series,  3 = N_ens,  4 = g, ...

ObservationSeries.observations[1].samples[1] is the y vector.
G[i].stored_data is the (n_obs × N_ens) output matrix.
"""
function load_eki_data()
    f   = JLD2.jldopen(EKI_PATH, "r")
    obj = f["single_stored_object"]
    fld = getfield(obj, :fields)

    # ── y_obs from observation_series ────────────────────────────────────────
    obs_series = fld[2]
    os_fld     = getfield(obs_series, :fields)
    obs_vec    = os_fld[1]
    ob1_fld    = getfield(obs_vec[1], :fields)
    y_obs      = ob1_fld[1][1] :: Vector{Float64}

    # ── G-matrices ────────────────────────────────────────────────────────────
    g_containers = fld[4]
    n_iters      = length(g_containers)
    G = [getfield(g_containers[i], :fields)[1] :: Matrix{Float64}
         for i in 1:n_iters]

    JLD2.close(f)
    return y_obs, G
end

# ── Load observation dates from JLD2 ──────────────────────────────────────────

function load_obs_dates()
    obs_data = JLD2.load(OBS_PATH)
    return obs_data["obs_dates"]
end

# ── Helpers ───────────────────────────────────────────────────────────────────

function linreg(x, y)
    mask = .!isnan.(x) .& .!isnan.(y)
    xm, ym = x[mask], y[mask]
    n   = length(xm)
    sx  = sum(xm); sy = sum(ym)
    sxx = dot(xm, xm); sxy = dot(xm, ym)
    b   = (n*sxy - sx*sy) / (n*sxx - sx^2)
    a   = (sy - b*sx) / n
    return a, b
end

function monthly_mean(dates, values)
    by_month = [Float64[] for _ in 1:12]
    for (d, v) in zip(dates, values)
        isnan(v) || push!(by_month[Dates.month(d)], v)
    end
    means = [isempty(v) ? NaN : mean(v) for v in by_month]
    stds  = [length(v) < 2 ? NaN : std(v) for v in by_month]
    return means, stds
end

function ensemble_rmse(G_mat, y)
    N_ens = size(G_mat, 2)
    [let d = G_mat[:, m] .- y; ok = .!isnan.(d); sqrt(mean(d[ok].^2)) end
     for m in 1:N_ens]
end

# ─────────────────────────────────────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────────────────────────────────────

println("Loading EKP data from eki_file.jld2…")
y_obs, G_all = load_eki_data()
n_iters      = length(G_all)
n_obs        = length(y_obs)
println("  y_obs length: $n_obs observation days")
println("  Completed iterations: $n_iters")

println("Loading observation dates…")
obs_dates = load_obs_dates()
@assert length(obs_dates) == n_obs "Date count mismatch: $(length(obs_dates)) vs $n_obs"

G_prior = G_all[1]
G_post  = G_all[end]

prior_mean = vec(mean(G_prior, dims=2))
post_mean  = vec(mean(G_post, dims=2))

iter_label_post = n_iters - 1

months_label = ["Jan","Feb","Mar","Apr","May","Jun",
                "Jul","Aug","Sep","Oct","Nov","Dec"]

# ─────────────────────────────────────────────────────────────────────────────
#  PLOT 1 — Timeseries: daily soil CO₂ obs / prior / posterior
# ─────────────────────────────────────────────────────────────────────────────

println("Plotting timeseries…")

fig1 = Figure(size=(1200, 500))
ax1 = Axis(fig1[1, 1];
    xlabel = "Day index",
    ylabel = "Soil CO₂ (ppm)",
    title  = "$SITE_ID: Daily Soil CO₂ at ~6 cm — Prior vs Posterior vs Obs",
)

lines!(ax1, 1:n_obs, y_obs;       color=:black,     linewidth=2.5, label="NEON Observations")
lines!(ax1, 1:n_obs, prior_mean;  color=:royalblue, linewidth=2,   label="Prior (iter 0)")
lines!(ax1, 1:n_obs, post_mean;   color=:firebrick, linewidth=2,
       label="Posterior (iter $iter_label_post)")

Legend(fig1[1, 2], ax1; framevisible=false)

ts_path = joinpath(OUTDIR, "plot_timeseries.png")
save(ts_path, fig1; px_per_unit=2)
println("  Saved → $ts_path")

# ─────────────────────────────────────────────────────────────────────────────
#  PLOT 1b — Monthly seasonality
# ─────────────────────────────────────────────────────────────────────────────

println("Plotting seasonality…")

fig1b = Figure(size=(800, 500))
ax1b = Axis(fig1b[1, 1];
    xlabel = "Month",
    ylabel = "Soil CO₂ (ppm)",
    title  = "$SITE_ID: Monthly Mean Soil CO₂ Seasonality",
    xticks = (1:12, months_label),
)

obs_m,   obs_s   = monthly_mean(obs_dates, y_obs)
prior_m, prior_s = monthly_mean(obs_dates, prior_mean)
post_m,  post_s  = monthly_mean(obs_dates, post_mean)
xs = 1:12

band!(ax1b, xs, obs_m   .- obs_s,   obs_m   .+ obs_s;   color=(:black,     0.12))
band!(ax1b, xs, prior_m .- prior_s, prior_m .+ prior_s; color=(:royalblue, 0.15))
band!(ax1b, xs, post_m  .- post_s,  post_m  .+ post_s;  color=(:firebrick, 0.15))

lines!(ax1b, xs, obs_m;   color=:black,     linewidth=2.5, label="Observations")
lines!(ax1b, xs, prior_m; color=:royalblue, linewidth=2,   label="Prior (iter 0)")
lines!(ax1b, xs, post_m;  color=:firebrick, linewidth=2,
       label="Posterior (iter $iter_label_post)")

scatter!(ax1b, xs, obs_m;   color=:black,     markersize=8)
scatter!(ax1b, xs, prior_m; color=:royalblue, markersize=6)
scatter!(ax1b, xs, post_m;  color=:firebrick, markersize=6)

Legend(fig1b[1, 2], ax1b; framevisible=false)

seasonality_path = joinpath(OUTDIR, "plot_seasonality.png")
save(seasonality_path, fig1b; px_per_unit=2)
println("  Saved → $seasonality_path")

# ─────────────────────────────────────────────────────────────────────────────
#  PLOT 2 — RMSE per iteration
# ─────────────────────────────────────────────────────────────────────────────

println("Plotting RMSE per iteration…")

all_rmse = [ensemble_rmse(G_all[i], y_obs) for i in 1:n_iters]
for (i, r) in enumerate(all_rmse)
    println("  Iter $(i-1): mean RMSE=$(round(mean(r), digits=2)), std=$(round(std(r), digits=2))")
end

fig2 = Figure(size=(max(600, 300 + 200*n_iters), 500))
ax2  = Axis(fig2[1, 1];
            xlabel = "EKI Iteration",
            ylabel = "RMSE (Soil CO₂ ppm)",
            title  = "$SITE_ID: Ensemble RMSE by Iteration",
            xticks = (0:n_iters-1, string.(0:n_iters-1)))

colors_iter = [i == 1       ? :royalblue :
               i == n_iters ? :firebrick : :gray
               for i in 1:n_iters]

for (i, rmse_v) in enumerate(all_rmse)
    xi     = i - 1
    col    = colors_iter[i]
    jitter = randn(Float32, length(rmse_v)) .* 0.05f0
    scatter!(ax2, fill(xi, length(rmse_v)) .+ jitter, rmse_v;
             color=(col, 0.55), markersize=7,
             label = i == 1 ? "Ensemble members" : "")
    μ = mean(rmse_v); σ = std(rmse_v)
    errorbars!(ax2, [xi], [μ], [σ]; color=col, linewidth=2, whiskerwidth=12)
    scatter!(ax2, [xi], [μ];
             color=col, markersize=14, marker=:diamond,
             label = i == 1       ? "Prior mean (iter 0)" :
                     i == n_iters ? "Posterior mean (iter $(i-1))" :
                                    "Mean (iter $(i-1))")
end

Legend(fig2[1, 2], ax2; framevisible=false)

rmse_path = joinpath(OUTDIR, "plot_rmse_per_iteration.png")
save(rmse_path, fig2; px_per_unit=2)
println("  Saved → $rmse_path")

# ─────────────────────────────────────────────────────────────────────────────
#  PLOT 3 — 1-to-1 scatter
# ─────────────────────────────────────────────────────────────────────────────

println("Plotting 1-to-1…")

fig3 = Figure(size=(600, 500))

ax3 = Axis(fig3[1, 1];
    xlabel = "Observed Soil CO₂ (ppm)",
    ylabel = "Modelled Soil CO₂ (ppm)",
    title  = "$SITE_ID: Soil CO₂ 1-to-1",
    aspect = DataAspect(),
)

vmin = min(minimum(filter(!isnan, y_obs)),
           minimum(filter(!isnan, vcat(prior_mean, post_mean))))
vmax = max(maximum(filter(!isnan, y_obs)),
           maximum(filter(!isnan, vcat(prior_mean, post_mean))))
pad  = 0.06 * (vmax - vmin)
xl   = (vmin - pad, vmax + pad)
lines!(ax3, collect(xl), collect(xl);
       color=:black, linewidth=1.5, linestyle=:dash, label="1:1")

# Prior
scatter!(ax3, y_obs, prior_mean; color=(:royalblue, 0.35), markersize=4)
a_p, b_p = linreg(y_obs, prior_mean)
xs_r = range(xl[1], xl[2]; length=200)
lines!(ax3, xs_r, a_p .+ b_p .* xs_r;
       color=:royalblue, linewidth=2,
       label="Prior  y=$(round(a_p,digits=1)) + $(round(b_p,digits=2))x")

# Posterior
scatter!(ax3, y_obs, post_mean; color=(:firebrick, 0.35), markersize=4)
a_q, b_q = linreg(y_obs, post_mean)
lines!(ax3, xs_r, a_q .+ b_q .* xs_r;
       color=:firebrick, linewidth=2,
       label="Posterior  y=$(round(a_q,digits=1)) + $(round(b_q,digits=2))x")

xlims!(ax3, xl); ylims!(ax3, xl)
Legend(fig3[2, 1], ax3; orientation=:horizontal, framevisible=false, tellwidth=false)

oto_path = joinpath(OUTDIR, "plot_1to1.png")
save(oto_path, fig3; px_per_unit=2)
println("  Saved → $oto_path")

println("\nAll done.")
