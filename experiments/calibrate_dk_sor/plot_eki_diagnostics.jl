"""
Diagnostic plots for DK-Sor EKI calibration results.

Reads eki_file.jld2 (at project root) for all observation and model data:
  - y_obs  — observation vector embedded in the EKP ObservationSeries
  - G      — G-matrices from each completed iteration

Date metadata (obs_dates) is reconstructed from the same NetCDF filter used
when building the observations, but only the dates are read — all observation
values come from the EKP file.

Three output figures:
  1. plot_seasonality.png        — monthly mean obs / prior / posterior
  2. plot_rmse_per_iteration.png — ensemble RMSE per iteration
  3. plot_1to1.png               — 1-to-1 scatter with regression lines

Usage:
    julia --project=.buildkite experiments/calibrate_dk_sor/plot_eki_diagnostics.jl
"""

using CairoMakie
using NCDatasets
using Dates
using Statistics
using LinearAlgebra
import JLD2
import ClimaLand

# ── Paths ─────────────────────────────────────────────────────────────────────

const climaland_dir = pkgdir(ClimaLand)
const EKI_PATH = joinpath(climaland_dir, "eki_file.jld2")
const FLUX_NC  = joinpath(climaland_dir, "DK_Sor",
                          "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")
const MET_NC   = joinpath(climaland_dir, "DK_Sor",
                          "DK-Sor_1997-2014_FLUXNET2015_Met.nc")
const OUTDIR   = joinpath(climaland_dir, "experiments/calibrate_dk_sor")

# ── Load all data from eki_file.jld2 ─────────────────────────────────────────

"""
Extract y_obs vector and G-matrices from the EKP's JLD2 file.

EKP field layout (positional):
  1 = u,  2 = observation_series,  3 = N_ens,  4 = g, ...

ObservationSeries.observations[1].samples[1] is the y vector.
G[i].stored_data is the (n_obs_total × N_ens) output matrix.
"""
function load_eki_data()
    f   = JLD2.jldopen(EKI_PATH, "r")
    obj = f["single_stored_object"]
    fld = getfield(obj, :fields)          # positional vector of EKP fields

    # ── y_obs from observation_series ────────────────────────────────────────
    obs_series = fld[2]
    os_fld     = getfield(obs_series, :fields)   # (:observations, :minibatcher, ...)
    obs_vec    = os_fld[1]                        # Vector of Observation objects
    ob1_fld    = getfield(obs_vec[1], :fields)    # (:samples, :covs, ...)
    y_obs      = ob1_fld[1][1] :: Vector{Float64} # samples[1]

    # ── G-matrices ────────────────────────────────────────────────────────────
    g_containers = fld[4]
    n_iters      = length(g_containers)
    G = [getfield(g_containers[i], :fields)[1] :: Matrix{Float64}
         for i in 1:n_iters]             # each (n_obs_total, N_ens)

    JLD2.close(f)
    return y_obs, G
end

# ── Obs dates: read only date/wind columns from NetCDF ────────────────────────

"""
Rebuild the ordered list of valid observation dates using the same mask
that was applied when generating the observation vector.  Only dates and
wind speed are read here; no flux values are loaded.
"""
function load_obs_dates(n_obs_expected::Int)
    flux_ds    = NCDataset(FLUX_NC, "r")
    flux_times = Date.(flux_ds["time"][:])
    # Need to know which days have valid flux data — re-read the validity flags
    nee_raw    = Float64.(coalesce.(flux_ds["NEE_daily"][:], NaN))
    qle_raw    = Float64.(coalesce.(flux_ds["Qle_daily"][:], NaN))
    qh_raw     = Float64.(coalesce.(flux_ds["Qh_daily"][:], NaN))
    close(flux_ds)

    met_ds     = NCDataset(MET_NC, "r")
    wind_raw   = Float64.(coalesce.(met_ds["Wind"][1, 1, :], NaN))
    met_times  = Date.(met_ds["time"][:])
    close(met_ds)

    wind_by_day = Dict{Date, Vector{Float64}}()
    for (i, d) in enumerate(met_times)
        isnan(wind_raw[i]) && continue
        push!(get!(wind_by_day, d, Float64[]), wind_raw[i])
    end
    daily_wind = Dict(d => mean(v) for (d, v) in wind_by_day)

    cal_mask   = Date(2004,1,1) .<= flux_times .<= Date(2013,12,31)
    valid_mask = cal_mask .&
        .!isnan.(nee_raw) .& .!isnan.(qle_raw) .& .!isnan.(qh_raw) .&
        (abs.(nee_raw) .< 1e10) .& (abs.(qle_raw) .< 1e10) .& (abs.(qh_raw) .< 1e10)
    for i in eachindex(valid_mask)
        if valid_mask[i]
            w = get(daily_wind, flux_times[i], NaN)
            valid_mask[i] = !(isnan(w) || w >= 5.0)
        end
    end

    obs_dates = flux_times[valid_mask]
    @assert length(obs_dates) == n_obs_expected """
        Date filter returned $(length(obs_dates)) days but EKP y_obs implies \
        $n_obs_expected.  Check that generate_observations.jl settings match.
    """
    return obs_dates
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
n_obs        = length(y_obs) ÷ 3
println("  y_obs length: $(length(y_obs)) → $n_obs observation days")
println("  Completed iterations: $n_iters")

println("Loading observation dates…")
obs_dates = load_obs_dates(n_obs)

# Split into three variables
y_nee = y_obs[1:n_obs]
y_qle = y_obs[n_obs+1:2*n_obs]
y_qh  = y_obs[2*n_obs+1:end]

G_prior = G_all[1]          # forward run from prior ensemble
G_post  = G_all[end]        # forward run from last updated ensemble

prior_nee_mean = vec(mean(G_prior[1:n_obs,         :], dims=2))
prior_qle_mean = vec(mean(G_prior[n_obs+1:2*n_obs, :], dims=2))
prior_qh_mean  = vec(mean(G_prior[2*n_obs+1:end,   :], dims=2))

post_nee_mean  = vec(mean(G_post[1:n_obs,         :], dims=2))
post_qle_mean  = vec(mean(G_post[n_obs+1:2*n_obs, :], dims=2))
post_qh_mean   = vec(mean(G_post[2*n_obs+1:end,   :], dims=2))

iter_label_post = n_iters - 1   # 0-based: G_all[1] is iter 0, G_all[end] is iter n_iters-1

months_label = ["Jan","Feb","Mar","Apr","May","Jun",
                "Jul","Aug","Sep","Oct","Nov","Dec"]

# ─────────────────────────────────────────────────────────────────────────────
#  PLOT 1 — Monthly seasonality
# ─────────────────────────────────────────────────────────────────────────────

println("Plotting seasonality…")

fig1 = Figure(size=(1200, 850))

var_rows = [
    (y_nee, prior_nee_mean, post_nee_mean, "NEE (gC m⁻² d⁻¹)"),
    (y_qle, prior_qle_mean, post_qle_mean, "Qle (W m⁻²)"),
    (y_qh,  prior_qh_mean,  post_qh_mean,  "Qh (W m⁻²)"),
]

for (row, (y_v, gp_v, gpost_v, ylabel)) in enumerate(var_rows)
    ax = Axis(fig1[row, 1];
              xlabel  = row == 3 ? "Month" : "",
              ylabel  = ylabel,
              xticks  = (1:12, months_label),
              title   = row == 1 ? "DK-Sor Monthly Seasonality (2004–2013)" : "")

    obs_m,   obs_s   = monthly_mean(obs_dates, y_v)
    prior_m, prior_s = monthly_mean(obs_dates, gp_v)
    post_m,  post_s  = monthly_mean(obs_dates, gpost_v)
    xs = 1:12

    band!(ax, xs, obs_m   .- obs_s,   obs_m   .+ obs_s;   color=(:black,      0.12))
    band!(ax, xs, prior_m .- prior_s, prior_m .+ prior_s; color=(:royalblue,  0.15))
    band!(ax, xs, post_m  .- post_s,  post_m  .+ post_s;  color=(:firebrick,  0.15))

    lines!(ax, xs, obs_m;   color=:black,      linewidth=2.5, label="Observations")
    lines!(ax, xs, prior_m; color=:royalblue,  linewidth=2,   label="Prior (iter 0)")
    lines!(ax, xs, post_m;  color=:firebrick,  linewidth=2,
           label="Posterior (iter $iter_label_post)")

    scatter!(ax, xs, obs_m;   color=:black,     markersize=8)
    scatter!(ax, xs, prior_m; color=:royalblue, markersize=6)
    scatter!(ax, xs, post_m;  color=:firebrick, markersize=6)

    row == 1 && Legend(fig1[row, 2], ax; framevisible=false)
end

seasonality_path = joinpath(OUTDIR, "plot_seasonality.png")
save(seasonality_path, fig1; px_per_unit=2)
println("  Saved → $seasonality_path")

# ─────────────────────────────────────────────────────────────────────────────
#  PLOT 2 — RMSE per iteration
# ─────────────────────────────────────────────────────────────────────────────

println("Plotting RMSE per iteration…")

# Compute RMSE for every G matrix
all_rmse = [ensemble_rmse(G_all[i], y_obs) for i in 1:n_iters]
for (i, r) in enumerate(all_rmse)
    println("  Iter $(i-1): mean RMSE=$(round(mean(r), digits=2)), std=$(round(std(r), digits=2))")
end

fig2 = Figure(size=(max(600, 300 + 200*n_iters), 500))
ax2  = Axis(fig2[1, 1];
            xlabel = "EKI Iteration",
            ylabel = "RMSE (combined NEE + Qle + Qh)",
            title  = "DK-Sor: Ensemble RMSE by Iteration",
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

fig3 = Figure(size=(1200, 480))

var_info = [
    (y_nee, prior_nee_mean, post_nee_mean, "NEE (gC m⁻² d⁻¹)"),
    (y_qle, prior_qle_mean, post_qle_mean, "Qle (W m⁻²)"),
    (y_qh,  prior_qh_mean,  post_qh_mean,  "Qh (W m⁻²)"),
]

for (col, (y_v, gp_v, gpost_v, varname)) in enumerate(var_info)
    ax = Axis(fig3[1, col];
              xlabel  = "Observed $varname",
              ylabel  = col == 1 ? "Modelled $varname" : "",
              title   = varname,
              aspect  = DataAspect())

    vmin = min(minimum(filter(!isnan, y_v)),
               minimum(filter(!isnan, vcat(gp_v, gpost_v))))
    vmax = max(maximum(filter(!isnan, y_v)),
               maximum(filter(!isnan, vcat(gp_v, gpost_v))))
    pad  = 0.06 * (vmax - vmin)
    xl   = (vmin - pad, vmax + pad)
    lines!(ax, collect(xl), collect(xl);
           color=:black, linewidth=1.5, linestyle=:dash, label="1:1")

    # Prior
    scatter!(ax, y_v, gp_v; color=(:royalblue, 0.35), markersize=4)
    a_p, b_p = linreg(y_v, gp_v)
    xs_r     = range(xl[1], xl[2]; length=200)
    lines!(ax, xs_r, a_p .+ b_p .* xs_r;
           color=:royalblue, linewidth=2,
           label="Prior  y=$(round(a_p,digits=1)) + $(round(b_p,digits=2))x")

    # Posterior
    scatter!(ax, y_v, gpost_v; color=(:firebrick, 0.35), markersize=4)
    a_q, b_q = linreg(y_v, gpost_v)
    lines!(ax, xs_r, a_q .+ b_q .* xs_r;
           color=:firebrick, linewidth=2,
           label="Posterior  y=$(round(a_q,digits=1)) + $(round(b_q,digits=2))x")

    xlims!(ax, xl); ylims!(ax, xl)
    col == 1 && Legend(fig3[2, 1:3], ax;
                       orientation=:horizontal, framevisible=false,
                       tellwidth=false)
end

oto_path = joinpath(OUTDIR, "plot_1to1.png")
save(oto_path, fig3; px_per_unit=2)
println("  Saved → $oto_path")

println("\nAll done.")
