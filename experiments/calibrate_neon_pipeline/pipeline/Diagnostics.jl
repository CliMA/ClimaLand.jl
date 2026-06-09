"""
Diagnostics.jl — `plot_eki_diagnostics(run; output_dir, eki_path, obs_filepath)`
produces the EKI diagnostic figures and returns the figures dir + final RMSE.

Function form of the original plot_eki_diagnostics.jl (frozen in
../calibrate_neon). All configuration is passed as arguments — no ENV, no
globals. Reads everything from eki_file.jld2 (y_obs + per-iteration G matrices)
and observations.jld2 (obs dates).

Returns `(; figures_dir, final_rmse)` where `final_rmse` is the posterior
ensemble-mean RMSE (fed into the master CSV).

This file is `include`d into Main by the driver.
"""

using CairoMakie
using Dates
using Statistics
using LinearAlgebra
import JLD2
import EnsembleKalmanProcesses as EKP

# ── EKP file readers (positional field access, as in the original) ───────────
"""
Extract y_obs and the per-iteration G matrices from an EKP JLD2 file.
Uses the EKP public API (get_obs / get_g) rather than positional field access,
which is robust to EKP struct changes.
"""
function _load_eki_data(eki_path)
    ekp = JLD2.load_object(eki_path)
    y_obs = Vector{Float64}(EKP.get_obs(ekp))
    # get_g(ekp) returns the (n_obs × N_ens) forward map for each completed
    # iteration as a vector of matrices.
    G = [Matrix{Float64}(g) for g in EKP.get_g(ekp; return_array = true)]
    return y_obs, G
end

function _linreg(x, y)
    mask = .!isnan.(x) .& .!isnan.(y)
    xm, ym = x[mask], y[mask]
    n = length(xm)
    sx = sum(xm); sy = sum(ym)
    sxx = dot(xm, xm); sxy = dot(xm, ym)
    b = (n * sxy - sx * sy) / (n * sxx - sx^2)
    a = (sy - b * sx) / n
    return a, b
end

function _monthly_mean(dates, values)
    by_month = [Float64[] for _ in 1:12]
    for (d, v) in zip(dates, values)
        isnan(v) || push!(by_month[Dates.month(d)], v)
    end
    means = [isempty(v) ? NaN : mean(v) for v in by_month]
    stds = [length(v) < 2 ? NaN : std(v) for v in by_month]
    return means, stds
end

function _ensemble_rmse(G_mat, y)
    N_ens = size(G_mat, 2)
    [let d = G_mat[:, m] .- y; ok = .!isnan.(d); sqrt(mean(d[ok] .^ 2)) end
     for m in 1:N_ens]
end

"""
    plot_eki_diagnostics(run; output_dir, eki_path, obs_filepath) -> (; figures_dir, final_rmse)
"""
function plot_eki_diagnostics(run; output_dir, eki_path, obs_filepath)
    site_id = run.site
    cal_depth_str = string(run.cal_depth)
    outdir = joinpath(output_dir, "figures_eki_diagnostics")
    mkpath(outdir)

    println("Loading EKP data from $eki_path …")
    y_obs, G_all = _load_eki_data(eki_path)
    n_iters = length(G_all)
    n_obs = length(y_obs)
    println("  y_obs length: $n_obs   completed iterations: $n_iters")

    obs_data = JLD2.load(obs_filepath)
    obs_dates = obs_data["obs_dates"]
    @assert length(obs_dates) == n_obs "Date count mismatch: $(length(obs_dates)) vs $n_obs"

    G_prior = G_all[1]
    G_post = G_all[end]
    prior_mean = vec(mean(G_prior, dims = 2))
    post_mean = vec(mean(G_post, dims = 2))
    iter_label_post = n_iters - 1
    months_label = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]

    # Posterior ensemble-mean RMSE (computed early so it can go in titles + file).
    final_rmse = mean(_ensemble_rmse(G_post, y_obs))
    rmse_str = round(final_rmse, digits = 2)

    # ── PLOT 1 — timeseries prior/post/obs ───────────────────────────────────
    fig1 = Figure(size = (1200, 500))
    ax1 = Axis(fig1[1, 1]; xlabel = "Day index", ylabel = "Soil CO₂ (ppm)",
        title = "$site_id: Daily Soil CO₂ at ~$cal_depth_str — Prior vs Posterior vs Obs")
    lines!(ax1, 1:n_obs, y_obs; color = :black, linewidth = 2.5, label = "NEON Observations")
    lines!(ax1, 1:n_obs, prior_mean; color = :royalblue, linewidth = 2, label = "Prior (iter 0)")
    lines!(ax1, 1:n_obs, post_mean; color = :firebrick, linewidth = 2,
        label = "Posterior (iter $iter_label_post)")
    Legend(fig1[1, 2], ax1; framevisible = false)
    save(joinpath(outdir, "plot_timeseries.png"), fig1; px_per_unit = 2)

    # ── PLOT 1a — posterior only ─────────────────────────────────────────────
    fig1a = Figure(size = (1200, 500))
    ax1a = Axis(fig1a[1, 1]; xlabel = "Day index", ylabel = "Soil CO₂ (ppm)",
        title = "$site_id: Daily Soil CO₂ — Posterior vs Obs (RMSE = $rmse_str ppm)")
    lines!(ax1a, 1:n_obs, y_obs; color = :black, linewidth = 2.5, label = "NEON Observations")
    lines!(ax1a, 1:n_obs, post_mean; color = :firebrick, linewidth = 2,
        label = "Posterior (iter $iter_label_post)")
    Legend(fig1a[1, 2], ax1a; framevisible = false)
    save(joinpath(outdir, "plot_timeseries_post.png"), fig1a; px_per_unit = 2)

    # ── PLOT 1b — monthly seasonality ────────────────────────────────────────
    fig1b = Figure(size = (800, 500))
    ax1b = Axis(fig1b[1, 1]; xlabel = "Month", ylabel = "Soil CO₂ (ppm)",
        title = "$site_id: Monthly Mean Soil CO₂ Seasonality",
        xticks = (1:12, months_label))
    obs_m, obs_s = _monthly_mean(obs_dates, y_obs)
    prior_m, prior_s = _monthly_mean(obs_dates, prior_mean)
    post_m, post_s = _monthly_mean(obs_dates, post_mean)
    xs = 1:12
    band!(ax1b, xs, obs_m .- obs_s, obs_m .+ obs_s; color = (:black, 0.12))
    band!(ax1b, xs, prior_m .- prior_s, prior_m .+ prior_s; color = (:royalblue, 0.15))
    band!(ax1b, xs, post_m .- post_s, post_m .+ post_s; color = (:firebrick, 0.15))
    lines!(ax1b, xs, obs_m; color = :black, linewidth = 2.5, label = "Observations")
    lines!(ax1b, xs, prior_m; color = :royalblue, linewidth = 2, label = "Prior (iter 0)")
    lines!(ax1b, xs, post_m; color = :firebrick, linewidth = 2,
        label = "Posterior (iter $iter_label_post)")
    scatter!(ax1b, xs, obs_m; color = :black, markersize = 8)
    scatter!(ax1b, xs, prior_m; color = :royalblue, markersize = 6)
    scatter!(ax1b, xs, post_m; color = :firebrick, markersize = 6)
    Legend(fig1b[1, 2], ax1b; framevisible = false)
    save(joinpath(outdir, "plot_seasonality.png"), fig1b; px_per_unit = 2)

    # ── PLOT 2 — RMSE per iteration ──────────────────────────────────────────
    all_rmse = [_ensemble_rmse(G_all[i], y_obs) for i in 1:n_iters]
    for (i, r) in enumerate(all_rmse)
        println("  Iter $(i-1): mean RMSE=$(round(mean(r), digits=2)), std=$(round(std(r), digits=2))")
    end
    fig2 = Figure(size = (max(600, 300 + 200 * n_iters), 500))
    ax2 = Axis(fig2[1, 1]; xlabel = "EKI Iteration", ylabel = "RMSE (Soil CO₂ ppm)",
        title = "$site_id: Ensemble RMSE by Iteration",
        xticks = (0:n_iters-1, string.(0:n_iters-1)))
    colors_iter = [i == 1 ? :royalblue : i == n_iters ? :firebrick : :gray for i in 1:n_iters]
    for (i, rmse_v) in enumerate(all_rmse)
        xi = i - 1
        col = colors_iter[i]
        jitter = randn(Float32, length(rmse_v)) .* 0.05f0
        scatter!(ax2, fill(xi, length(rmse_v)) .+ jitter, rmse_v;
            color = (col, 0.55), markersize = 7, label = i == 1 ? "Ensemble members" : "")
        μ = mean(rmse_v); σ = std(rmse_v)
        errorbars!(ax2, [xi], [μ], [σ]; color = col, linewidth = 2, whiskerwidth = 12)
        scatter!(ax2, [xi], [μ]; color = col, markersize = 14, marker = :diamond,
            label = i == 1 ? "Prior mean (iter 0)" :
                    i == n_iters ? "Posterior mean (iter $(i-1))" : "Mean (iter $(i-1))")
    end
    Legend(fig2[1, 2], ax2; framevisible = false)
    save(joinpath(outdir, "plot_rmse_per_iteration.png"), fig2; px_per_unit = 2)

    # ── PLOT 3 — 1-to-1 ──────────────────────────────────────────────────────
    fig3 = Figure(size = (600, 500))
    ax3 = Axis(fig3[1, 1]; xlabel = "Observed Soil CO₂ (ppm)",
        ylabel = "Modelled Soil CO₂ (ppm)", title = "$site_id: Soil CO₂ 1-to-1",
        aspect = DataAspect())
    vmin = min(minimum(filter(!isnan, y_obs)),
               minimum(filter(!isnan, vcat(prior_mean, post_mean))))
    vmax = max(maximum(filter(!isnan, y_obs)),
               maximum(filter(!isnan, vcat(prior_mean, post_mean))))
    pad = 0.06 * (vmax - vmin)
    xl = (vmin - pad, vmax + pad)
    lines!(ax3, collect(xl), collect(xl); color = :black, linewidth = 1.5,
        linestyle = :dash, label = "1:1")
    scatter!(ax3, y_obs, prior_mean; color = (:royalblue, 0.35), markersize = 4)
    a_p, b_p = _linreg(y_obs, prior_mean)
    xs_r = range(xl[1], xl[2]; length = 200)
    lines!(ax3, xs_r, a_p .+ b_p .* xs_r; color = :royalblue, linewidth = 2,
        label = "Prior  y=$(round(a_p,digits=1)) + $(round(b_p,digits=2))x")
    scatter!(ax3, y_obs, post_mean; color = (:firebrick, 0.35), markersize = 4)
    a_q, b_q = _linreg(y_obs, post_mean)
    lines!(ax3, xs_r, a_q .+ b_q .* xs_r; color = :firebrick, linewidth = 2,
        label = "Posterior  y=$(round(a_q,digits=1)) + $(round(b_q,digits=2))x")
    xlims!(ax3, xl); ylims!(ax3, xl)
    Legend(fig3[2, 1], ax3; orientation = :horizontal, framevisible = false, tellwidth = false)
    save(joinpath(outdir, "plot_1to1.png"), fig3; px_per_unit = 2)

    # Append the posterior RMSE to final_parameter_means.txt (written by the
    # calibration step; the RMSE is only known here, after the diagnostics).
    final_params_file = joinpath(output_dir, "final_parameter_means.txt")
    if isfile(final_params_file)
        contents = read(final_params_file, String)
        if !occursin("final_rmse", contents)
            open(final_params_file, "a") do io
                println(io, "final_rmse = $(round(final_rmse, digits = 3))  # posterior ensemble-mean RMSE (ppm)")
            end
        end
    end

    println("  Posterior ensemble-mean RMSE = $(round(final_rmse, digits=3)) ppm")
    println("All EKI diagnostic figures saved to $outdir")
    return (; figures_dir = outdir, final_rmse = final_rmse)
end
