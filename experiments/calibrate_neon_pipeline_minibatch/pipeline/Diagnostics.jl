"""
Diagnostics.jl — `plot_eki_diagnostics(run; output_dir, eki_path, obs_filepath)`
produces the EKI RMSE-per-iteration figure and returns the figures dir + final
RMSE.

MINIBATCH: under minibatching each EKI iteration sees a DIFFERENT subset of
windows, so there is no single flat `y_obs` aligned across iterations. The old
per-depth timeseries / seasonality / 1-to-1 plots (which assumed one fixed obs
vector) are therefore dropped. Only the RMSE-per-iteration plot is kept: for each
iteration `i`, the ensemble G is compared against THAT iteration's own minibatch
observation, `EKP.get_obs(ekp, i)` — a pure, cursor-free lookup. RMSE points may
step between iterations when different windows are drawn; that is expected.

Returns `(; figures_dir, final_rmse)` where `final_rmse` is the posterior
ensemble-mean RMSE on the FINAL iteration's minibatch (fed into the master CSV).

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
Extract the EKP object and the per-iteration G matrices from an EKP JLD2 file.
Uses the EKP public API (get_g) rather than positional field access. MINIBATCH:
we return the `ekp` itself (not a single y_obs) so the caller can fetch each
iteration's own minibatch observation via `EKP.get_obs(ekp, i)`.
"""
function _load_eki_data(eki_path)
    ekp = JLD2.load_object(eki_path)
    # get_g(ekp) returns the (n_obs × N_ens) forward map for each completed
    # iteration as a vector of matrices (n_obs = that iteration's minibatch length).
    G = [Matrix{Float64}(g) for g in EKP.get_g(ekp; return_array = true)]
    return ekp, G
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

MINIBATCH: keeps only the RMSE-per-iteration plot. Each iteration's ensemble G is
scored against THAT iteration's own minibatch observation (`EKP.get_obs(ekp, i)`),
so G and y always have matching length even though the drawn windows change from
iteration to iteration. `final_rmse` is the posterior ensemble-mean RMSE on the
last iteration's minibatch.
"""
function plot_eki_diagnostics(run; output_dir, eki_path, obs_filepath)
    site_id = run.site
    outdir = joinpath(output_dir, "figures_eki_diagnostics")
    mkpath(outdir)

    println("Loading EKP data from $eki_path …")
    ekp, G_all = _load_eki_data(eki_path)
    n_iters = length(G_all)
    println("  completed iterations: $n_iters")

    # Per-iteration ensemble RMSE, each against its OWN minibatch observation.
    # get_obs(ekp, i) is a pure, cursor-free lookup into the ObservationSeries.
    all_rmse = Vector{Vector{Float64}}(undef, n_iters)
    for i in 1:n_iters
        y_i = Vector{Float64}(EKP.get_obs(ekp, i))
        if length(y_i) != size(G_all[i], 1)
            @warn "iter $i: minibatch obs length $(length(y_i)) ≠ G rows " *
                  "$(size(G_all[i], 1)); skipping (NaN)."
            all_rmse[i] = fill(NaN, size(G_all[i], 2))
        else
            all_rmse[i] = _ensemble_rmse(G_all[i], y_i)
        end
    end
    for (i, r) in enumerate(all_rmse)
        println("  Iter $(i-1): mean RMSE=$(round(mean(filter(!isnan, r)), digits=2)), ",
            "std=$(round(std(filter(!isnan, r)), digits=2))")
    end

    # Posterior ensemble-mean RMSE = final iteration's mean (fed into the CSV).
    final_rmse = mean(filter(!isnan, all_rmse[end]))

    # ── RMSE per iteration ───────────────────────────────────────────────────
    fig2 = Figure(size = (max(600, 300 + 200 * n_iters), 500))
    ax2 = Axis(fig2[1, 1]; xlabel = "EKI Iteration",
        ylabel = "RMSE (Soil CO₂ ppm, per-iter minibatch)",
        title = "$site_id: Ensemble RMSE by Iteration (minibatch)",
        xticks = (0:n_iters-1, string.(0:n_iters-1)))
    colors_iter = [i == 1 ? :royalblue : i == n_iters ? :firebrick : :gray for i in 1:n_iters]
    for (i, rmse_v) in enumerate(all_rmse)
        xi = i - 1
        col = colors_iter[i]
        jitter = randn(Float32, length(rmse_v)) .* 0.05f0
        scatter!(ax2, fill(xi, length(rmse_v)) .+ jitter, rmse_v;
            color = (col, 0.55), markersize = 7, label = i == 1 ? "Ensemble members" : "")
        μ = mean(filter(!isnan, rmse_v)); σ = std(filter(!isnan, rmse_v))
        errorbars!(ax2, [xi], [μ], [σ]; color = col, linewidth = 2, whiskerwidth = 12)
        scatter!(ax2, [xi], [μ]; color = col, markersize = 14, marker = :diamond,
            label = i == 1 ? "Prior mean (iter 0)" :
                    i == n_iters ? "Posterior mean (iter $(i-1))" : "Mean (iter $(i-1))")
    end
    Legend(fig2[1, 2], ax2; framevisible = false)
    save(joinpath(outdir, "plot_rmse_per_iteration.png"), fig2; px_per_unit = 2)

    # Append the posterior RMSE to final_parameter_means.txt (written by the
    # calibration step; the RMSE is only known here, after the diagnostics).
    final_params_file = joinpath(output_dir, "final_parameter_means.txt")
    if isfile(final_params_file)
        contents = read(final_params_file, String)
        if !occursin("final_rmse", contents)
            open(final_params_file, "a") do io
                println(io, "final_rmse = $(round(final_rmse, digits = 3))  # posterior ensemble-mean RMSE on final minibatch (ppm)")
            end
        end
    end

    println("  Posterior ensemble-mean RMSE (final minibatch) = $(round(final_rmse, digits=3)) ppm")
    println("EKI RMSE figure saved to $outdir")
    return (; figures_dir = outdir, final_rmse = final_rmse)
end
