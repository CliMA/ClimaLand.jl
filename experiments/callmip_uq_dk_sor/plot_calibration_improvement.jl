"""
Plot calibration improvement: RMSE convergence + monthly NEE comparison.
"""

import JLD2, Statistics, Dates, NCDatasets
using Plots, Printf
gr()

cal_dir = joinpath(@__DIR__, "..", "calibrate_dk_sor")
obs_nc  = joinpath(@__DIR__, "..", "..", "DK_Sor",
    "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")
out_dir = joinpath(@__DIR__, "output_evaluation")
isdir(out_dir) || mkdir(out_dir)

println("Loading observations..."); flush(stdout)
ds = NCDatasets.NCDataset(obs_nc)
nee_obs_v  = Float64.(NCDatasets.coalesce.(ds["NEE_daily"][:], NaN))
qle_obs_v  = Float64.(NCDatasets.coalesce.(ds["Qle_daily"][:], NaN))
qh_obs_v   = Float64.(NCDatasets.coalesce.(ds["Qh_daily"][:], NaN))
obs_dates  = Dates.Date.(ds["time"][:])
NCDatasets.close(ds)
obs_nee = Dict(zip(obs_dates, nee_obs_v))
obs_qle = Dict(zip(obs_dates, qle_obs_v))
obs_qh  = Dict(zip(obs_dates, qh_obs_v ))

function load_diag(iter, member)
    p = joinpath(cal_dir, "output",
        "iteration_$(lpad(iter,3,'0'))",
        "member_$(lpad(member,3,'0'))",
        "daily_diagnostics.jld2")
    isfile(p) || return nothing
    d = JLD2.load(p)
    (dates = Dates.Date.(d["dates"]),
     nee   = d["nee"] .* 12.0 .* 86400.0,
     qle   = d["qle"],
     qh    = d["qh"])
end

function rmse_flux(dates, vals, obs)
    e2 = [(vals[i] - obs[dates[i]])^2
           for i in eachindex(dates)
           if haskey(obs, dates[i]) && !isnan(obs[dates[i]])]
    isempty(e2) ? NaN : sqrt(Statistics.mean(e2))
end

function monthly_mean(dates, vals)
    [Statistics.mean(filter(!isnan,
        [vals[i] for i in eachindex(dates) if Dates.month(dates[i]) == mo]))
     for mo in 1:12]
end

# ── Convergence table ─────────────────────────────────────────────────────────
println("Computing RMSE convergence..."); flush(stdout)
iters   = Int[]
nee_rmse = Float64[]
qle_rmse = Float64[]
qh_rmse  = Float64[]

for iter in 0:9
    nr = Float64[]; qr = Float64[]; hr = Float64[]
    for m in 1:33
        r = load_diag(iter, m); isnothing(r) && continue
        v = rmse_flux(r.dates, r.nee, obs_nee); isnan(v) || push!(nr, v)
        v = rmse_flux(r.dates, r.qle, obs_qle); isnan(v) || push!(qr, v)
        v = rmse_flux(r.dates, r.qh,  obs_qh);  isnan(v) || push!(hr, v)
    end
    isempty(nr) && continue
    push!(iters,    iter)
    push!(nee_rmse, Statistics.mean(nr))
    push!(qle_rmse, Statistics.mean(qr))
    push!(qh_rmse,  Statistics.mean(hr))
end

# ── Best posterior member ─────────────────────────────────────────────────────
function find_best(obs_nee)
    br = Inf; bm = 1; bi = 0
    for iter in 0:9, m in 1:33
        r = load_diag(iter, m); isnothing(r) && continue
        v = rmse_flux(r.dates, r.nee, obs_nee)
        if v < br; br = v; bm = m; bi = iter; end
    end
    br, bm, bi
end
best_rmse, best_m, best_i = find_best(obs_nee)
r_best  = load_diag(best_i, best_m)
r_prior = load_diag(0, 1)

# ── Monthly means ─────────────────────────────────────────────────────────────
obs_mn     = [Statistics.mean(filter(!isnan,
                  [get(obs_nee, d, NaN) for d in obs_dates if Dates.month(d) == mo]))
              for mo in 1:12]
prior_mn   = monthly_mean(r_prior.dates, r_prior.nee)
post_mn    = monthly_mean(r_best.dates,  r_best.nee)
months     = 1:12
mn_labels  = ["J","F","M","A","M","J","J","A","S","O","N","D"]

# ─── Plot ─────────────────────────────────────────────────────────────────────
println("Generating figure..."); flush(stdout)

# Panel 1 — RMSE convergence (log scale, NEE only since it changes dramatically)
p1 = plot(iters, nee_rmse;
    yscale    = :log10,
    marker    = :circle, ms = 6,
    lw        = 2.5, color = :steelblue,
    xlabel    = "EKI iteration",
    ylabel    = "NEE RMSE (gC m⁻² d⁻¹)",
    title     = "EKI Convergence — NEE",
    legend    = false,
    grid      = true, gridalpha = 0.3,
    xticks    = iters,
    ylims     = (1.5, 1.5e4),
    yticks    = ([10^0, 10^1, 10^2, 10^3, 10^4],
                 ["1","10","100","1000","10 000"]),
    titlefontsize = 11)
# Annotate starting and ending values
annotate!(p1, iters[1]   + 0.2, nee_rmse[1]   * 2.5,
    text(@sprintf("%.0f", nee_rmse[1]),   8, :left, :darkgray))
annotate!(p1, iters[end] + 0.2, nee_rmse[end] * 2.5,
    text(@sprintf("%.2f", nee_rmse[end]), 8, :left, :steelblue))

# Panel 2 — Qle / Qh convergence (linear)
p2 = plot(iters, qle_rmse;
    marker = :circle, ms = 5, lw = 2, color = :darkorange,
    label  = "Qle", xlabel = "EKI iteration",
    ylabel = "RMSE (W m⁻²)", title = "EKI Convergence — Qle & Qh",
    legend = :topright, grid = true, gridalpha = 0.3,
    xticks = iters, titlefontsize = 11)
plot!(p2, iters, qh_rmse;
    marker = :square, ms = 5, lw = 2, color = :crimson, label = "Qh")

# Panel 3 — Monthly NEE seasonal cycle
p3 = plot(months, prior_mn;
    color = :gray, lw = 2, ls = :dash, marker = :circle, ms = 4,
    label = "Prior",
    xlabel = "Month", ylabel = "NEE (gC m⁻² d⁻¹)",
    title  = "Mean Seasonal Cycle — NEE",
    xticks = (months, mn_labels), legend = :bottomleft,
    grid = true, gridalpha = 0.3, titlefontsize = 11)
plot!(p3, months, post_mn;
    color = :steelblue, lw = 2.5, marker = :circle, ms = 4,
    label = @sprintf("Posterior (iter %d, m%03d)", best_i, best_m))
plot!(p3, months, obs_mn;
    color = :black, lw = 2, ls = :dashdot, marker = :diamond, ms = 4,
    label = "FLUXNET obs")
hline!(p3, [0.0]; color = :black, lw = 0.8, ls = :dot, label = "")

# Assemble
pfig = plot(p1, p2, p3;
    layout  = (1, 3),
    size    = (1400, 420),
    leftmargin  = 5Plots.mm,
    bottommargin = 8Plots.mm,
    plot_title  = "ClimaLand DK-Sor — EKI Calibration Improvement",
    plot_titlefontsize = 12)

out_path = joinpath(out_dir, "calibration_improvement.png")
savefig(pfig, out_path)
println("Saved → $out_path")

# ── Text summary ──────────────────────────────────────────────────────────────
println()
println("="^60)
println("  Calibration Summary — DK-Sor (ClimaLand Phase 1a)")
println("="^60)
@printf("  NEE RMSE  : %.2f → %.2f  gC/m²/d  (−%.0f%%)\n",
    nee_rmse[1], nee_rmse[end],
    100*(nee_rmse[1] - nee_rmse[end]) / nee_rmse[1])
@printf("  EKI best  : iter %d, member %03d, NEE RMSE = %.4f\n",
    best_i, best_m, best_rmse)
@printf("  Qle RMSE  : %.1f → %.1f  W/m²  (Δ%.1f)\n",
    qle_rmse[1], qle_rmse[end], qle_rmse[end] - qle_rmse[1])
@printf("  Qh  RMSE  : %.1f → %.1f  W/m²  (Δ%.1f)\n",
    qh_rmse[1],  qh_rmse[end],  qh_rmse[end]  - qh_rmse[1])
println("="^60)
println("  Seasonal NEE (gC/m²/d): Prior / Post / Obs")
for mo in 1:12
    @printf("    %-3s  %+6.2f  %+6.2f  %+6.2f\n",
        Dates.monthabbr(Dates.Date(2000, mo, 1)),
        prior_mn[mo], post_mn[mo], obs_mn[mo])
end
println("="^60)
println("Done.")
