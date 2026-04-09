"""
Multi-year flux check: plot NEE, Qle, Qh for 2004-2013 (model vs obs).
Uses the posterior UKI iteration to check for drift and reasonableness.
Run with:
    julia --project=.buildkite experiments/calibrate_dk_sor/plot_multiyear_check.jl
"""

using Statistics, Printf, Dates
import JLD2, CairoMakie

const climaland_dir = abspath(joinpath(@__DIR__, "..", ".."))

# ── Load calibration diagnostics (use median member from final iteration) ──────
# Member 016 = central/mean sigma-point in 33-member UKI ensemble
iter_dir = joinpath(@__DIR__, "output", "iteration_009")
member   = "member_016"
diag = JLD2.load(joinpath(iter_dir, member, "daily_diagnostics.jld2"))

dates_all = diag["dates"]
nee_all   = diag["nee"] .* (12 * 86400)   # mol/m²/s → gC/m²/d
qle_all   = diag["qle"]                    # W/m²
qh_all    = diag["qh"]                     # W/m²

cal_mask    = Date(2004,1,1) .<= dates_all .<= Date(2013,12,31)
dates_mod   = dates_all[cal_mask]
nee_mod     = nee_all[cal_mask]
qle_mod     = qle_all[cal_mask]
qh_mod      = qh_all[cal_mask]

# ── Load observations ──────────────────────────────────────────────────────────
obs = JLD2.load(joinpath(@__DIR__, "observations.jld2"))
obs_dates = obs["obs_dates"]
y_obs     = obs["y_obs"]
n_obs     = length(obs_dates)
nee_obs_d = y_obs[1:n_obs]
qle_obs_d = y_obs[n_obs+1:2n_obs]
qh_obs_d  = y_obs[2n_obs+1:3n_obs]

# ── Align model to obs dates ───────────────────────────────────────────────────
# Build a dict for fast lookup of model day → index
date_to_idx = Dict(d => i for (i,d) in enumerate(dates_mod))

nee_mod_matched = [get(date_to_idx, d, nothing) !== nothing ? nee_mod[date_to_idx[d]] : NaN for d in obs_dates]
qle_mod_matched = [get(date_to_idx, d, nothing) !== nothing ? qle_mod[date_to_idx[d]] : NaN for d in obs_dates]
qh_mod_matched  = [get(date_to_idx, d, nothing) !== nothing ? qh_mod[date_to_idx[d]] : NaN for d in obs_dates]

# ── Per-year comparison table ──────────────────────────────────────────────────
println("\nPer-year annual means (obs days only, 2004-2013):")
println("Year | NEE_obs | NEE_mod | Qle_obs | Qle_mod | Qh_obs | Qh_mod")
for yr in 2004:2013
    ym = year.(obs_dates) .== yr
    no = mean(nee_obs_d[ym]); nm = mean(filter(!isnan, nee_mod_matched[ym]))
    lo = mean(qle_obs_d[ym]); lm = mean(filter(!isnan, qle_mod_matched[ym]))
    ho = mean(qh_obs_d[ym]);  hm = mean(filter(!isnan, qh_mod_matched[ym]))
    @printf(" %d | %+5.2f   | %+5.2f   | %5.1f   | %5.1f   | %5.1f  | %5.1f\n",
        yr, no, nm, lo, lm, ho, hm)
end

# Overall RMSE (obs-matched days)
valid = .!isnan.(nee_mod_matched)
nee_rmse = sqrt(mean((nee_obs_d[valid] .- nee_mod_matched[valid]).^2))
qle_rmse = sqrt(mean((qle_obs_d[valid] .- qle_mod_matched[valid]).^2))
qh_rmse  = sqrt(mean((qh_obs_d[valid]  .- qh_mod_matched[valid]).^2))
@printf("\nOverall RMSE: NEE=%.2f gC/m²/d, Qle=%.1f W/m², Qh=%.1f W/m²\n",
    nee_rmse, qle_rmse, qh_rmse)

# ── Figure: full 10-year timeseries ───────────────────────────────────────────
fig = CairoMakie.Figure(size=(1400, 900), fontsize=13)

ax1 = CairoMakie.Axis(fig[1,1],
    title = "DK-Sor 2004–2013 — NEE (gC m⁻² d⁻¹) | UKI iter 9, member 016",
    ylabel = "NEE (gC/m²/d)", xticklabelsvisible = false)
ax2 = CairoMakie.Axis(fig[2,1],
    title = "Latent Heat (W/m²)",
    ylabel = "Qle (W/m²)", xticklabelsvisible = false)
ax3 = CairoMakie.Axis(fig[3,1],
    title = "Sensible Heat (W/m²)",
    ylabel = "Qh (W/m²)")

# x-axis as numeric day index
xs_mod = 1:length(dates_mod)
# Build matching x for obs
obs_x = Float64[]
for d in obs_dates
    idx = get(date_to_idx, d, nothing)
    idx !== nothing && push!(obs_x, Float64(idx))
end

CairoMakie.lines!(ax1, xs_mod, nee_mod, color=:steelblue, linewidth=0.8, label="Model")
CairoMakie.scatter!(ax1, obs_x, nee_obs_d, color=(:forestgreen,0.6), markersize=3, label="Obs")
CairoMakie.hlines!(ax1, [0], color=:black, linewidth=0.5, linestyle=:dash)

CairoMakie.lines!(ax2, xs_mod, qle_mod, color=:steelblue, linewidth=0.8, label="Model")
CairoMakie.scatter!(ax2, obs_x, qle_obs_d, color=(:forestgreen,0.6), markersize=3, label="Obs")

CairoMakie.lines!(ax3, xs_mod, qh_mod, color=:steelblue, linewidth=0.8, label="Model")
CairoMakie.scatter!(ax3, obs_x, qh_obs_d, color=(:forestgreen,0.6), markersize=3, label="Obs")

# Add year separators and labels
year_starts = [findfirst(d -> d >= Date(yr,1,1), dates_mod) for yr in 2005:2013]
for (yr, xs) in zip(2005:2013, year_starts)
    xs === nothing && continue
    for ax in [ax1, ax2, ax3]
        CairoMakie.vlines!(ax, [Float64(xs)], color=(:gray, 0.4), linewidth=1)
    end
    CairoMakie.text!(ax1, Float64(xs)+5, 3.5, text="$yr", fontsize=10, color=:gray30)
end

CairoMakie.axislegend(ax1, position=:lt, framevisible=false)

# ── Figure 2: annual mean bar chart (obs vs model) ────────────────────────────
yrs = 2004:2013
nee_obs_yr = [mean(nee_obs_d[year.(obs_dates) .== yr]) for yr in yrs]
nee_mod_yr = [mean(filter(!isnan, nee_mod_matched[year.(obs_dates) .== yr])) for yr in yrs]
qle_obs_yr = [mean(qle_obs_d[year.(obs_dates) .== yr]) for yr in yrs]
qle_mod_yr = [mean(filter(!isnan, qle_mod_matched[year.(obs_dates) .== yr])) for yr in yrs]
qh_obs_yr  = [mean(qh_obs_d[year.(obs_dates) .== yr]) for yr in yrs]
qh_mod_yr  = [mean(filter(!isnan, qh_mod_matched[year.(obs_dates) .== yr])) for yr in yrs]

fig2 = CairoMakie.Figure(size=(1000, 750), fontsize=13)
yr_labels = string.(collect(yrs))

ax_n = CairoMakie.Axis(fig2[1,1], ylabel="Mean NEE (gC/m²/d)",
    xticks=(1:10, yr_labels), title="DK-Sor Annual Means: obs vs UKI posterior (iter 9)")
ax_l = CairoMakie.Axis(fig2[2,1], ylabel="Mean Qle (W/m²)",
    xticks=(1:10, yr_labels))
ax_h = CairoMakie.Axis(fig2[3,1], ylabel="Mean Qh (W/m²)",
    xticks=(1:10, yr_labels))

w = 0.35
xs = 1:10
CairoMakie.barplot!(ax_n, xs .- w/2, nee_obs_yr, width=w, color=:forestgreen, label="Obs")
CairoMakie.barplot!(ax_n, xs .+ w/2, nee_mod_yr, width=w, color=:steelblue, label="Model")
CairoMakie.hlines!(ax_n, [0], color=:black, linewidth=0.5)

CairoMakie.barplot!(ax_l, xs .- w/2, qle_obs_yr, width=w, color=:forestgreen, label="Obs")
CairoMakie.barplot!(ax_l, xs .+ w/2, qle_mod_yr, width=w, color=:steelblue, label="Model")

CairoMakie.barplot!(ax_h, xs .- w/2, qh_obs_yr, width=w, color=:forestgreen, label="Obs")
CairoMakie.barplot!(ax_h, xs .+ w/2, qh_mod_yr, width=w, color=:steelblue, label="Model")

CairoMakie.axislegend(ax_n, position=:rt, framevisible=false)

# ── Save ──────────────────────────────────────────────────────────────────────
out1 = joinpath(@__DIR__, "calibrate_dk_sor_output", "multiyear_timeseries.png")
out2 = joinpath(@__DIR__, "calibrate_dk_sor_output", "multiyear_annual_means.png")
CairoMakie.save(out1, fig)
CairoMakie.save(out2, fig2)
println("\nSaved → $out1")
println("Saved → $out2")
