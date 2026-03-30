"""
Single forward run with Alexis's 12-param EKI posterior parameters
(extracted from iteration_000_12param_backup/eki_file.jld2).

Parameters are stored in:
    experiments/calibrate_dk_sor/output/iteration_999/member_001/parameters.toml

Runs the DK-Sor model for 2003-01-01 to 2014-01-01, saves daily_diagnostics.jld2,
then computes RMSE vs observations and plots seasonal cycle.

Usage:
    julia --project=.buildkite experiments/calibrate_dk_sor/run_alexis_test.jl
"""

using Dates, Statistics
import ClimaLand
import ClimaCalibrate
import JLD2
using CairoMakie

const SITE_ID      = "DK-Sor"
const OUTPUT_DIR   = abspath(joinpath(@__DIR__, "..", "..", "experiments/calibrate_dk_sor/output"))
const OBS_FILEPATH = abspath(joinpath(@__DIR__, "..", "..", "experiments/calibrate_dk_sor/observations.jld2"))
const OUT_DIR      = abspath(joinpath(@__DIR__, "..", "..", "experiments/calibrate_dk_sor"))
const DT           = Float64(450)

# Single-worker serial run — no Distributed needed.
# model_interface.jl defines ClimaCalibrate.forward_model(iteration, member).
include(
    abspath(joinpath(@__DIR__, "..", "..", "experiments/calibrate_dk_sor/model_interface.jl"))
)

println("Running forward model with Alexis posterior params (iteration=999, member=1)…")
ClimaCalibrate.forward_model(999, 1)
println("Forward run complete.")

# ── Load results and compute RMSE ─────────────────────────────────────────────

diag_path = ClimaCalibrate.path_to_ensemble_member(OUTPUT_DIR, 999, 1)
d    = JLD2.load(joinpath(diag_path, "daily_diagnostics.jld2"))
dts  = d["dates"]::Vector{Date}
nee  = d["nee"]::Vector{Float64}
qle  = d["qle"]::Vector{Float64}
qh   = d["qh"]::Vector{Float64}

obs  = JLD2.load(OBS_FILEPATH)
y    = obs["y_obs"]::Vector{Float64}
obs_dates = obs["obs_dates"]::Vector{Date}
n_obs     = length(obs_dates)
y_nee = y[1:n_obs]; y_qle = y[n_obs+1:2n_obs]; y_qh = y[2n_obs+1:end]

# Match model to obs dates
mod_dict_nee = Dict(zip(dts, nee))
mod_dict_qle = Dict(zip(dts, qle))
mod_dict_qh  = Dict(zip(dts, qh))

nee_mod = [get(mod_dict_nee, d, NaN) for d in obs_dates] .* 12.0 .* 86400.0  # mol CO2/m²/s → gC/m²/d
qle_mod = [get(mod_dict_qle, d, NaN) for d in obs_dates]
qh_mod  = [get(mod_dict_qh,  d, NaN) for d in obs_dates]

rmse(a, b) = let ok=.!isnan.(a).&&.!isnan.(b); sqrt(mean((a[ok].-b[ok]).^2)) end
println("\nRMSE vs observations:")
println("  NEE: $(round(rmse(nee_mod, y_nee), digits=3)) gC/m²/d")
println("  Qle: $(round(rmse(qle_mod, y_qle), digits=2)) W/m²")
println("  Qh:  $(round(rmse(qh_mod,  y_qh),  digits=2)) W/m²")

# ── Plot seasonal cycle ────────────────────────────────────────────────────────

function monthly_mean_vals(dates, vals; cal_start=Date(2004,1,1), cal_end=Date(2013,12,31))
    by_month = [Float64[] for _ in 1:12]
    for (dt, v) in zip(dates, vals)
        (isnan(v) || dt < cal_start || dt > cal_end) && continue
        push!(by_month[Dates.month(dt)], v)
    end
    mn = [isempty(v) ? NaN : mean(v) for v in by_month]
    sd = [length(v) < 2 ? NaN : std(v) for v in by_month]
    return mn, sd
end

obs_cal_dates = obs_dates  # already filtered to 2004-2013
nee_obs_mn, nee_obs_sd = monthly_mean_vals(obs_cal_dates, y_nee)
qle_obs_mn, qle_obs_sd = monthly_mean_vals(obs_cal_dates, y_qle)
qh_obs_mn,  qh_obs_sd  = monthly_mean_vals(obs_cal_dates, y_qh)

nee_mod_mn, _ = monthly_mean_vals(obs_cal_dates, nee_mod)
qle_mod_mn, _ = monthly_mean_vals(obs_cal_dates, qle_mod)
qh_mod_mn,  _ = monthly_mean_vals(obs_cal_dates, qh_mod)

fig = Figure(size=(900, 950), fontsize=13)
MO_LAB = ["J","F","M","A","M","J","J","A","S","O","N","D"]
x = collect(1.0:12.0)

vars = [
    (nee_obs_mn, nee_obs_sd, nee_mod_mn, "NEE (gC m⁻² d⁻¹)", rmse(nee_mod, y_nee)),
    (qle_obs_mn, qle_obs_sd, qle_mod_mn, "Qle (W m⁻²)",       rmse(qle_mod, y_qle)),
    (qh_obs_mn,  qh_obs_sd,  qh_mod_mn,  "Qh (W m⁻²)",        rmse(qh_mod,  y_qh)),
]

for (row, (obs_mn, obs_sd, mod_mn, ylabel, r)) in enumerate(vars)
    ax = Axis(fig[row, 1];
              title   = "$ylabel  |  RMSE = $(round(r, digits=2))",
              ylabel  = ylabel,
              xlabel  = row == 3 ? "Month" : "",
              xticks  = (1:12, MO_LAB))
    band!(ax, x, obs_mn .- obs_sd, obs_mn .+ obs_sd; color=(:black, 0.15))
    lines!(ax, x, obs_mn;  color=:black,     linewidth=2.5, label="Observed")
    scatter!(ax, x, obs_mn; color=:black,    markersize=7)
    lines!(ax, x, mod_mn;  color=:firebrick, linewidth=2.5, label="Alexis posterior")
    scatter!(ax, x, mod_mn; color=:firebrick, markersize=7)
    axislegend(ax; position=:lt, framevisible=false, labelsize=11)
end

Label(fig[0, 1]; text="DK-Sor: Alexis 12-param Posterior  |  2004–2013",
      fontsize=15, font=:bold)
rowgap!(fig.layout, 8)

out_path = joinpath(OUT_DIR, "alexis_test_seasonality.png")
save(out_path, fig; px_per_unit=2)
println("\nPlot saved → $out_path")
