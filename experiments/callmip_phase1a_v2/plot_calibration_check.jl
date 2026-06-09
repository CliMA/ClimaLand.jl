"""
Plot calibration diagnostics:
  1. EKI convergence (observation misfit vs iteration)
  2. Seasonal climatology model vs obs for NEE, Qle, Qh

Usage:
    julia --project=experiments/callmip_phase1a_v2 \\
          experiments/callmip_phase1a_v2/plot_calibration_check.jl
"""

using Dates
using Statistics
import JLD2
import ClimaCalibrate
using NCDatasets
using Plots

const climaland_dir   = abspath(joinpath(@__DIR__, "..", ".."))
const exp_dir         = @__DIR__
const flux_nc_path    = joinpath(climaland_dir, "DK_Sor",
                                  "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")
const cal_output_dir  = joinpath(exp_dir, "output_eki")
const callmip_sim_dir = joinpath(exp_dir, "output_callmip_sims")
const fig_dir         = joinpath(exp_dir, "figures")
isdir(fig_dir) || mkpath(fig_dir)

const mol_CO2_to_gC_day = 12.0 * 86400.0
const CAL_STOP = Date(2012, 12, 31)

# ── 1. EKI convergence ────────────────────────────────────────────────────────
iter = 0
rmse_vec = Float64[]
while true
    ipath = joinpath(cal_output_dir, "iteration_$(lpad(iter, 3, '0'))", "eki_file.jld2")
    isfile(ipath) || break
    d = JLD2.load(ipath)
    if haskey(d, "rmse")
        push!(rmse_vec, d["rmse"])
    end
    global iter += 1
end

if !isempty(rmse_vec)
    fig_conv = plot(0:length(rmse_vec)-1, rmse_vec;
                    xlabel = "EKI iteration",
                    ylabel = "Observation misfit (normalised RMSE)",
                    title  = "EKI convergence — DK-Sor CalLMIP",
                    marker = :circle, lw = 2, legend = false)
    savefig(fig_conv, joinpath(fig_dir, "01_eki_convergence.png"))
    @info "Saved EKI convergence plot"
else
    @warn "No RMSE data found in EKI output — skipping convergence plot"
end

# ── 2. Seasonal climatology ────────────────────────────────────────────────────
ds = NCDatasets.NCDataset(flux_nc_path, "r")
flux_dates  = Date.(ds["time"][:])
nee_obs_all = Float64.(coalesce.(ds["NEE_daily"][:], NaN))
qle_obs_all = Float64.(coalesce.(ds["Qle_daily"][:], NaN))
qh_obs_all  = Float64.(coalesce.(ds["Qh_daily"][:],  NaN))
close(ds)
d2i_flux = Dict(dt => i for (i, dt) in enumerate(flux_dates))

function load_sim_diag(m)
    path = joinpath(
        ClimaCalibrate.path_to_ensemble_member(callmip_sim_dir, 0, m),
        "callmip_diagnostics.jld2",
    )
    isfile(path) || error("Missing $path — run run_callmip_simulations.jl first")
    d = JLD2.load(path)
    d["dates"], d["surface_data"]
end

prior_dates, prior_sd = load_sim_diag(1)
post_dates,  post_sd  = load_sim_diag(2)

function monthly_mean(vals, dates, obs_raw, obs_d2i; cal_only=false)
    model_by_mon = Dict(m => Float64[] for m in 1:12)
    obs_by_mon   = Dict(m => Float64[] for m in 1:12)
    for (k, dt) in enumerate(dates)
        cal_only && dt > CAL_STOP && continue
        k > length(vals) && continue
        mv = vals[k]; isnan(mv) && continue
        push!(model_by_mon[month(dt)], mv)
        if haskey(obs_d2i, dt)
            ov = obs_raw[obs_d2i[dt]]
            !isnan(ov) && push!(obs_by_mon[month(dt)], ov)
        end
    end
    return ([mean(get(model_by_mon, m, [NaN])) for m in 1:12],
            [mean(get(obs_by_mon,   m, [NaN])) for m in 1:12])
end

months = 1:12
mnames = ["J","F","M","A","M","J","J","A","S","O","N","D"]

for (varname, prior_vals, post_vals, obs_raw, ylabel, scale) in [
    ("NEE", get(prior_sd, "nee", Float64[]) .* mol_CO2_to_gC_day,
            get(post_sd,  "nee", Float64[]) .* mol_CO2_to_gC_day,
            nee_obs_all, "NEE (gC/m²/d)", 1.0),
    ("Qle", get(prior_sd, "lhf", Float64[]),
            get(post_sd,  "lhf", Float64[]),
            qle_obs_all, "Qle (W/m²)", 1.0),
    ("Qh",  get(prior_sd, "shf", Float64[]),
            get(post_sd,  "shf", Float64[]),
            qh_obs_all,  "Qh (W/m²)", 1.0),
]
    pr_mon, obs_mon = monthly_mean(prior_vals, prior_dates, obs_raw, d2i_flux; cal_only=true)
    po_mon, _       = monthly_mean(post_vals,  post_dates,  obs_raw, d2i_flux; cal_only=true)

    fig = plot(months, obs_mon;
               label = "FLUXNET obs", lw = 2, lc = :black, marker = :circle,
               xlabel = "Month", ylabel = ylabel,
               title  = "DK-Sor $varname seasonal climatology (1997-2012)")
    plot!(fig, months, pr_mon; label = "Prior",     lw = 2, lc = :dodgerblue,  ls = :dash)
    plot!(fig, months, po_mon; label = "Posterior", lw = 2, lc = :orangered,   ls = :solid)
    xticks!(fig, months, mnames)

    fpath = joinpath(fig_dir, "02_seasonal_$(varname).png")
    savefig(fig, fpath)
    @info "Saved $fpath"
end

println("Plots saved to $fig_dir")
