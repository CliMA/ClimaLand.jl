"""
Evaluate calibration: compute RMSE/bias table and seasonal plots for NEE, Qle, Qh.

Usage:
    julia --project=experiments/callmip_phase1a_v2 \\
          experiments/callmip_phase1a_v2/evaluate_calibration.jl
"""

using Dates
using Statistics
import JLD2
import ClimaCalibrate
using NCDatasets
using Printf

const climaland_dir  = abspath(joinpath(@__DIR__, "..", ".."))
const exp_dir        = @__DIR__
const flux_nc_path   = joinpath(climaland_dir, "DK_Sor",
                                 "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")
const callmip_sim_dir = joinpath(exp_dir, "output_callmip_sims")

const mol_CO2_to_gC_day = 12.0 * 86400.0
const START_DATE = Date(1997, 1, 1)
const CAL_STOP   = Date(2012, 12, 31)

# ── Load observations ─────────────────────────────────────────────────────────
ds = NCDatasets.NCDataset(flux_nc_path, "r")
flux_dates  = Date.(ds["time"][:])
nee_obs_all = Float64.(coalesce.(ds["NEE_daily"][:], NaN))
qle_obs_all = Float64.(coalesce.(ds["Qle_daily"][:], NaN))
qh_obs_all  = Float64.(coalesce.(ds["Qh_daily"][:],  NaN))
close(ds)

d2i_flux = Dict(dt => i for (i, dt) in enumerate(flux_dates))
obs_lookup = (
    nee = d2i_flux, qle = d2i_flux, qh = d2i_flux,
    nee_raw = nee_obs_all, qle_raw = qle_obs_all, qh_raw = qh_obs_all,
)

# ── Load model diagnostics (prior = member 1, posterior = member 2) ───────────
function load_member(m)
    path = joinpath(
        ClimaCalibrate.path_to_ensemble_member(callmip_sim_dir, 0, m),
        "callmip_diagnostics.jld2",
    )
    isfile(path) || error("Diagnostics not found: $path")
    d = JLD2.load(path)
    dates = d["dates"]
    sd    = d["surface_data"]
    nee_model = get(sd, "nee", Float64[]) .* mol_CO2_to_gC_day
    qle_model = get(sd, "lhf", Float64[])
    qh_model  = get(sd, "shf", Float64[])
    return dates, nee_model, qle_model, qh_model
end

prior_dates, prior_nee, prior_qle, prior_qh  = load_member(1)
post_dates,  post_nee,  post_qle,  post_qh   = load_member(2)

# ── Helper: RMSE/bias on matched days ─────────────────────────────────────────
function metrics(model_vals, model_dates, obs_raw, obs_d2i; cal_only = false)
    pairs_m = Float64[]
    pairs_o = Float64[]
    for (k, dt) in enumerate(model_dates)
        cal_only && dt > CAL_STOP && continue
        haskey(obs_d2i, dt) || continue
        ov = obs_raw[obs_d2i[dt]]
        isnan(ov) && continue
        k <= length(model_vals) || continue
        mv = model_vals[k]
        isnan(mv) && continue
        push!(pairs_m, mv); push!(pairs_o, ov)
    end
    isempty(pairs_m) && return (rmse = NaN, bias = NaN, n = 0)
    rmse = sqrt(mean((pairs_m .- pairs_o) .^ 2))
    bias = mean(pairs_m .- pairs_o)
    return (rmse = rmse, bias = bias, n = length(pairs_m))
end

println("\n╔═══════════════════════════════════════════════════════════════════════╗")
println("║          CalLMIP Phase 1a DK-Sor — Calibration Evaluation            ║")
println("╠══════════╦═══════════════════════════════════════╦═══════════════════╣")
println("║ Variable ║           Calibration (1997-2012)     ║   Full 1997-2013  ║")
println("║          ║  Prior RMSE / bias  |  Post RMSE/bias ║  Post RMSE/bias   ║")
println("╠══════════╬═══════════════════════════════════════╬═══════════════════╣")

for (varname, prior_vals, prior_dt, post_vals, post_dt, obs_raw, units) in [
    ("NEE",  prior_nee, prior_dates, post_nee, post_dates, nee_obs_all, "gC/m²/d"),
    ("Qle",  prior_qle, prior_dates, post_qle, post_dates, qle_obs_all, "W/m²"),
    ("Qh",   prior_qh,  prior_dates, post_qh,  post_dates, qh_obs_all,  "W/m²"),
]
    pr_cal  = metrics(prior_vals, prior_dt, obs_raw, d2i_flux; cal_only=true)
    po_cal  = metrics(post_vals,  post_dt,  obs_raw, d2i_flux; cal_only=true)
    po_full = metrics(post_vals,  post_dt,  obs_raw, d2i_flux; cal_only=false)

    @printf "║ %-8s ║ %5.2f / %+5.2f %5s  | %5.2f / %+5.2f %5s ║ %5.2f / %+5.2f %5s║\n" \
        varname pr_cal.rmse pr_cal.bias units \
        po_cal.rmse po_cal.bias units \
        po_full.rmse po_full.bias units
end
println("╚══════════╩═══════════════════════════════════════╩═══════════════════╝")

# ── Monthly climatology ───────────────────────────────────────────────────────
function monthly_climatology(vals, dates, obs_raw, obs_d2i; cal_only=false)
    model_by_mon = Dict(m => Float64[] for m in 1:12)
    obs_by_mon   = Dict(m => Float64[] for m in 1:12)
    for (k, dt) in enumerate(dates)
        cal_only && dt > CAL_STOP && continue
        k > length(vals) && continue
        mv = vals[k]
        isnan(mv) && continue
        push!(model_by_mon[month(dt)], mv)
        if haskey(obs_d2i, dt)
            ov = obs_raw[obs_d2i[dt]]
            !isnan(ov) && push!(obs_by_mon[month(dt)], ov)
        end
    end
    ([mean(model_by_mon[m]) for m in 1:12],
     [mean(obs_by_mon[m])   for m in 1:12])
end

println("\nMonthly NEE climatology (prior / posterior / obs, gC/m²/d):")
prior_nee_mon, obs_nee_mon = monthly_climatology(prior_nee, prior_dates, nee_obs_all, d2i_flux)
post_nee_mon,  _           = monthly_climatology(post_nee,  post_dates,  nee_obs_all, d2i_flux)
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
@printf "%-5s  %6s  %6s  %6s\n" "Month" "Prior" "Post" "Obs"
for m in 1:12
    @printf "%-5s  %+6.2f  %+6.2f  %+6.2f\n" months[m] prior_nee_mon[m] post_nee_mon[m] obs_nee_mon[m]
end

# Winter NEE deficit
djf_prior = mean(prior_nee_mon[[12, 1, 2]])
djf_post  = mean(post_nee_mon[[12, 1, 2]])
djf_obs   = mean(obs_nee_mon[[12, 1, 2]])
println("\nDJF mean NEE — Prior: $(round(djf_prior; digits=3)), Post: $(round(djf_post; digits=3)), Obs: $(round(djf_obs; digits=3)) gC/m²/d")
if abs(djf_post - djf_obs) > 0.5
    @warn "Winter (DJF) NEE deficit: model $(round(djf_post; digits=3)) vs obs $(round(djf_obs; digits=3)) gC/m²/d — $(round(abs(djf_post-djf_obs)/abs(djf_obs)*100; digits=1))% error"
end
