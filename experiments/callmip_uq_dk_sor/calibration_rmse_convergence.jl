"""
Quick calibration RMSE convergence analysis.
Reads daily_diagnostics.jld2 from calibrate_dk_sor/output and FLUXNET NetCDF.
Prints per-iteration ensemble-mean RMSE and monthly NEE comparison table.
"""

import JLD2, Statistics, Dates, NCDatasets
using Printf

cal_dir = joinpath(@__DIR__, "..", "calibrate_dk_sor")
obs_nc  = joinpath(@__DIR__, "..", "..", "DK_Sor",
    "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")

println("Loading FLUXNET observations...")
ds = NCDatasets.NCDataset(obs_nc)
nee_obs_v = Float64.(NCDatasets.coalesce.(ds["NEE_daily"][:], NaN))
qle_obs_v = Float64.(NCDatasets.coalesce.(ds["Qle_daily"][:], NaN))
qh_obs_v  = Float64.(NCDatasets.coalesce.(ds["Qh_daily"][:], NaN))
obs_dates  = Dates.Date.(ds["time"][:])
NCDatasets.close(ds)

obs_nee = Dict(zip(obs_dates, nee_obs_v))
obs_qle = Dict(zip(obs_dates, qle_obs_v))
obs_qh  = Dict(zip(obs_dates, qh_obs_v))

function load_diag(iter, member)
    p = joinpath(cal_dir, "output",
        "iteration_$(lpad(iter,3,'0'))",
        "member_$(lpad(member,3,'0'))",
        "daily_diagnostics.jld2")
    isfile(p) || return nothing
    d = JLD2.load(p)
    return (
        dates = Dates.Date.(d["dates"]),
        nee   = d["nee"] .* 12.0 .* 86400.0,   # mol/m²/s → gC/m²/d
        qle   = d["qle"],
        qh    = d["qh"],
    )
end

function rmse_flux(dates, vals, obs)
    e2 = [(vals[i] - obs[dates[i]])^2
           for i in eachindex(dates)
           if haskey(obs, dates[i]) && !isnan(obs[dates[i]])]
    isempty(e2) ? NaN : sqrt(Statistics.mean(e2))
end

println("\nIter  |  NEE RMSE  |  Qle RMSE  |  Qh RMSE  | n")
println("-"^55)
for iter in 0:10
    nr = Float64[];  qr = Float64[];  hr = Float64[]
    for m in 1:33
        r = load_diag(iter, m);  isnothing(r) && continue
        v = rmse_flux(r.dates, r.nee, obs_nee);  isnan(v) || push!(nr, v)
        v = rmse_flux(r.dates, r.qle, obs_qle);  isnan(v) || push!(qr, v)
        v = rmse_flux(r.dates, r.qh,  obs_qh);   isnan(v) || push!(hr, v)
    end
    isempty(nr) && continue
    @printf("  %2d  |  %7.4f   |  %7.3f   |  %7.3f  | %d\n",
        iter,
        Statistics.mean(nr),
        Statistics.mean(qr),
        Statistics.mean(hr),
        length(nr))
end

# ── Prior vs best posterior head-to-head ────────────────────────────────────
println()
r0 = load_diag(0, 1)

function find_best_member(obs_nee)
    best_r = Inf;  bm = 1;  best_iter = 0
    for iter in 0:10, m in 1:33
        r = load_diag(iter, m);  isnothing(r) && continue
        v = rmse_flux(r.dates, r.nee, obs_nee)
        if v < best_r
            best_r = v; bm = m; best_iter = iter
        end
    end
    return best_r, bm, best_iter
end
best_r, bm, best_iter = find_best_member(obs_nee)
r10 = load_diag(best_iter, bm)

@printf("Prior  (iter 0,  m001): NEE=%.4f  Qle=%.3f  Qh=%.3f\n",
    rmse_flux(r0.dates,  r0.nee,  obs_nee),
    rmse_flux(r0.dates,  r0.qle,  obs_qle),
    rmse_flux(r0.dates,  r0.qh,   obs_qh))
@printf("Post   (iter %2d, m%03d): NEE=%.4f  Qle=%.3f  Qh=%.3f\n",
    best_iter, bm,
    best_r,
    rmse_flux(r10.dates, r10.qle, obs_qle),
    rmse_flux(r10.dates, r10.qh,  obs_qh))

# ── Monthly NEE comparison ────────────────────────────────────────────────────
println()
println("Monthly NEE: Prior / Post / FLUXNET (gC/m²/d)")
println("  Mon    Prior    Post     Obs")
println("  ─────────────────────────────")
for mo in 1:12
    om = Statistics.mean(filter(!isnan,
        [get(obs_nee, d, NaN) for d in obs_dates if Dates.month(d) == mo]))
    pm = Statistics.mean(
        [v for (d, v) in zip(r0.dates,  r0.nee)  if Dates.month(d) == mo])
    qm = Statistics.mean(
        [v for (d, v) in zip(r10.dates, r10.nee) if Dates.month(d) == mo])
    @printf("  %-3s  %+7.3f  %+7.3f  %+7.3f\n",
        Dates.monthabbr(mo), pm, qm, om)
end
println()
println("Done.")
