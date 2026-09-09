#=
CONTINUOUS CalLMIP output — one task = one continuous 1997-2014 integration
(preserves inter-annual soil/carbon memory; no per-year cold-start). Driven by a
SLURM array via CALLMIP_TASK (= SLURM_ARRAY_TASK_ID):

  task 0          -> PRIOR     (default params)              -> output_callmip_sims/prior/callmip_diagnostics.jld2
  task 1          -> POSTERIOR  (calibrated mean, posterior_mean.jld2) -> output_callmip_sims/posterior/callmip_diagnostics.jld2
  task 2..(N+1)   -> ensemble member k=task-2 (posterior_samples.jld2) -> output_callmip_sims/members/member_<k>.jld2

The POSTERIOR run reads the deliverable posterior mean from posterior_mean.jld2
(the EKI/UTKI mean written by make_eki_posterior.jl) and the ensemble members from
posterior_samples.jld2; the draw uses a fixed RNG (MersenneTwister(2024)) so the
member set is reproducible. assemble_continuous_ensemble.jl then stacks the members
into ensemble_diagnostics.jld2 and write_callmip_netcdf.jl writes the 2 NetCDFs.
=#
import ClimaComms; ClimaComms.@import_required_backends
using ClimaLand, ClimaDiagnostics, ClimaUtilities, ClimaCore, NCDatasets, JLD2
using Dates, Statistics, Random
include(joinpath(@__DIR__, "model_runner.jl"))   # run_prior_year (supports stop_year)

const Y0     = 1997
const Y1     = 2014
const SIMDIR = joinpath(@__DIR__, "output_callmip_sims")
const CALDIR = joinpath(@__DIR__, "output_calibration")
const NMEM   = parse(Int, get(ENV, "CALLMIP_N_MEMBERS", "30"))
const TASK   = parse(Int, get(ENV, "CALLMIP_TASK", get(ENV, "SLURM_ARRAY_TASK_ID", "0")))

# Spin up the WATER + TEMPERATURE states before the output run: cycle the forcing
# for SPINUP_YEARS, take the equilibrated final state, then run 1997-2014 from it.
# Removes the ~2-3 yr soil-moisture (mrso) cold-start transient. Carbon (SOC/CO2/O2)
# is NOT spun up — it is prescribed from ClimaLand defaults. Same params for both.
# NOTE: 2014-12-31 is still fill (forcing ends 2014-12-31 → final :daily window can't close).
const SPINUP_YEARS = parse(Int, get(ENV, "CALLMIP_SPINUP_YEARS", "5"))
function runcont(ov)
    Yspin = run_prior_year(Y0; stop_year = Y0 + SPINUP_YEARS - 1, with_co2 = true,
                           dt = 900.0, param_overrides = ov, return_state = true)
    run_prior_year(Y0; stop_year = Y1, callmip = true, with_co2 = true,
                   dt = 900.0, param_overrides = ov, init_state = Yspin)
end

if TASK == 0
    @info "=== PRIOR continuous $Y0-$Y1 ==="
    surf, col, z, dates = runcont(Dict{String,Float64}())
    out = joinpath(SIMDIR, "prior"); mkpath(out)
    jldsave(joinpath(out, "callmip_diagnostics.jld2");
        dates = collect(dates), surface_data = surf, column_data = col, z_soil = z)
    @info "PRIOR done: $(length(dates)) days"

elseif TASK == 1
    @info "=== POSTERIOR MEAN continuous $Y0-$Y1 ==="
    pm = JLD2.load(joinpath(CALDIR, "posterior_mean.jld2"))
    ov = Dict(string.(pm["param_names"]) .=> Float64.(pm["posterior_mean"]))
    surf, col, z, dates = runcont(ov)
    out = joinpath(SIMDIR, "posterior"); mkpath(out)
    jldsave(joinpath(out, "callmip_diagnostics.jld2");
        dates = collect(dates), surface_data = surf, column_data = col, z_soil = z)
    @info "POSTERIOR done: $(length(dates)) days"

else
    k = TASK - 2                       # 0-based member index
    ps = JLD2.load(joinpath(CALDIR, "posterior_samples.jld2"))
    cp = ps["constrained_posterior"]; pn = string.(ps["param_names"])
    n_mcmc = size(cp, 2)
    draw = sort(randperm(MersenneTwister(2024), n_mcmc)[1:min(NMEM, n_mcmc)])
    @assert 0 <= k < length(draw) "member task $TASK out of range (N=$(length(draw)))"
    j = draw[k + 1]
    @info "=== ENSEMBLE member $(k+1)/$(length(draw)) continuous $Y0-$Y1 ==="
    ov = Dict(pn .=> Float64.(cp[:, j]))
    surf, _, _, dates = runcont(ov)
    out = joinpath(SIMDIR, "members"); mkpath(out)
    jldsave(joinpath(out, "member_$(k).jld2");
        nee = surf["nee"], lhf = surf["lhf"], shf = surf["shf"],
        dates = collect(dates), member = k + 1, param_names = pn)
    @info "member $(k+1) done: $(length(dates)) days"
end
@info "=== task $TASK complete ==="
