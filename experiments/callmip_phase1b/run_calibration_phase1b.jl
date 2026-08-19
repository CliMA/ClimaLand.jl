"""
UTKI calibration for CalLMIP Phase 1b, Scenario 1 (multisite).

Faithful port of experiments/callmip_dksor/run_calibration.jl to 21 sites:
UTKI (TransformUnscented, impose_prior, 2p+1 = 23 members for the 11 Phase 1a
carbon parameters), biennial consecutive minibatches over the 2003–2014
observation windows (G dim = 2 × 756 = 1512 per iteration), member-parallel
pmap on CPU workers (Phase 3 decision: one 21-column ColumnEnsemble per member
per core; a member run = 1 spin-up year + 2 window years ≈ 15 min).

Test gate (Phase 4): serial, 2 pseudo-members, 1 iteration, single group:
  julia --project=.buildkite experiments/callmip_phase1b/run_calibration_phase1b.jl --test
Full run:
  julia -p 23 --project=.buildkite experiments/callmip_phase1b/run_calibration_phase1b.jl
Checkpointed per iteration (ekp_checkpoint.jld2); rerun to resume.
"""

using Distributed

const IS_TEST = "--test" in ARGS

if !IS_TEST && nworkers() == 1
    const N_WORKERS = parse(Int, get(ENV, "CALLMIP1B_WORKERS", "23"))
    addprocs(N_WORKERS; exeflags = "--project=$(Base.active_project())")
end
@info IS_TEST ? "TEST mode (serial)" : "Using $(nworkers()) workers"

@everywhere import ClimaComms
@everywhere ClimaComms.@import_required_backends
@everywhere include(joinpath(@__DIR__, "forward_model_phase1b.jl"))

import EnsembleKalmanProcesses as EKP
using JLD2
using Random
using Statistics
using LinearAlgebra
using Dates

const OUTDIR =
    get(ENV, "CALLMIP1B_OUTDIR", joinpath(@__DIR__, "output_calibration"))
const CHECKPOINT_PATH = joinpath(OUTDIR, "ekp_checkpoint.jld2")
const N_ITERATIONS =
    parse(Int, get(ENV, "CALLMIP1B_N_ITER", IS_TEST ? "1" : "10"))

include(joinpath(@__DIR__, "..", "callmip_dksor", "priors.jl")) # build_dk_sor_priors: the 11-param set (D2)

# ── Observations ─────────────────────────────────────────────────────────────
obs_data = jldopen(joinpath(OUTDIR, "observations.jld2"), "r") do f
    (;
        calib_years = f["calib_years"],
        site_ids = f["site_ids"],
        obs_series = f["obs_series"],
    )
end
calib_years = obs_data.calib_years
@assert obs_data.site_ids == CALIBRATION_SITE_IDS "site order mismatch obs vs forward model"

# Biennial consecutive groups (Phase 1a convention; drop trailing singleton).
group_idx = [
    collect(i:min(i + 1, length(calib_years))) for i in 1:2:length(calib_years)
]
if length(group_idx[end]) == 1
    @warn "Dropping trailing singleton year $(calib_years[group_idx[end][1]])"
    pop!(group_idx)
end
group_years = [calib_years[g] for g in group_idx]
group_obs =
    [EKP.combine_observations(obs_data.obs_series[g]) for g in group_idx]
IS_TEST && (group_years = group_years[1:1]; group_obs = group_obs[1:1])
obs_series_ekp =
    EKP.ObservationSeries(group_obs, EKP.RandomFixedSizeMinibatcher(1))
@info "Minibatch groups: $(group_years)"

# ── Prior + EKP (UTKI) ───────────────────────────────────────────────────────
prior, _ = build_dk_sor_priors()
rng = Random.MersenneTwister(42)
build_ekp() = EKP.EnsembleKalmanProcess(
    obs_series_ekp,
    EKP.TransformUnscented(prior; impose_prior = true);
    scheduler = EKP.DataMisfitController(terminate_at = 100),
    verbose = true,
    rng,
)

function load_checkpoint()
    isfile(CHECKPOINT_PATH) || return nothing, 0
    d = JLD2.load(CHECKPOINT_PATH)
    @info "Resuming from iteration $(d["iteration"])"
    return d["ekp"], d["iteration"]
end

function run_ensemble(ekp, years)
    ϕ_ens = EKP.get_ϕ_final(prior, ekp)
    cols = 1:size(ϕ_ens, 2)
    # test mode: 2 members serially, no EKP update afterwards (a partial-NaN
    # ensemble would exercise the failure handler, not the pipeline)
    IS_TEST && return hcat(map(m -> run_member(ϕ_ens[:, m], years), 1:2)...)
    return hcat(pmap(m -> run_member(ϕ_ens[:, m], years), cols)...)
end

# ── Loop ─────────────────────────────────────────────────────────────────────
mkpath(OUTDIR)
ekp, start_iter = load_checkpoint()
if isnothing(ekp)
    ekp = build_ekp()
    start_iter = 0
    @info "UTKI ensemble size = $(EKP.get_N_ens(ekp))"
end

for iter in (start_iter + 1):N_ITERATIONS
    mb = EKP.get_current_minibatch(ekp)
    years = reduce(vcat, group_years[mb])
    @info "Iteration $iter/$N_ITERATIONS — window $years"
    t = @elapsed G_ens = run_ensemble(ekp, years)
    G_ens[.!isfinite.(G_ens)] .= NaN
    nfail = count(j -> any(isnan, G_ens[:, j]), 1:size(G_ens, 2))
    @info "  G_ens $(size(G_ens)); NaN members: $nfail; wall $(round(t / 60; digits = 1)) min"
    if IS_TEST
        finite_frac = count(isfinite, G_ens) / length(G_ens)
        obs = EKP.get_obs(group_obs[1])
        um = diag(EKP.get_obs_noise_cov(group_obs[1])) .< 1.0e11 # unmasked
        d = G_ens[:, 1] .- obs
        @info "TEST gate" finite_frac rms_misfit_unmasked =
            sqrt(mean(abs2, d[um .& isfinite.(d)]))
        @info "TEST mode: exiting before EKP update"
        exit(0)
    end
    EKP.update_ensemble!(ekp, G_ens)
    jldsave(CHECKPOINT_PATH; ekp, iteration = iter)
    ϕm = EKP.get_ϕ_mean_final(prior, ekp)
    for (n, v) in zip(PARAM_NAMES, ϕm)
        @info "    $(rpad(n, 34)) = $(round(v; sigdigits = 4))"
    end
end

ϕ_final = EKP.get_ϕ_final(prior, ekp)
jldsave(
    joinpath(OUTDIR, "posterior_ensemble.jld2");
    param_names = PARAM_NAMES,
    ϕ_final,
    ϕ_mean = EKP.get_ϕ_mean_final(prior, ekp),
    calib_years,
    group_years,
)
@info "Saved posterior ensemble → $(joinpath(OUTDIR, "posterior_ensemble.jld2"))"
