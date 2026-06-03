"""
ClimaCalibrate driver for DK-Sor single-site calibration.

Runs the PR1693 DK-Sor workflow with Alexis-style minibatching: an
ObservationSeries over fixed-size windows and a FixedMinibatcher with
batch size = 2.

Usage:
    julia --project=.buildkite experiments/calibrate_dk_sor/run_calibration.jl
"""

using Dates
using Distributed
import Random
import ClimaCalibrate
import ClimaLand
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2
using LinearAlgebra

# ── Configuration ────────────────────────────────────────────────────────────

const SITE_ID = "DK-Sor"
const N_ITERATIONS = 10
const MINIBATCH_SIZE = 2
const DT = Float64(900)

const climaland_dir = abspath(joinpath(@__DIR__, "..", ".."))
const OUTPUT_DIR = joinpath(climaland_dir, "experiments/calibrate_dk_sor/output")
const OBS_FILEPATH =
    joinpath(climaland_dir, "experiments/calibrate_dk_sor/observations.jld2")

# ── Priors ───────────────────────────────────────────────────────────────────
# Loaded from the shared priors.jl file; emulate_sample.jl uses the same function
# so both scripts are guaranteed to use identical prior definitions.
include(joinpath(@__DIR__, "priors.jl"))
prior, priors = build_dk_sor_priors()

# ── Load Observations ────────────────────────────────────────────────────────

obs_data = JLD2.load(OBS_FILEPATH)
observation_vector = obs_data["observation_vector"]
window_names = obs_data["window_names"]
window_dates = obs_data["window_dates"]

obs_series = EKP.ObservationSeries(
    Dict(
        "observations" => observation_vector,
        "names" => window_names,
        "minibatcher" => ClimaCalibrate.minibatcher_over_samples(
            length(observation_vector),
            MINIBATCH_SIZE,
        ),
    ),
)

println("Loaded $(length(observation_vector)) observation windows ($(first(window_names)) … $(last(window_names)))")
println("Minibatch size: $MINIBATCH_SIZE")

# ── Create EKP ───────────────────────────────────────────────────────────────

rng = Random.MersenneTwister(1234)
ekp = EKP.EnsembleKalmanProcess(
    obs_series,
    EKP.TransformUnscented(prior; impose_prior = true),
    verbose = true,
    rng = rng,
    scheduler = EKP.DataMisfitController(terminate_at = 100),
)

N_ens = EKP.get_N_ens(ekp)
println("Ensemble size: $N_ens (for $(length(priors)) parameters)")

# ── Add SLURM Workers ───────────────────────────────────────────────────────

curr_backend = ClimaCalibrate.get_backend()
if curr_backend == ClimaCalibrate.CaltechHPCBackend ||
   curr_backend == ClimaCalibrate.GCPBackend
    @info "Check your slurm script that the number of tasks is the same as the ensemble size ($N_ens)"
    addprocs(ClimaCalibrate.SlurmManager())
elseif nworkers() == 1
    @info "No SLURM detected and only 1 worker — running serially (JuliaBackend fallback)"
end

# ── Broadcast Configuration to Workers ───────────────────────────────────────

@everywhere using Distributed
@everywhere import ClimaLand
@everywhere const SITE_ID = $SITE_ID
@everywhere const OUTPUT_DIR = $OUTPUT_DIR
@everywhere const OBS_FILEPATH = $OBS_FILEPATH
@everywhere const DT = $DT
@everywhere include(
    joinpath(abspath(joinpath(@__DIR__, "..", "..")), "experiments/calibrate_dk_sor/model_interface.jl"),
)

# ── Run Calibration ──────────────────────────────────────────────────────────

println("\n=== Starting ClimaCalibrate calibration ===")
println("  Backend: WorkerBackend")
println("  Iterations: $N_ITERATIONS")
println("  Ensemble size: $N_ens")
println("  Windows: $(length(window_names)) ($(first(window_names))–$(last(window_names))), minibatch=$MINIBATCH_SIZE")
println("  Days per window: $(length(window_dates[1]))")
println("  Parameters: $(length(priors))")
println("  LAI: Copernicus")
println("  Output: $OUTPUT_DIR")

eki = ClimaCalibrate.calibrate(
    ClimaCalibrate.WorkerBackend,
    ekp,
    N_ITERATIONS,
    prior,
    OUTPUT_DIR,
)

# ── Print Final Results ──────────────────────────────────────────────────────

final_params = EKP.get_ϕ_mean_final(prior, eki)
param_names = [only(PD.get_name(d)) for d in priors]
println("\n=== Calibration Complete ===")
println("Final parameter means:")
for (name, val) in zip(param_names, final_params)
    println("  $name = $(round(val, sigdigits = 5))")
end
