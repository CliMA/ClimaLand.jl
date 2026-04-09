"""
ClimaCalibrate driver for DK-Sor single-site calibration.

Calibrates 16 parameters (9 canopy + 3 DAMM soilCO2 + 3 respiration + ac_canopy) against daily NEE, Qle, Qh
using TransformUnscented Kalman Inversion (N_ens = 33). All ~10 years of observations
(2004-2013, wind-filtered) are used at each iteration (no minibatching).

Copernicus LAI is used for the vegetation forcing. DT = 900 s.

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
y_obs = obs_data["y_obs"]
noise_cov = obs_data["noise_cov"]
obs_dates = obs_data["obs_dates"]
cal_years = obs_data["cal_years"]

n_obs = length(obs_dates)
println("Loaded $n_obs valid observation days ($(first(cal_years))-$(last(cal_years)))")
println("Observation vector length: $(length(y_obs)) (3 × $n_obs)")

# ── Create EKP ───────────────────────────────────────────────────────────────

rng = Random.MersenneTwister(1234)
ekp = EKP.EnsembleKalmanProcess(
    y_obs,
    noise_cov,
    EKP.TransformUnscented(prior; impose_prior = true),
    verbose = true,
    rng = rng,
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
println("  Observation days: $n_obs (wind-filtered, $(first(cal_years))-$(last(cal_years)))")
println("  Parameters: $(length(priors)) (9 canopy + 3 DAMM soilCO2)")
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
