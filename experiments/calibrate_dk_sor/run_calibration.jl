"""
ClimaCalibrate driver for DK-Sor single-site calibration.

Calibrates 12 parameters (9 canopy + 3 DAMM soilCO2) against daily NEE, Qle, Qh
using TransformUnscented Kalman Inversion with ObservationSeries (minibatched by
year) and WorkerBackend + SlurmManager for parallel ensemble members.

10 years of observations (2004-2013), minibatch_size=2 → 5 batches/epoch.

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
const DT = Float64(450)

const climaland_dir = pkgdir(ClimaLand)
const OUTPUT_DIR = joinpath(climaland_dir, "experiments/calibrate_dk_sor/output")
const OBS_FILEPATH =
    joinpath(climaland_dir, "experiments/calibrate_dk_sor/observations.jld2")

# ── Priors ───────────────────────────────────────────────────────────────────
# Prior names MUST match ClimaParams TOML keys, since ClimaCalibrate writes
# parameter TOMLs using these names and LP.create_toml_dict reads them.

priors = [
    # Canopy parameters (names already match TOML keys)
    PD.constrained_gaussian("moisture_stress_c", 0.5, 0.3, 0.01, 5.0),
    PD.constrained_gaussian("pmodel_cstar", 0.43, 0.15, 0.05, 2.0),
    PD.constrained_gaussian("pmodel_β", 51.0, 20.0, 5.0, 500.0),
    PD.constrained_gaussian("leaf_Cd", 0.1, 0.05, 0.005, 1.0),
    PD.constrained_gaussian("canopy_z_0m_coeff", 0.05, 0.03, 0.001, 0.3),
    PD.constrained_gaussian("canopy_z_0b_coeff", 0.001, 0.0005, 1e-5, 0.01),
    PD.constrained_gaussian("canopy_d_coeff", 0.1, 0.05, 0.001, 0.95),
    PD.constrained_gaussian("canopy_K_lw", 0.85, 0.25, 0.1, 2.0),
    PD.constrained_gaussian("canopy_emissivity", 0.97, 0.02, 0.9, 1.0),
    # DAMM soilCO2 parameters (renamed to match TOML keys)
    PD.constrained_gaussian(
        "soilCO2_pre_exponential_factor",
        25000.0,
        10000.0,
        1000.0,
        200000.0,
    ),
    PD.constrained_gaussian("michaelis_constant", 0.01, 0.005, 1e-4, 0.1),
    PD.constrained_gaussian("O2_michaelis_constant", 0.01, 0.005, 1e-4, 0.1),
]
prior = PD.combine_distributions(priors)

# ── Load Observations & Build ObservationSeries ──────────────────────────────

obs_data = JLD2.load(OBS_FILEPATH)
observation_vector = obs_data["observation_vector"]
cal_years = obs_data["cal_years"]

n_samples = length(observation_vector)
println("Loaded $n_samples yearly observations ($(first(cal_years))-$(last(cal_years)))")

obs_series = EKP.ObservationSeries(
    Dict(
        "observations" => observation_vector,
        "names" => [string(yr) for yr in cal_years],
        "minibatcher" => ClimaCalibrate.minibatcher_over_samples(
            n_samples,
            MINIBATCH_SIZE,
        ),
    ),
)

# ── Create EKP ───────────────────────────────────────────────────────────────

rng = Random.MersenneTwister(1234)
ekp = EKP.EnsembleKalmanProcess(
    obs_series,
    EKP.TransformUnscented(prior; impose_prior = true),
    verbose = true,
    rng = rng,
)

N_ens = EKP.get_N_ens(ekp)
println("Ensemble size: $N_ens (for $(length(priors)) parameters)")
println("Minibatch size: $MINIBATCH_SIZE years → $(div(n_samples, MINIBATCH_SIZE)) batches/epoch")

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
    joinpath(pkgdir(ClimaLand), "experiments/calibrate_dk_sor/model_interface.jl"),
)

# ── Run Calibration ──────────────────────────────────────────────────────────

println("\n=== Starting ClimaCalibrate calibration ===")
println("  Backend: WorkerBackend")
println("  Iterations: $N_ITERATIONS")
println("  Ensemble size: $N_ens")
println("  Calibration years: $(first(cal_years))-$(last(cal_years))")
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
