"""
ClimaCalibrate driver for NEON site soil CO₂ calibration.

Calibrates 3 DAMM soil CO₂ parameters (soilCO2_pre_exponential_factor,
michaelis_constant, O2_michaelis_constant) against daily soil CO₂ concentration
observations at NEON sites using Unscented Kalman Inversion.

Configuration via environment variables:
    NEON_SITE_ID      — NEON site ID (default: "NEON-srer")
    NEON_SPINUP_DAYS  — Number of spinup days (default: 20)
    NEON_N_ITERATIONS — Number of UKI iterations (default: 10)

Usage:
    julia --project=.buildkite experiments/calibrate_neon/run_calibration.jl
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
outdir = "/kiwi-data/Data/groupMembers/evametz/Neon/Neon_data/dataframes_Neon/outputrun"
const SITE_ID = get(ENV, "NEON_SITE_ID", "NEON-srer")
const N_ITERATIONS = parse(Int, get(ENV, "NEON_N_ITERATIONS", "10"))
const DT = Float64(450)
SITE_ID = "NEON-cper"
const climaland_dir = pkgdir(ClimaLand)
const OUTPUT_DIR = joinpath(outdir, "experiments/calibrate_neon/output")
const OBS_FILEPATH =
    joinpath(climaland_dir, "experiments/calibrate_neon/observations.jld2")

# ── Priors ───────────────────────────────────────────────────────────────────
# Prior names MUST match ClimaParams TOML keys, since ClimaCalibrate writes
# parameter TOMLs using these names and LP.create_toml_dict reads them.

priors = [
    PD.constrained_gaussian(
        "soilCO2_pre_exponential_factor",
        2000.0,
        1000.0,
        100.0,
        20000.0,
    ),
    PD.constrained_gaussian("michaelis_constant", 0.046, 0.020, 1e-5, 0.1),
    PD.constrained_gaussian("O2_michaelis_constant", 0.066, 0.03, 1e-5, 0.12),
]
prior = PD.combine_distributions(priors)

# ── Load Observations ────────────────────────────────────────────────────────

obs_data = JLD2.load(OBS_FILEPATH)
y_obs = obs_data["y_obs"]
noise_cov = obs_data["noise_cov"]
obs_dates = obs_data["obs_dates"]

n_obs = length(obs_dates)
println("Loaded $n_obs valid observation days for site $SITE_ID")
println("Observation vector length: $n_obs")

# ── Create EKP ───────────────────────────────────────────────────────────────

rng = Random.MersenneTwister(1234)
ekp = EKP.EnsembleKalmanProcess(
    y_obs,
    noise_cov,
    EKP.Unscented(prior; α_reg = 1.0, update_freq = 1);
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
    joinpath(pkgdir(ClimaLand), "experiments/calibrate_neon/model_interface.jl"),
)

# ── Run Calibration ──────────────────────────────────────────────────────────

println("\n=== Starting ClimaCalibrate calibration ===")
println("  Backend: WorkerBackend")
println("  Site: $SITE_ID")
println("  Iterations: $N_ITERATIONS")
println("  Ensemble size: $N_ens")
println("  Observation days: $n_obs")
println("  Parameters: $(length(priors)) (3 DAMM soilCO2)")
println("  LAI: MODIS")
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
