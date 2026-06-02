"""
ClimaCalibrate driver for NEON site soil CO₂ calibration with porosity scaling.

Calibrates 4 parameters (soilCO2_pre_exponential_factor, michaelis_constant,
O2_michaelis_constant, nu_scaling_factor) against daily soil CO₂ concentration
observations at NEON sites using Unscented Kalman Inversion.

nu_scaling_factor is a scalar multiplier applied to the spatially-varying soil
porosity profile (ν) loaded from Gupta et al. 2020. It is clamped so that
all porosity values remain in [0.05, 0.95].

Configuration via environment variables:
    NEON_SITE_ID      — NEON site ID (default: "NEON-srer")
    NEON_SPINUP_DAYS  — Number of spinup days (default: 20)
    NEON_N_ITERATIONS — Number of UKI iterations (default: 10)

Usage:
    julia --project=.buildkite experiments/calibrate_neon/run_calibration_wPorosity.jl
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

# ── Helper Functions ────────────────────────────────────────────────────────

function find_available_dir(base_dir)
    """Find an available directory by appending _1, _2, etc. if base_dir exists."""
    if !isdir(base_dir)
        return base_dir
    end
    i = 1
    while isdir("$(base_dir)_$i")
        i += 1
    end
    return "$(base_dir)_$i"
end

# ── Configuration ────────────────────────────────────────────────────────────
const SITE_ID = get(ENV, "NEON_SITE_ID", "NEON-srer")
const N_ITERATIONS = parse(Int, get(ENV, "NEON_N_ITERATIONS", "10"))
const SPINUP_DAYS = parse(Int, get(ENV, "NEON_SPINUP_DAYS", "20"))
const outputpath = get(ENV, "CALL_OUTPUT_PATH", "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_siteruns/$(SITE_ID)/")

start_date =
    DateTime(get(ENV, "NEON_START_DATE", string(Date(2009,1,1))))
stop_date = DateTime(get(ENV, "NEON_STOP_DATE", string(Date(2009,12,31))))

const DT = Float64(180)
const climaland_dir = pkgdir(ClimaLand)
output_base = joinpath(outputpath, "output_wPorosity")
const OUTPUT_DIR = find_available_dir(output_base)
const PRIOR_VALUES_FILE = joinpath(OUTPUT_DIR, "prior_values.txt")
const FINAL_PARAMS_FILE = joinpath(OUTPUT_DIR, "final_parameter_means.txt")
const OBS_FILEPATH = joinpath(outputpath, "observations.jld2")
const Caldepthnum = get(ENV, "CALL_DEPTH", "0.00")

# ── Priors ───────────────────────────────────────────────────────────────────
# Prior names MUST match ClimaParams TOML keys for ClimaParams-registered
# parameters. nu_scaling_factor is NOT a ClimaParams key — it is read
# directly from the TOML via Julia's TOML parser in model_interface_wPorosity.jl.

soilCO2_reference_rate = [2.526e-7, 1.0e-7, 5.0e-8, 1.0e-6]
soilCO2_activation_energy = [40000.0, 15000.0, 20000.0, 100000.0,]
michaelis_constant = [0.3, 0.2, 0.01, 1.0]
O2_michaelis_constant = [0.005, 0.003, 5.0e-4, 5.0e-2]
nu_scaling_factor = [1.0, 0.15, 0.5, 1.5]

priors = [
    PD.constrained_gaussian(
        "soilCO2_reference_rate",
        soilCO2_reference_rate[1],
        soilCO2_reference_rate[2],
        soilCO2_reference_rate[3],
        soilCO2_reference_rate[4],
    ),
    PD.constrained_gaussian(
        "soilCO2_activation_energy",
        soilCO2_activation_energy[1],
        soilCO2_activation_energy[2],
        soilCO2_activation_energy[3],
        soilCO2_activation_energy[4],
    ),
    PD.constrained_gaussian(
        "michaelis_constant",
        michaelis_constant[1],
        michaelis_constant[2],
        michaelis_constant[3],
        michaelis_constant[4],
    ),
    PD.constrained_gaussian(
        "O2_michaelis_constant",
        O2_michaelis_constant[1],
        O2_michaelis_constant[2],
        O2_michaelis_constant[3],
        O2_michaelis_constant[4],
    ),
    PD.constrained_gaussian(
        "nu_scaling_factor",
        nu_scaling_factor[1],
        nu_scaling_factor[2],
        nu_scaling_factor[3],
        nu_scaling_factor[4],
    ),
]
prior = PD.combine_distributions(priors)

# Create output directory
mkpath(OUTPUT_DIR)
open(PRIOR_VALUES_FILE, "w") do io
    println(io, "soilCO2_reference_rate = $(soilCO2_reference_rate)")
    println(io, "soilCO2_activation_energy = $(soilCO2_activation_energy)")  # neu
    println(io, "michaelis_constant = $(michaelis_constant)")
    println(io, "O2_michaelis_constant = $(O2_michaelis_constant)")
    println(io, "nu_scaling_factor = $(nu_scaling_factor)")
    println(io, "Means only:")
    println(io, "soilCO2_reference_rate = $(soilCO2_reference_rate[1])")
    println(io, "soilCO2_activation_energy = $(soilCO2_activation_energy[1])")  # neu
    println(io, "michaelis_constant = $(michaelis_constant[1])")
    println(io, "O2_michaelis_constant = $(O2_michaelis_constant[1])")
    println(io, "nu_scaling_factor = $(nu_scaling_factor[1])")
end
scripts_src = joinpath(climaland_dir, "experiments/calibrate_neon")
scripts_dst = joinpath(OUTPUT_DIR, "model_scripts")
cp(scripts_src, scripts_dst; force=true)
scripts_src = joinpath(climaland_dir, "src/standalone/Soil/Biogeochemistry")
cp(scripts_src, joinpath(scripts_dst, "Biogeochemistry"); force=true)

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

# ── Add Workers ──────────────────────────────────────────────────────────────

curr_backend = ClimaCalibrate.get_backend()
if curr_backend == ClimaCalibrate.CaltechHPCBackend ||
   curr_backend == ClimaCalibrate.GCPBackend
    @info "Check your slurm script that the number of tasks is the same as the ensemble size ($N_ens)"
    addprocs(ClimaCalibrate.SlurmManager())
else
    current_workers = nprocs() == 1 ? 0 : nworkers()
    n_add = max(0, N_ens - current_workers)
    if n_add > 0
        @info "No SLURM detected — adding $n_add local workers for $N_ens ensemble members"
        addprocs(n_add)
    end
end

# ── Broadcast Configuration to Workers ───────────────────────────────────────

@everywhere using Distributed
@everywhere import ClimaLand
@everywhere const SITE_ID = $SITE_ID
@everywhere const OUTPUT_DIR = $OUTPUT_DIR
@everywhere const OBS_FILEPATH = $OBS_FILEPATH
@everywhere const DT = $DT
@everywhere const Caldepthnum = $Caldepthnum

@everywhere include(    joinpath(pkgdir(ClimaLand), "experiments/calibrate_neon/site_metadata.jl"),
)
@everywhere include(    joinpath(pkgdir(ClimaLand), "experiments/calibrate_neon/model_interface_wPorosity.jl"),
)

# ── Run Calibration ──────────────────────────────────────────────────────────

println("\n=== Starting ClimaCalibrate calibration (with porosity scaling) ===")
println("  Backend: WorkerBackend")
println("  Site: $SITE_ID")
println("  Iterations: $N_ITERATIONS")
println("  Ensemble size: $N_ens")
println("  Observation days: $n_obs")
println("  Parameters: $(length(priors)) (4 DAMM soilCO2 + nu_scaling_factor)")
println("  LAI: MODIS")
println("  Output: $OUTPUT_DIR")

eki = ClimaCalibrate.calibrate(
    ClimaCalibrate.WorkerBackend(),
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

open(FINAL_PARAMS_FILE, "w") do io
    println(io, "=== Calibration Complete ===")
    println(io, "Site: $SITE_ID")
    println(io, "Final parameter means:")
    for (name, val) in zip(param_names, final_params)
        println(io, "  $name = $(round(val, sigdigits = 5))")
    end
end
println("Saved final parameter means to $FINAL_PARAMS_FILE")

rmprocs(workers())
