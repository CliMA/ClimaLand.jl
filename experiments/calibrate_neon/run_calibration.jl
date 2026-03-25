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
start_date =
    DateTime(get(ENV, "NEON_START_DATE", string(Date(2009,1,1))))
stop_date = DateTime(get(ENV, "NEON_STOP_DATE", string(Date(2009,12,31))))

const DT = Float64(450)
#SITE_ID = "NEON-cper"
const climaland_dir = pkgdir(ClimaLand)
outdir = "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_siteruns/$(SITE_ID)/$(SITE_ID)_$(Date(start_date))_$(Date(stop_date))_SpinUp$(SPINUP_DAYS)/"
output_base = joinpath(outdir, "output")
const OUTPUT_DIR = find_available_dir(output_base)
const PRIOR_VALUES_FILE = joinpath(OUTPUT_DIR, "prior_values.txt")
const OBS_FILEPATH =
    "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_siteruns/$(SITE_ID)/$(SITE_ID)_$(Date(start_date))_$(Date(stop_date))_SpinUp$(SPINUP_DAYS)/observations.jld2"
    #joinpath(climaland_dir, "experiments/calibrate_neon/observations.jld2")


# ── Priors ───────────────────────────────────────────────────────────────────
# Prior names MUST match ClimaParams TOML keys, since ClimaCalibrate writes
# parameter TOMLs using these names and LP.create_toml_dict reads them.

preexp = [2000.0, 1000.0, 100.0, 20000.0]
michaelis = [0.046, 0.020, 1e-5, 0.1]
o2_michaelis = [0.066, 0.03, 1e-5, 0.12]

priors = [
    PD.constrained_gaussian(
        "soilCO2_pre_exponential_factor",
        preexp[1],
        preexp[2],
        preexp[3],
        preexp[4],
    ),
    PD.constrained_gaussian(
        "michaelis_constant",
        michaelis[1],
        michaelis[2],
        michaelis[3],
        michaelis[4],
    ),
    PD.constrained_gaussian(
        "O2_michaelis_constant",
        o2_michaelis[1],
        o2_michaelis[2],
        o2_michaelis[3],
        o2_michaelis[4],
    ),
]
prior = PD.combine_distributions(priors)

# Create output directory (find_available_dir ensures unique path)
mkpath(OUTPUT_DIR)
open(PRIOR_VALUES_FILE, "w") do io
    println(io, "soilCO2_pre_exponential_factor = $(preexp)")
    println(io, "michaelis_constant = $(michaelis)")
    println(io, "O2_michaelis_constant = $(o2_michaelis)")
end

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
#elseif nworkers() == 1
#    @info "No SLURM detected and only 1 worker — running serially (JuliaBackend fallback)"
else
    # nworkers() returns 1 when nprocs()==1 (main process counted as worker),
    # but after addprocs() the main process stops counting as a worker.
    # Always add N_ens workers so every ensemble member runs in parallel.
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
