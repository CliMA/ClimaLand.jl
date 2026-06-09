"""
ClimaCalibrate driver for NEON site soil CO₂ calibration.

Calibrates 4 DAMM soil CO₂ parameters (soilCO2_reference_rate,
soilCO2_activation_energy, michaelis_constant, O2_michaelis_constant) against daily soil CO₂ concentration
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
SITE_ID = get(ENV, "NEON_SITE_ID", "NEON-srer")
N_ITERATIONS = parse(Int, get(ENV, "NEON_N_ITERATIONS", "10"))
SPINUP_DAYS = parse(Int, get(ENV, "NEON_SPINUP_DAYS", "20"))
outputpath = get(ENV, "CALL_OUTPUT_PATH", "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_siteruns/$(SITE_ID)/")

start_date =
    DateTime(get(ENV, "NEON_START_DATE", string(Date(2009,1,1))))
stop_date = DateTime(get(ENV, "NEON_STOP_DATE", string(Date(2009,12,31))))

DT = Float64(900)#450)
#SITE_ID = "NEON-cper"
climaland_dir = pkgdir(ClimaLand)
output_base = joinpath(outputpath, "output")
OUTPUT_DIR = find_available_dir(output_base)
ENV["CALL_OUTPUT_1"] = OUTPUT_DIR
PRIOR_VALUES_FILE = joinpath(OUTPUT_DIR, "prior_values.txt")
FINAL_PARAMS_FILE = joinpath(OUTPUT_DIR, "final_parameter_means.txt")
OBS_FILEPATH = joinpath(outputpath, "observations.jld2")
    #"/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_siteruns/$(SITE_ID)/$(SITE_ID)_$(Date(start_date))_$(Date(stop_date))_SpinUp$(SPINUP_DAYS)/observations.jld2"
    #joinpath(climaland_dir, "experiments/calibrate_neon/observations.jld2")
Caldepthnum = get(ENV, "CALL_DEPTH", "0.00")


# ── Priors ───────────────────────────────────────────────────────────────────
# Prior names MUST match ClimaParams TOML keys, since ClimaCalibrate writes
# parameter TOMLs using these names

# Prior format: [mean, std, lower_bound, upper_bound]
#
# Davidson et al. (2011) DAMM defaults (ClimaParams) for reference:
#   soilCO2_pre_exponential_factor = 194e3  (V_ref_sx; not directly used here)
#   soilCO2_activation_energy      = 61e3   (J/mol)
#   michaelis_constant             = 5e-3   (kg C / m³)
#   O2_michaelis_constant          = 4e-3   (m³ / m³)
#
# Rationale for the values below:
# - soilCO2_reference_rate (V_ref_sx): controls absolute respiration magnitude.
#   Centered at the previously calibrated NEON value (~5e-9) but with broader
#   spread so the calibration can move up/down by an order of magnitude.
# - soilCO2_activation_energy (Ea_sx): centered near Davidson default 61 kJ/mol
#   with moderate spread covering typical Q10 behavior across sites.
# - michaelis_constant (kM_sx): centered at Davidson default 5e-3. At this value
#   MM_sx stays substrate-limited so the θ_l^3 pulse signal propagates into
#   respiration. Previous prior (3e-4) saturated MM_sx and damped the pulses.
# - O2_michaelis_constant (kM_o2): centered at Davidson default 4e-3 so O2
#   limitation only kicks in for genuinely wet conditions.
soilCO2_reference_rate    = [1.0e-7, 5.0e-8, 1.0e-12, 1.0e-6]
soilCO2_activation_energy = [110550.0, 40000.0, 10.0, 200000.0]
michaelis_constant        = [0.04, 2.0e-2, 0.0, Inf]
O2_michaelis_constant     = [0.04, 2.0e-2, 0.0, Inf]

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
]
prior = PD.combine_distributions(priors)

# Create output directory (find_available_dir ensures unique path)
mkpath(OUTPUT_DIR)

# Snapshot model scripts IMMEDIATELY — before any include/load of the run-time
# scripts and before the long calibration loop, so the saved copy reflects
# the code that this run actually uses (not edits made afterwards while the
# run is still going).
scripts_src = joinpath(climaland_dir, "experiments/calibrate_neon")
scripts_dst = joinpath(OUTPUT_DIR, "model_scripts")
cp(scripts_src, scripts_dst; force=true)
scripts_src = joinpath(climaland_dir, "src/standalone/Soil/Biogeochemistry")
cp(scripts_src, joinpath(scripts_dst, "Biogeochemistry"); force=true)

open(PRIOR_VALUES_FILE, "w") do io
    println(io, "soilCO2_reference_rate = $(soilCO2_reference_rate)")
    println(io, "soilCO2_activation_energy = $(soilCO2_activation_energy)")  # neu
    println(io, "michaelis_constant = $(michaelis_constant)")
    println(io, "O2_michaelis_constant = $(O2_michaelis_constant)")
    println(io, "Means only:")
    println(io, "soilCO2_reference_rate = $(soilCO2_reference_rate[1])")
    println(io, "soilCO2_activation_energy = $(soilCO2_activation_energy[1])")  # neu
    println(io, "michaelis_constant = $(michaelis_constant[1])")
    println(io, "O2_michaelis_constant = $(O2_michaelis_constant[1])")
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
# Globals (no `const`) so the same Julia session can re-`include` this script
# with different values across iterations of a multi-run pipeline.
@everywhere SITE_ID = $SITE_ID
@everywhere OUTPUT_DIR = $OUTPUT_DIR
@everywhere OBS_FILEPATH = $OBS_FILEPATH
@everywhere DT = $DT
@everywhere Caldepthnum = $Caldepthnum

@everywhere include(    joinpath(pkgdir(ClimaLand), "experiments/calibrate_neon/site_metadata.jl"),
)
@everywhere include(    joinpath(pkgdir(ClimaLand), "experiments/calibrate_neon/model_interface.jl"),
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
    ClimaCalibrate.WorkerBackend(),
    ekp,
    NeonModelInterface(),
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

let
    iter_dirs = filter(d -> isdir(joinpath(OUTPUT_DIR, d)) && startswith(d, "iteration_"),
                       readdir(OUTPUT_DIR))
    if !isempty(iter_dirs)
        last_iter = maximum(sort(iter_dirs))
        eki_path = joinpath(OUTPUT_DIR, last_iter, "eki_file.jld2")
        if isfile(eki_path)
            ENV["CALL_EKI_PATH"] = eki_path
            println("Set CALL_EKI_PATH = $eki_path")
        end
    end
end

rmprocs(workers())