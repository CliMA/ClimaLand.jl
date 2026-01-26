using Dates
using Distributed
using LinearAlgebra: Diagonal
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaComms
import ClimaLand
import EnsembleKalmanProcesses as EKP
import JLD2

#include(joinpath(pkgdir(ClimaLand), "experiments/calibration/api.jl"))
include(joinpath(@__DIR__, "api.jl"))

#NOTE: must run generate_observations.jl FIRST (slurm file already does this)

# -------------------- CALIBRATION CONFIG --------------------
# Using SMAP soil moisture observations only
const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["sm_surface"],  # SMAP soil moisture only
    minibatch_size = 1,
    n_iterations = 3,
    # Time period with SMAP data coverage
    sample_date_ranges = [("2016-12-1", "2019-9-1")],
    extend = Dates.Month(3),
    spinup = Dates.Month(3),
    # Grid (horizontal elements, vertical levels)
    nelements = (101, 15),
    output_dir = "experiments/calibration/land_model",
    rng_seed = 42,
    # SMAP observation file (same format as ERA5)
    obs_vec_filepath = "experiments/calibration/land_observation_vector.jld2",
)

if abspath(PROGRAM_FILE) == @__FILE__
    # -------------------- PRIORS: 6 uSPAC trait parameters --------------------
    # θ = [βkx_base, βkx_coord, βψx50_base, βψx50_slope, βΠR_base, βΠR_slope]
    # Recentered on iteration_002 mean values (best performing iteration)
    # Transformations ensure physically realistic ranges:
    # - kx: log-normal → absolute conductivity (m day⁻¹ MPa⁻¹)
    # - ψx50: exponential → negative water tension (MPa)
    # - ΠR: logistic → [0, 1] (isohydric to anisohydric)
    priors = EKP.ParameterDistribution[  
        EKP.constrained_gaussian("βkx_base", -3.477, 0.45, -10.0, -2.0),     # Baseline conductivity
        EKP.constrained_gaussian("βkx_coord", -0.014, 0.24, -2.0, 0.5),      # kx-P50 coordination
        EKP.constrained_gaussian("βψx50_base", 1.540, 0.64, -4.0, 4.0),      # ψx50 baseline
        EKP.constrained_gaussian("βψx50_slope", 1.525, 0.51, 0.0, 10.0),     # ψx50 slope
        EKP.constrained_gaussian("βΠR_base", -5.034, 0.71, -8.0, -2.0),      # ΠR baseline
        EKP.constrained_gaussian("βΠR_slope", 6.795, 1.21, 2.0, 15.0)        # ΠR slope (isohydric-anisohydric transition)
    ]
    prior = EKP.combine_distributions(priors)

    # -------------------- OBSERVATIONS --------------------
    # Load observation vector
    # Check if it's the new joint format (raw samples) or old format (EKP.Observation)
    obs_filepath = CALIBRATE_CONFIG.obs_vec_filepath
    
    if endswith(obs_filepath, "joint.jld2")
        # New format: raw scaled samples, create EKP.Observation on-the-fly
        include(joinpath(@__DIR__, "load_joint_observations.jl"))
        observation_vector = load_joint_observations(obs_filepath)
    else
        # SMAP format: file contains multiple objects, extract observation_vector
        data = JLD2.load(obs_filepath)
        observation_vector = data["observation_vector"]
    end
    
    @info "Loaded observation vector" filepath=obs_filepath n_periods=length(observation_vector)

    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    minibatch_size = CALIBRATE_CONFIG.minibatch_size
    
    # Create ObservationSeries matching run_calibration.jl pattern
    obs_series = EKP.ObservationSeries(
        Dict(
            "observations" => observation_vector,
            "names" => [
                string(Dates.year(DateTime(start_date))) for
                (start_date, stop_date) in sample_date_ranges
            ],
            "minibatcher" => ClimaCalibrate.minibatcher_over_samples(
                length(observation_vector),
                minibatch_size,
            ),
        ),
    )

    # -------------------- EKP SETUP --------------------
    rng_seed = CALIBRATE_CONFIG.rng_seed
    rng = Random.MersenneTwister(rng_seed)

    # Note: You should check that the ensemble size is the same as the number of
    # tasks in the batch script
    # For example, if you are calibrating 3 parameters and are using
    # EKP.TransformUnscented, then the number of tasks should be 7, since
    # 3 * 2 + 1 = 7
    # TransformUnscented → N_ens = 2P + 1 = 2*6 + 1 = 13
    ekp = EKP.EnsembleKalmanProcess(
        obs_series,
        EKP.TransformUnscented(prior; impose_prior=true),
        verbose = true,
        rng     = rng,
        scheduler = EKP.DataMisfitController(terminate_at = 100),
    )

    # Update ensemble size for TransformUnscented: N_ens = 2P + 1 = 2*6 + 1 = 13
    N_ens = EKP.get_N_ens(ekp)  # Should be 13 now
    
    # Detect the backend based on environment/hostname
    hostname = gethostname()
    curr_backend = if occursin("derecho", hostname) || occursin("deg", hostname)
        ClimaCalibrate.DerechoBackend
    elseif haskey(ENV, "SLURM_JOB_ID")
        ClimaCalibrate.CaltechHPCBackend
    else
        ClimaCalibrate.GCPBackend
    end
    
    if curr_backend == ClimaCalibrate.DerechoBackend
        addprocs(
            ClimaCalibrate.PBSManager(N_ens),
            q = "main", A = "UCIT0011",
            l_select = "1:ngpus=1:ncpus=4",
            l_walltime = "11:30:00",
        )
    elseif (curr_backend == ClimaCalibrate.CaltechHPCBackend) ||
           (curr_backend == ClimaCalibrate.GCPBackend)
        @info "Ensure your Slurm job launches exactly $N_ens tasks (one per ensemble member)."
        addprocs(ClimaCalibrate.SlurmManager())
    end

    # -------------------- DISTRIBUTE MAPPERS --------------------
    include(joinpath(pkgdir(ClimaLand), "experiments/calibration/observation_map.jl"))

    @everywhere import ClimaLand
    @everywhere experiment_dir = joinpath(pkgdir(ClimaLand), "experiments")
    @everywhere include(joinpath(pkgdir(ClimaLand), "experiments/calibration/api.jl"))
    @everywhere CALIBRATE_CONFIG = $CALIBRATE_CONFIG
    @everywhere include(joinpath(experiment_dir, "calibration", "model_interface.jl"))

    # -------------------- RUN CALIBRATION --------------------
    
    # Print configuration summary BEFORE running (not after!)
    @info "Calibration Configuration:" CALIBRATE_CONFIG.short_names CALIBRATE_CONFIG.n_iterations CALIBRATE_CONFIG.minibatch_size
    @info "Sample date ranges ($(length(CALIBRATE_CONFIG.sample_date_ranges)) periods):"
    for (i, (start, stop)) in enumerate(CALIBRATE_CONFIG.sample_date_ranges)
        println("  $i: $start to $stop")
    end
    
    eki = ClimaCalibrate.calibrate(
        ClimaCalibrate.WorkerBackend,
        ekp,
        CALIBRATE_CONFIG.n_iterations,
        prior,
        CALIBRATE_CONFIG.output_dir,
    )

end
