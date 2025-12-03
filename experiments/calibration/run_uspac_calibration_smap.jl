using Dates
using Distributed
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaComms
import ClimaLand
import EnsembleKalmanProcesses as EKP
import JLD2

#include(joinpath(pkgdir(ClimaLand), "experiments/calibration/api.jl"))
include(joinpath(@__DIR__, "api.jl"))

#NOTE: must run generate_observations.jl AND generate_smap_observations.jl FIRST (slurm file already does this)

# Generate seasonal date ranges programmatically
function generate_seasonal_date_ranges(start_year, end_year)
    ranges = []
    for year in start_year:end_year
        if year == 2015
            # SMAP starts April 2015
            push!(ranges, ("2015-06-01", "2015-06-01"))  # Summer
            push!(ranges, ("2015-09-01", "2015-09-01"))  # Fall
            push!(ranges, ("2015-12-01", "2015-12-01"))  # Winter
        else
            push!(ranges, ("$year-03-01", "$year-03-01"))  # Spring
            push!(ranges, ("$year-06-01", "$year-06-01"))  # Summer
            push!(ranges, ("$year-09-01", "$year-09-01"))  # Fall
            push!(ranges, ("$year-12-01", "$year-12-01"))  # Winter
        end
    end
    return ranges
end

# -------------------- CALIBRATION CONFIG --------------------
# For JOINT calibration against LH flux AND SMAP soil moisture
const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["lhf", "sm_surface"],  # ← Both variables
    minibatch_size = 35,
    n_iterations = 10,
    sample_date_ranges = generate_seasonal_date_ranges(2015, 2023),
    extend = Dates.Month(0),
    spinup = Dates.Month(3),
    nelements = (101, 15),
    output_dir = "experiments/calibration/land_model_joint",  # ← Joint calibration output
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector_joint.jld2",  # ← Combined obs
)

# For latent heat flux only (original):
# const CALIBRATE_CONFIG = CalibrateConfig(;
#     short_names = ["lhf"],
#     obs_vec_filepath = "experiments/calibration/land_observation_vector.jld2",
#     output_dir = "experiments/calibration/land_model",
# )

# For SMAP soil moisture only:
# const CALIBRATE_CONFIG = CalibrateConfig(;
#     short_names = ["sm_surface"],
#     obs_vec_filepath = "experiments/calibration/land_observation_vector_SMAP.jld2",
#     output_dir = "experiments/calibration/land_model_smap",
# )

if abspath(PROGRAM_FILE) == @__FILE__
    # -------------------- PRIORS: 8 uSPAC-Π parameters --------------------
    # θ = [alpha_R, beta_R, alpha_F, beta_F, alpha_T, beta_Ts, alpha_S, beta_Ss]
    # Changed prior means from 0.0 to 0.1 to avoid division-by-zero in stomatalconductance.jl
    # When all π-groups = 0, the model has mathematical singularities that produce zero conductance
    priors = EKP.ParameterDistribution[
        EKP.constrained_gaussian("alpha_R", 0.1, 0.30, -Inf, Inf),
        EKP.constrained_gaussian("beta_R",  0.1, 0.30, -Inf, Inf),
        EKP.constrained_gaussian("alpha_F", 0.1, 0.30, -Inf, Inf),
        EKP.constrained_gaussian("beta_F",  0.1, 0.30, -Inf, Inf),
        EKP.constrained_gaussian("alpha_T", 0.1, 0.30, -Inf, Inf),
        EKP.constrained_gaussian("beta_Ts", 0.1, 0.30, -Inf, Inf),
        EKP.constrained_gaussian("alpha_S", 0.1, 0.30, -Inf, Inf),
        EKP.constrained_gaussian("beta_Ss", 0.1, 0.30, -Inf, Inf),
    ]
    prior = EKP.combine_distributions(priors)

    # -------------------- OBSERVATIONS --------------------
    obs_data = JLD2.load(CALIBRATE_CONFIG.obs_vec_filepath)
    
    # **JOINT OBSERVATION VECTOR**
    # Check if this is a joint calibration file (has both lhf_vector and sm_vector)
    if "lhf_vector_scaled" in keys(obs_data) && "sm_vector_scaled" in keys(obs_data)
        @info "Loading JOINT observation vector (LH + SMAP) with pre-applied scaling"
        
        # Use the SCALED versions (this is what we want for calibration)
        lhf_obs = obs_data["lhf_vector_scaled"]
        sm_obs = obs_data["sm_vector_scaled"]
        
        # The observation_vector in the file is already concatenated and scaled
        observation_vector = obs_data["observation_vector"]
        
        @info "Joint observations loaded" n_lhf=length(lhf_obs) n_sm=length(sm_obs) total=length(observation_vector)
        
        # Print weighting info
        if haskey(obs_data, "lhf_weight") && haskey(obs_data, "sm_weight")
            @info "Observation weighting" lhf_weight=obs_data["lhf_weight"] sm_weight=obs_data["sm_weight"] normalize_by_variance=obs_data["normalize_by_variance"]
        end
        
        # Print SMAP coverage info
        if "sm_metadata" in keys(obs_data)
            @info "SMAP observation metadata:"
            for (i, meta) in enumerate(obs_data["sm_metadata"])
                println("  Period $i: $(meta["start_date"]) to $(meta["stop_date"])")
                println("    Observations: $(meta["n_observations"])")
                if haskey(meta, "coverage_stats")
                    println("    Coverage: $(meta["coverage_stats"]["land_coverage_pct"])%")
                    if haskey(meta["coverage_stats"], "avg_model_points_per_smap_pixel")
                        println("    Avg model pts/SMAP pixel: $(meta["coverage_stats"]["avg_model_points_per_smap_pixel"])")
                    end
                end
            end
        end
    else
        # Fallback for single-variable files (SMAP-only or LH-only)
        @info "Loading single-variable observation vector"
        observation_vector = obs_data["observation_vector"]
        
        # Print SMAP info if available
        if "metadata" in keys(obs_data) && !isempty(obs_data["metadata"])
            meta = obs_data["metadata"]
            if isa(meta, AbstractVector)
                @info "SMAP-only calibration" n_periods=length(meta)
                for (i, m) in enumerate(meta)
                    if haskey(m, "n_observations")
                        println("  Period $i: $(m["n_observations"]) observations")
                    end
                end
            end
        end
    end

    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    minibatch_size     = CALIBRATE_CONFIG.minibatch_size

    obs_series = EKP.ObservationSeries(Dict(
        "observations" => observation_vector,
        "names" => [string(Dates.year(start_date)) for (start_date, _) in sample_date_ranges],
        "minibatcher" => ClimaCalibrate.minibatcher_over_samples(
            length(observation_vector), minibatch_size),
    ))

    # -------------------- EKP SETUP --------------------
    rng_seed = CALIBRATE_CONFIG.rng_seed
    rng = Random.MersenneTwister(rng_seed)

    # Note: You should check that the ensemble size is the same as the number of
    # tasks in the batch script
    # For example, if you are calibrating 3 parameters and are using
    # EKP.TransformUnscented, then the number of tasks should be 7, since
    # 3 * 2 + 1 = 7 # so, 8 * 2 + 1 = 17
    # TransformUnscented → N_ens = 2P + 1 = 17
    ekp = EKP.EnsembleKalmanProcess(
        obs_series,
        EKP.TransformUnscented(prior; impose_prior=true),
        verbose = true,
        rng     = rng,
        scheduler = EKP.DataMisfitController(terminate_at = 100),
    )

    # If you prefer a larger stochastic ensemble, swap for:
    # ekp = EKP.EnsembleKalmanProcess(
    #     obs_series, EKP.Inversion(); rng, scheduler = EKP.DataMisfitController(terminate_at=100)
    # )
    # and set N_ens via initial ensemble size in ClimaCalibrate.calibrate (handled internally).

    # -------------------- BACKEND / WORKERS --------------------
    curr_backend = ClimaCalibrate.get_backend()
    N_ens = EKP.get_N_ens(ekp)
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
