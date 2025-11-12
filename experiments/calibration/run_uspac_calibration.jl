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

#NOTE: must run generate_observations.jl FIRST (slurm file already does this)

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
# Change "short_names" to the model diagnostic(s) you want to match.
# For global LE calibration, use "lhf" (total latent heat flux).
const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["lhf"],
    minibatch_size = 35,
    n_iterations = 10,
    sample_date_ranges = generate_seasonal_date_ranges(2015, 2023),
    extend = Dates.Month(0),
    spinup = Dates.Month(3),
    nelements = (101, 15),  # ← One 'e' not two!
    output_dir = "experiments/calibration/land_model",
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector.jld2",
)

# Later: add SMAP SWC later
# const CALIBRATE_CONFIG = CalibrateConfig(;
#     short_names = ["sm_surface"],   # <<< SMAP variable key
#     obs_vec_filepath = "experiments/calibration/land_observation_vector_SMAP_SM.jld2",
#     # keep other fields (date ranges, nelements, etc.)
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
    # This is a *vector* you’ve prebuilt from your RS LE (e.g., FLUXCOM/GLEAM),
    # using the same grid/mask/time aggregation that observation_map.jl will extract
    # from the model diagnostics for each minibatch.
    observation_vector = JLD2.load_object(CALIBRATE_CONFIG.obs_vec_filepath)

    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    minibatch_size     = CALIBRATE_CONFIG.minibatch_size

    # Note: names are purely cosmetic labels for EKP’s ObservationSeries
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
