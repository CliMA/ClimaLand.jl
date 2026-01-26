using Dates
using Distributed
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaComms
import ClimaLand
import EnsembleKalmanProcesses as EKP
import JLD2

include(joinpath(pkgdir(ClimaLand), "experiments/calibration/api.jl"))

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["lwu", "shf", "lhf"],
    minibatch_size = 1,
    n_iterations = 12,
    sample_date_ranges = [
        ("$(2000 + 1*i)-12-1", "$(2001 + 1*i)-9-1") for i in 0:12
    ], # 2000 to 2020
    extend = Dates.Month(3),
    spinup = Dates.Month(3),
    nelements = (180, 360, 15),
    output_dir = "/glade/derecho/scratch/kdeck/recalibrate_saturated_K_lw_year",
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector_year.jld2",
    model_type = ClimaLand.LandModel,
)


if abspath(PROGRAM_FILE) == @__FILE__
    priors = [
        EKP.constrained_gaussian("moisture_stress_c", 0.5, 0.1, 0, 2),
        EKP.constrained_gaussian("pmodel_cstar", 0.41, 0.11, 0, Inf),
        EKP.constrained_gaussian("pmodel_β", 146, 10, 0, Inf),
        EKP.constrained_gaussian("leaf_Cd", 0.01, 0.005, 0, Inf),
        EKP.constrained_gaussian("canopy_K_lw", 0.5, 0.25, 0, Inf),
    ]
    prior = EKP.combine_distributions(priors)

    observation_vector = JLD2.load_object(CALIBRATE_CONFIG.obs_vec_filepath)

    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    minibatch_size = CALIBRATE_CONFIG.minibatch_size
    obs_series = EKP.ObservationSeries(
        Dict(
            "observations" => observation_vector,
            "names" => [
                string(Dates.year(start_date)) for
                (start_date, stop_date) in sample_date_ranges
            ],
            "minibatcher" => ClimaCalibrate.minibatcher_over_samples(
                length(observation_vector),
                minibatch_size,
            ),
        ),
    )

    rng_seed = CALIBRATE_CONFIG.rng_seed
    rng = Random.MersenneTwister(rng_seed)

    # Note: You should check that the ensemble size is the same as the number of
    # tasks in the batch script
    # For example, if you are calibrating 3 parameters and are using
    # EKP.TransformUnscented, then the number of tasks should be 7, since
    # 3 * 2 + 1 = 7
    ekp = EKP.EnsembleKalmanProcess(
        obs_series,
        EKP.TransformUnscented(prior, impose_prior = true),
        verbose = true,
        rng = rng,
        scheduler = EKP.DataMisfitController(terminate_at = 100),
    )

    # Note: Using this script on Derecho requires changes to addprocs to use
    # the PBSManager
    curr_backend = ClimaCalibrate.get_backend()
    N_ens = EKP.get_N_ens(ekp)
    if curr_backend == ClimaCalibrate.DerechoBackend
        addprocs(
            ClimaCalibrate.PBSManager(N_ens),
            q = "main",
            A = "UCIT0011",
            l_select = "1:ngpus=1:ncpus=4",
            l_walltime = "11:30:00",
        )
    elseif (curr_backend == ClimaCalibrate.CaltechHPCBackend) ||
           (curr_backend == ClimaCalibrate.GCPBackend)
        @info "Check your slurm script that the number of tasks is the same as the ensemble size ($N_ens)"
        addprocs(ClimaCalibrate.SlurmManager())
    end

    include(
        joinpath(
            pkgdir(ClimaLand),
            "experiments/calibration/observation_map.jl",
        ),
    )

    @everywhere import ClimaLand
    @everywhere experiment_dir = joinpath(pkgdir(ClimaLand), "experiments")
    @everywhere include(
        joinpath(pkgdir(ClimaLand), "experiments/calibration/api.jl"),
    )
    @everywhere CALIBRATE_CONFIG = $CALIBRATE_CONFIG
    @everywhere include(
        joinpath(experiment_dir, "calibration", "model_interface.jl"),
    )

    eki = ClimaCalibrate.calibrate(
        ClimaCalibrate.WorkerBackend,
        ekp,
        CALIBRATE_CONFIG.n_iterations,
        prior,
        CALIBRATE_CONFIG.output_dir,
    )
end
