using Dates
using Distributed
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaComms
import ClimaLand
import EnsembleKalmanProcesses as EKP
import JLD2

addprocs(ClimaCalibrate.SlurmManager())

# TODO: Remove this after api/config is okay
# include(
#     joinpath(pkgdir(ClimaLand), "experiments/better_calibration/getters.jl"),
# )

include(
    joinpath(pkgdir(ClimaLand), "experiments/better_calibration/config.jl"),
)

include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments/better_calibration/observation_map.jl",
    ),
)

@everywhere begin
    import ClimaLand
    experiment_dir = joinpath(pkgdir(ClimaLand), "experiments")
    include(
        joinpath(experiment_dir, "better_calibration", "model_interface.jl"),
    )
    CALIBRATE_CONFIG
end

priors = [
        EKP.constrained_gaussian("wet_1", 0.375, 0.2, 0.0, 1.0),
        EKP.constrained_gaussian("dry_1", 0.485, 0.2, 0.0, 1.0),
        EKP.constrained_gaussian("wet_2", 0.345, 0.1, 0.0, 1.0),
        EKP.constrained_gaussian("dry_2", 0.455, 0.2, 0.0, 1.0),
        EKP.constrained_gaussian("wet_3", 0.315, 0.1, 0.0, 1.0),
        EKP.constrained_gaussian("dry_3", 0.425, 0.2, 0.0, 1.0),
        EKP.constrained_gaussian("wet_4", 0.3, 0.1, 0.0, 1.0),
        EKP.constrained_gaussian("dry_4", 0.41, 0.2, 0.0, 1.0),
    ]
prior = EKP.combine_distributions(priors)

observation_vector = JLD2.load_object(
        "experiments/better_calibration/land_observation_vector.jld2",
    )

sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
minibatch_size = CALIBRATE_CONFIG.minibatch_size
obs_series = EKP.ObservationSeries(
        Dict(
            "observations" => observation_vector,
            "names" => [
                string(Dates.year(start_date)) for
                (start_date, end_date) in sample_date_ranges
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
# For example, if you calibrating 3 parameters and are using
# EKP.TransformUnscented, then the number of tasks should be 7, since
# 3 * 2 + 1 = 7
ekp = EKP.EnsembleKalmanProcess(
        obs_series,
        EKP.TransformUnscented(prior, impose_prior = true),
        verbose = true,
        rng = rng,
        scheduler = EKP.DataMisfitController(terminate_at = 100),
    )

eki = ClimaCalibrate.calibrate(
    ClimaCalibrate.WorkerBackend,
    ekp,
    CALIBRATE_CONFIG.n_iterations,
    prior,
    CALIBRATE_CONFIG.output_dir,
)
