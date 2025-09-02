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
    short_names = ["swu"],
    minibatch_size = 1,
    n_iterations = 7,
    sample_date_ranges = [("2008-12-1", "2009-9-1"), ("2009-12-1", "2010-9-1"), ("2010-12-1", "2011-9-1"), ("2011-12-1", "2012-9-1"), ("2012-12-1", "2013-9-1"), ("2013-12-1", "2014-9-1")],
    extend = Dates.Month(3),
    spinup = Dates.Year(1),
    nelements = (101, 15),
    output_dir = "experiments/calibration/land_model",
    rng_seed = 42,
)


if abspath(PROGRAM_FILE) == @__FILE__
    # Note: Using this script on Derecho requires changes to addprocs to use
    # the PBSManager
    addprocs(ClimaCalibrate.SlurmManager())

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

    # true solution is at 0.96
    priors =
        [EKP.constrained_gaussian("alpha_0", 0.6, 0.1, 0.1, 0.8),
        EKP.constrained_gaussian("delta_alpha", 0.2, 0.05, 0.0, 0.3),
        EKP.constrained_gaussian("k", 10.0, 2.0, 1.0, 20.0),
        EKP.constrained_gaussian("beta", 0.4, 0.2, 0.05, 0.9),
        EKP.constrained_gaussian("x0", 0.2, 0.05, 0.1, 0.8)]
    prior = EKP.combine_distributions(priors)

    observation_vector =
        JLD2.load_object("experiments/calibration/land_observation_vector.jld2")

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

    eki = ClimaCalibrate.calibrate(
        ClimaCalibrate.WorkerBackend,
        ekp,
        CALIBRATE_CONFIG.n_iterations,
        prior,
        CALIBRATE_CONFIG.output_dir,
    )
end
