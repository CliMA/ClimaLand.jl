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
    short_names = ["lwu"],
    minibatch_size = 1,
    n_iterations = 3,
    sample_date_ranges = [("2007-12-1", "2008-9-1")],
    extend = Dates.Month(3),
    spinup = Dates.Month(3),
    nelements = (101, 15),
    output_dir = "experiments/calibration/land_model",
    rng_seed = 42,
)

# true solution is at 0.96
priors =
    [EKP.constrained_gaussian("emissivity_bare_soil", 0.82, 0.12, 0.0, 2.0)]
prior = EKP.combine_distributions(priors)
ekp_process = EKP.Unscented(prior)
ensemble_size = ekp_process.N_ens

if abspath(PROGRAM_FILE) == @__FILE__
    # addprocs is slightly different on derecho or central. Could add more.
    if contains(gethostname(), "derecho") # Derecho
        addprocs(
            ClimaCalibrate.PBSManager(ensemble_size),
            q = "main",
            A = "UCIT0011",
            l_select = "1:ngpus=1:ncpus=1",
            l_walltime = "05:30:00",
        )
    elseif contains(gethostname(), "cluster") # Caltech central
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


    observation_vector =
        JLD2.load_object("experiments/calibration/land_observation_vector.jld2")

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
            "metadata" => [
                prior,
                CALIBRATE_CONFIG.short_names,
                CALIBRATE_CONFIG.nelements,
            ],
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
