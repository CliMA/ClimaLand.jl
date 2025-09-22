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
    short_names = ["lhf"],# later add shf, lwu
    minibatch_size = 1,
    n_iterations = 10,
    sample_date_ranges = [("$(2000 + 2*i)-12-1", "$(2002 + 2*i)-9-1") for i in 0:9], # 2000 to 2020
    extend = Dates.Month(3),
    spinup = Dates.Month(3),
    nelements = (101, 15),
    output_dir = "/glade/derecho/scratch/kdeck/p-model-cal",
    rng_seed = 42,
)


if abspath(PROGRAM_FILE) == @__FILE__
    # Note: Using this script on Derecho requires changes to addprocs to use
    # the PBSManager
    if ClimaCalibrate.get_backend() == ClimaCalibrate.DerechoBackend
        addprocs(
            ClimaCalibrate.PBSManager(21),
            q = "main",
            A = "UCIT0011",
            l_select = "1:ngpus=1:ncpus=4",
            l_walltime = "11:30:00",
        )
    elseif ClimaCalibrate.get_backend() == ClimaCalibrate.CaltechHPCBackend
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

    priors = [
        EKP.constrained_gaussian("low_water_pressure_sensitivity", 5e-6, 4e-6, 0, Inf),
        EKP.constrained_gaussian("moisture_stress_ref_water_pressure", -2e6, 1e6, -Inf, 0),
        EKP.constrained_gaussian("a", 0.00196, 0.0007, 0, Inf),
        EKP.constrained_gaussian("K_sat_plant", 7e-8, 3e-8, 0, 1e-6),
        EKP.constrained_gaussian("psi_63", -408.16, 100.0, -1000.0, -50.0),
        EKP.constrained_gaussian("Weibull_c", 4, 1, 0.2, 6),
        EKP.constrained_gaussian("pmodel_cstar", 0.41, 0.11, 0, Inf),
        EKP.constrained_gaussian("pmodel_β", 146, 10, 0, Inf),
        EKP.constrained_gaussian("pmodel_ϕ0_c3", 0.052, 0.02, 0, Inf),
        EKP.constrained_gaussian("pmodel_ϕ0_c4", 0.057, 0.02, 0, Inf),
#       EKP.constrained_gaussian("emissivity_bare_soil", 0.96, 0.03, 0.0, 1.0),
#       EKP.constrained_gaussian("canopy_emissivity", 0.96, 0.03, 0.0, 1.0),
#       EKP.constrained_gaussian("ac_canopy", 7.5e3, 5e3, 1e3, Inf)
    ]
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
