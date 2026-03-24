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
    short_names = ["gpp"],
    minibatch_size = 4,
    n_iterations = 10,
    # 10 yearly samples: each covers DJF through SON (Dec 1 → Sep 1)
    # with extend = Month(3), simulation runs through Nov 30
    sample_date_ranges = [
        ("$(year)-12-1", "$(year+1)-9-1") for year in 2000:2009
    ],
    extend = Dates.Month(3),
    spinup = Dates.Year(1),
    nelements = (180, 360, 15),
    output_dir = "experiments/calibration/land_model",
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector.jld2",
    model_type = ClimaLand.LandModel,
)


if abspath(PROGRAM_FILE) == @__FILE__
    # 4 P-model parameters + 1 soil moisture stress parameter
    # Ensemble size for TransformUnscented: 5 * 2 + 1 = 11 members
    # Note: ϕ0_c3/ϕ0_c4 are not calibrated because temperature_dep_yield = true
    # uses the quadratic coefficients (ϕa0, ϕa1, ϕa2) instead.
    priors = [
        EKP.constrained_gaussian("pmodel_cstar", 0.41, 0.05, 0.2, 0.7),
        EKP.constrained_gaussian("pmodel_β_c3", 146.0, 40.0, 50.0, 300.0),
        EKP.constrained_gaussian("pmodel_β_c4", 16.222, 5.0, 5.0, 40.0),
        EKP.constrained_gaussian("pmodel_α", 0.933, 0.02, 0.85, 0.999),
        EKP.constrained_gaussian("moisture_stress_c", 0.27, 0.15, 0.05, 1.0),
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
