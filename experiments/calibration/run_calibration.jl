# Calibration of the global land model

# The code sets up and runs a calibration of the global land model or bucket
# model depending on the `model_type` in the calibration config. To start a
# simulation on Derecho or GCP, you run
# `bash experiments/calibration/run_calibration.sh` in the root directory.

using Dates
using Distributed
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaLand
import EnsembleKalmanProcesses as EKP
import JLD2

include(joinpath(pkgdir(ClimaLand), "experiments", "calibration", "api.jl"))
include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments",
        "calibration",
        "observation_map.jl",
    ),
)

model_interface = joinpath(
    pkgdir(ClimaLand),
    "experiments",
    "calibration",
    "model_interface.jl",
)

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["lwu"],
    minibatch_size = 1,
    n_iterations = 1,
    sample_date_ranges = [("2007-12-1", "2007-12-1")],
    extend = Dates.Month(3),
    spinup = Dates.Month(0),
    nelements = (180, 360, 15),
    output_dir = "experiments/calibration/land_model",
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector.jld2",
    model_type = ClimaLand.LandModel,
)


if abspath(PROGRAM_FILE) == @__FILE__
    # true solution is at 0.96
    priors =
        [EKP.constrained_gaussian("emissivity_bare_soil", 0.82, 0.12, 0.0, 2.0)]
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
        EKP.TransformUnscented(prior, impose_prior = true);
        verbose = true,
        rng,
        scheduler = EKP.DataMisfitController(terminate_at = 100),
    )

    # Note: Using this script on Derecho requires changes to addprocs to use
    # the PBSManager
    curr_backend = ClimaCalibrate.get_backend()
    N_ens = EKP.get_N_ens(ekp)

    if curr_backend == ClimaCalibrate.DerechoBackend
        # The HPC keyword arguments specify the resources for running a single
        # ensemble member
        backend = ClimaCalibrate.DerechoBackend(;
            model_interface,
            verbose = true,
            hpc_kwargs = Dict(
                :cpus_per_task => 4,
                :gpus_per_task => 1,
                # Options include "premium", "regular", "economy", "preempt"
                :job_priority => "regular",
                :ntasks => 1,
                # 180 minutes is 3 hours
                :time => 180,
            ),
        )
    elseif curr_backend == ClimaCalibrate.GCPBackend
        backend = ClimaCalibrate.GCPBackend(;
            model_interface,
            verbose = true,
            hpc_kwargs = Dict(
                :cpus_per_task => 4,
                :gpus_per_task => 1,
                :ntasks => 1,
                :partition => "a3mega",
                :time => 180,
            ),
        )
    elseif curr_backend == ClimaCalibrate.ClimaGPUBackend
        # This backend is used for testing. Since it is not easy to specify only
        # one worker with ClimaGPUBackend, we use addprocs and the WorkerBackend
        # instead.
        @info "Check your slurm script that the number of tasks is the same as the ensemble size ($N_ens)"
        addprocs(ClimaCalibrate.SlurmManager())
        @everywhere import ClimaLand
        @everywhere experiment_dir = joinpath(pkgdir(ClimaLand), "experiments")
        @everywhere include(
            joinpath(pkgdir(ClimaLand), "experiments", "calibration", "api.jl"),
        )
        @everywhere const CALIBRATE_CONFIG = $CALIBRATE_CONFIG
        @everywhere include(
            joinpath(experiment_dir, "calibration", "model_interface.jl"),
        )
        (; n_iterations, output_dir) = CALIBRATE_CONFIG
        eki = ClimaCalibrate.calibrate(
            ClimaCalibrate.WorkerBackend(),
            ekp,
            n_iterations,
            prior,
            output_dir;
        )
        return nothing
    else
        error("Unsupported backend: $(ClimaCalibrate.get_backend())")
    end

    if curr_backend != ClimaCalibrate.ClimaGPUBackend
        (; n_iterations, output_dir) = CALIBRATE_CONFIG
        eki = ClimaCalibrate.calibrate(
            backend,
            ekp,
            n_iterations,
            prior,
            output_dir;
        )
    end
end
