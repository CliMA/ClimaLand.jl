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
model_interface_filepath = joinpath(
    pkgdir(ClimaLand),
    "experiments",
    "calibration",
    "model_interface.jl",
)
include(model_interface_filepath)

# Note: This has only been tested with the WorkerBackend
const TEST_CALIBRATION = haskey(ENV, "TEST_CALIBRATION")

# Optional CLI override for the calibration output directory. Pass the target
# path as the first positional argument to run_calibration.jl /
# generate_observations.jl (e.g. bash run_calibration.sh /scratch/foo).
#
# Per-member worker jobs run a generated `model_run.jl` directly, so ARGS is
# empty there. Recover the output_dir from PROGRAM_FILE in that case — the
# worker script always lives at $output_dir/iteration_NNN/member_MMM/model_run.jl.
const OUTPUT_DIR = if length(ARGS) >= 1
    ARGS[1]
elseif basename(PROGRAM_FILE) == "model_run.jl"
    dirname(dirname(dirname(abspath(PROGRAM_FILE))))
else
    "experiments/calibration/land_model"
end

# Include the calibration configuration. This defines CALIBRATE_CONFIG,
# get_calibration_prior(), and NOISE_SCALARS. When TEST_CALIBRATION is set,
# load the single-parameter test config; otherwise load the production config.
# To run a different production calibration, set the CALIBRATION_CONFIG env
# var before invoking run_calibration.sh, e.g.
#   CALIBRATION_CONFIG=gpp.jl bash experiments/calibration/run_calibration.sh
const CONFIG_FILE =
    TEST_CALIBRATION ? "test.jl" :
    get(ENV, "CALIBRATION_CONFIG", "energy_fluxes.jl")
include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments",
        "calibration",
        "configs",
        CONFIG_FILE,
    ),
)

# Workaround for ClimaCalibrate v0.3.0: its job script inlines a multi-line
# program into `julia -e`, which Derecho's set_gpu_rank wrapper word-splits on
# newlines so julia only parses `import`. Restore v0.2.2's behavior of writing
# model_run.jl to disk and passing the path. Remove once fixed upstream.
function ClimaCalibrate.Calibration.generate_job_script_for_ensemble_member(
    backend::ClimaCalibrate.Backend.HPCBackend,
    iter,
    member,
    output_dir,
    model_interface_filepath,
    experiment_dir,
    exeflags,
)
    job_body = """
    import ClimaCalibrate
    iteration = $iter; member = $member
    model_interface_filepath = "$model_interface_filepath"
    include(model_interface_filepath)
    interface = ClimaCalibrate._load(joinpath("$output_dir", "interface.jld2"))
    ClimaCalibrate.forward_model(interface, iteration, member)
    ClimaCalibrate.write_model_completed("$output_dir", iteration, member)
    """
    member_path =
        ClimaCalibrate.Calibration.path_to_ensemble_member(output_dir, iter, member)
    mkpath(member_path)
    julia_filepath = joinpath(member_path, "model_run.jl")
    write(julia_filepath, job_body)

    julia_command = "julia --project=$experiment_dir $exeflags $julia_filepath"

    member_log = ClimaCalibrate.Calibration.path_to_model_log(output_dir, iter, member)
    return ClimaCalibrate.Backend.make_job_script(
        backend,
        julia_command;
        job_name = "run_$(iter)_$(member)",
        output = member_log,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    prior = get_calibration_prior()

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

    (; rng_seed) = CALIBRATE_CONFIG
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

    model_interface = LandModelInterface(CALIBRATE_CONFIG)

    # Run test calibration (only support the WorkerBackend and clusters with
    # Slurm)
    if TEST_CALIBRATION
        # This backend is used for testing. Since it is not easy to specify only
        # one worker with ClimaGPUBackend, we use addprocs and the WorkerBackend
        # instead.
        @info "Check your slurm script that the number of tasks is the same as the ensemble size ($N_ens)"
        addprocs(ClimaCalibrate.SlurmManager())
        @everywhere include($model_interface_filepath)
        (; n_iterations, output_dir) = CALIBRATE_CONFIG
        eki = ClimaCalibrate.calibrate(
            ClimaCalibrate.WorkerBackend(),
            ekp,
            model_interface,
            n_iterations,
            prior,
            output_dir,
        )
        return nothing
    end

    # The directives specify the resources for running a single
    # ensemble member
    directives = [
        :cpus_per_task => 4,
        :gpus_per_task => 1,
        :ntasks => 1,
        # 180 minutes is 3 hours
        :time => 180,
    ]
    # Pin the module version so forward-model jobs match the project Manifests;
    # unversioned `climacommon` resolves to whatever is current.
    modules = ["climacommon/2025_02_25"]
    # HDF5_USE_FILE_LOCKING=FALSE is required on Derecho's Lustre scratch.
    env_vars = [
        "CLIMACOMMS_CONTEXT" => "SINGLETON",
        "CLIMACOMMS_DEVICE" => "CUDA",
        "HDF5_USE_FILE_LOCKING" => "FALSE",
    ]

    # Determine which backend to submit job scripts to
    if curr_backend == ClimaCalibrate.DerechoBackend
        derecho_directives = [:job_priority => "regular"]
        backend = ClimaCalibrate.DerechoBackend(;
            directives = vcat(directives, derecho_directives),
            modules,
            env_vars,
        )
    elseif curr_backend == ClimaCalibrate.GCPBackend
        gcp_directives = [:partition => "a3mega"]
        backend = ClimaCalibrate.GCPBackend(;
            directives = vcat(directives, gcp_directives),
            modules,
            env_vars,
        )
    elseif curr_backend == ClimaCalibrate.CaltechHPCBackend
        central_directives = [:partition => "gpu"]
        backend = ClimaCalibrate.CaltechHPCBackend(;
            directives = vcat(directives, central_directives),
            modules,
            env_vars,
        )
    elseif curr_backend == ClimaCalibrate.ClimaGPUBackend
        backend =
            ClimaCalibrate.ClimaGPUBackend(; directives, modules, env_vars)
    else
        error("Unsupported backend: $(ClimaCalibrate.get_backend())")
    end

    (; n_iterations, output_dir) = CALIBRATE_CONFIG
    eki = ClimaCalibrate.calibrate(
        backend,
        ekp,
        model_interface,
        n_iterations,
        prior,
        output_dir,
    )
end
