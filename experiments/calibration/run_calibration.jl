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

"""
    _loaded_climacommon()

Return the `climacommon/<version>` module name currently loaded in the parent
shell, by scanning the `LOADEDMODULES` environment variable (colon-separated).
Falls back to `"climacommon"` (default version) if none is found.

This lets forward-model jobs pick up whichever `climacommon` the driver script
loaded, instead of hardcoding a version per backend.
"""
function _loaded_climacommon()
    loaded = get(ENV, "LOADEDMODULES", "")
    # Expect entries of the form `climacommon/YYYY_MM_DD`. If multiple are
    # loaded, the last one wins (matches `module`'s own resolution order).
    pattern = r"^climacommon/\d{4}_\d{2}_\d{2}$"
    matches = filter(mod -> occursin(pattern, mod), split(loaded, ':'))
    isempty(matches) && return "climacommon"
    return String(last(matches))
end

"""
    _forwarded_env_exports()

Return shell-export lines for env vars that must be propagated from the driver
shell into each forward-model PBS/Slurm job. PBS does not inherit the submit
environment without `-V`, so worker jobs would otherwise re-read `ENV` and see
an empty `CALIBRATION_CONFIG`, silently falling back to the default config in
`run_calibration.jl` and producing sims that disagree with the obs vector.
"""
function _forwarded_env_exports()
    forwarded = ("CALIBRATION_CONFIG", "HDF5_USE_FILE_LOCKING")
    lines = String[]
    for name in forwarded
        haskey(ENV, name) || continue
        push!(lines, "export $name=$(ENV[name])")
    end
    return join(lines, "\n")
end

"""
    module_load_string(::ClimaCalibrate.ClimaGPUBackend)
    module_load_string(::ClimaCalibrate.DerechoBackend)

Return the `module load` string injected into forward-model job scripts. The
`climacommon` version is read from the driver shell's `LOADEDMODULES` so the
version only needs to be specified once (in `run_calibration.sh`).
"""
function ClimaCalibrate.module_load_string(::ClimaCalibrate.ClimaGPUBackend)
    return """module purge
    module load $(_loaded_climacommon())
    $(_forwarded_env_exports())"""
end

function ClimaCalibrate.module_load_string(::ClimaCalibrate.DerechoBackend)
    return """module purge
    module load $(_loaded_climacommon())
    $(_forwarded_env_exports())"""
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

    # Run test calibration (only support the WorkerBackend and clusters with
    # Slurm)
    if TEST_CALIBRATION
        # This backend is used for testing. Since it is not easy to specify only
        # one worker with ClimaGPUBackend, we use addprocs and the WorkerBackend
        # instead.
        @info "Check your slurm script that the number of tasks is the same as the ensemble size ($N_ens)"
        addprocs(ClimaCalibrate.SlurmManager())
        @everywhere include($model_interface)
        (; n_iterations, output_dir) = CALIBRATE_CONFIG
        eki = ClimaCalibrate.calibrate(
            ClimaCalibrate.WorkerBackend(),
            ekp,
            n_iterations,
            prior,
            output_dir,
        )
        return nothing
    end

    # The HPC keyword arguments specify the resources for running a single
    # ensemble member
    hpc_kwargs = Dict(
        :cpus_per_task => 4,
        :gpus_per_task => 1,
        :ntasks => 1,
        # 180 minutes is 3 hours
        :time => 180,
    )

    # Determine which backend to submit job scripts to
    if curr_backend == ClimaCalibrate.DerechoBackend
        derecho_hpc_kwargs = Dict(
            # Options include "premium", "regular", "economy", "preempt"
            :job_priority => "regular",
        )
        backend = ClimaCalibrate.DerechoBackend(;
            model_interface,
            verbose = true,
            hpc_kwargs = merge(hpc_kwargs, derecho_hpc_kwargs),
        )
    elseif curr_backend == ClimaCalibrate.GCPBackend
        gcp_hpc_kwargs = Dict(:partition => "a3mega")
        backend = ClimaCalibrate.GCPBackend(;
            model_interface,
            verbose = true,
            hpc_kwargs = merge(hpc_kwargs, gcp_hpc_kwargs),
        )
    elseif curr_backend == ClimaCalibrate.CaltechHPCBackend
        central_hpc_kwargs = Dict(:partition => "gpu")
        backend = ClimaCalibrate.CaltechHPCBackend(;
            model_interface,
            verbose = true,
            hpc_kwargs = merge(hpc_kwargs, central_hpc_kwargs),
        )
    elseif curr_backend == ClimaCalibrate.ClimaGPUBackend
        clima_hpc_kwargs = Dict()
        backend = ClimaCalibrate.ClimaGPUBackend(;
            model_interface,
            verbose = true,
            hpc_kwargs = merge(hpc_kwargs, clima_hpc_kwargs),
        )
    else
        error("Unsupported backend: $(ClimaCalibrate.get_backend())")
    end

    (; n_iterations, output_dir) = CALIBRATE_CONFIG
    eki =
        ClimaCalibrate.calibrate(backend, ekp, n_iterations, prior, output_dir)
end
