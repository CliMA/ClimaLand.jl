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

if !TEST_CALIBRATION
    const CALIBRATE_CONFIG = CalibrateConfig(;
        short_names = ["gpp"],
        minibatch_size = 4,
        n_iterations = 10,
        # 10 yearly samples: each covers DJF, MAM, JJA, SON (Dec 1 → Dec 1)
        # with extend = Month(3), simulation runs through Feb 28 of year+2
        sample_date_ranges = [
            ("$(year)-12-1", "$(year+1)-12-1") for year in 2001:2010
        ],
        extend = Dates.Month(3),
        spinup = Dates.Year(1),
        nelements = (180, 360, 15),
        output_dir = "/glade/derecho/scratch/arenchon/calibration_gpp",
        rng_seed = 42,
        obs_vec_filepath = "experiments/calibration/land_observation_vector.jld2",
        model_type = ClimaLand.LandModel,
    )
else
    @info "Using calibration config for test calibration"
    const CALIBRATE_CONFIG = CalibrateConfig(;
        short_names = ["gpp"],
        minibatch_size = 1,
        n_iterations = 1,
        sample_date_ranges = [("2007-12-1", "2007-12-1")],
        extend = Dates.Month(3),
        spinup = Dates.Month(0),
        nelements = (180, 360, 15),
        output_dir = "/glade/derecho/scratch/arenchon/calibration_gpp",
        rng_seed = 42,
        obs_vec_filepath = "experiments/calibration/land_observation_vector.jld2",
        model_type = ClimaLand.LandModel,
    )
end

"""
    module_load_string(::ClimaCalibrate.ClimaGPUBackend)

Load the appropriate module for `clima`.

This is needed to load the right `climacommon` version on each cluster. This
function will be removed after support is added in ClimaCalibrate.jl for
choosing which version of `climacommon` to load.
"""
function ClimaCalibrate.module_load_string(::ClimaCalibrate.ClimaGPUBackend)
    return """module purge
    module load climacommon/2026_02_18"""
end

function ClimaCalibrate.module_load_string(::ClimaCalibrate.DerechoBackend)
    return """module purge
    module load climacommon/2025_02_25"""
end

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
