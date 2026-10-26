# Calibrate the land model of the global calibration (models/snowy_land.jl),
# forced by ERA5, on columns representing global land, with the calibrated
# variables, noise, and prior of configs/energy_fluxes.jl (ERA5 energy fluxes,
# FLUXCOM GPP, and the MODIS snow cover fraction) and TransformUnscented EKI.
#
# Usage, from the root of the repository (`TopoSnowCoverFractionModel` reads
# topographic_data.nc in the working directory):
#   julia --project=.buildkite experiments/calibration/era5_columns/run_calibration.jl [OUTPUT_DIR]
#
# The forward runs (the ensemble members, and the runs with the default and the
# calibrated parameters) run one after another in this process or, with
# `CALIBRATION_BACKEND=caltech`, each as a Slurm job with one GPU on the Caltech
# cluster. With `ERA5_LOWRES_FORCING` set, the columns are forced by the 8° ERA5
# of 2008 instead of the hourly 1° ERA5. The calibrated parameters
# (`era5_parameters.toml`), figures, and a summary are written to
# `OUTPUT_DIR/results`. Rerunning with the same `OUTPUT_DIR` resumes an
# interrupted calibration (Caltech backend only).

import ClimaCalibrate
import EnsembleKalmanProcesses as EKP
import JLD2
import Random

const MODEL_INTERFACE_FILE = joinpath(@__DIR__, "model_interface.jl")
include(MODEL_INTERFACE_FILE)
include(joinpath(@__DIR__, "analysis.jl"))

"""Slurm directives of a forward run on the Caltech cluster, where GPUs are
requested by type with `--gres`."""
const CALTECH_DIRECTIVES = [
    :partition => "gpu",
    :reservation => "clima",
    :gres => "gpu:p100:1",
    :ntasks => 1,
    :cpus_per_task => 4,
    :mem => "32G",
    :time => 180,
]

"""
    caltech_backend()

Return a `CaltechHPCBackend` that runs each forward run as a Slurm job with one
GPU, with the `climacommon` module and the Julia depot of this process.
"""
function caltech_backend()
    modules = filter(
        m -> startswith(m, "climacommon/"),
        split(get(ENV, "LOADEDMODULES", ""), ':'),
    )
    env_vars =
        ["CLIMACOMMS_DEVICE" => "CUDA", "CLIMACOMMS_CONTEXT" => "SINGLETON"]
    for name in ("JULIA_DEPOT_PATH", "JULIA_CUDA_MEMORY_POOL")
        haskey(ENV, name) && push!(env_vars, name => ENV[name])
    end
    backend = ClimaCalibrate.CaltechHPCBackend(;
        directives = CALTECH_DIRECTIVES,
        modules = isempty(modules) ? String[] : [String(last(modules))],
        env_vars,
    )
    # `--gpus-per-task`, which `SlurmConfig` always sets, conflicts with `--gres`
    delete!(backend.hpc_config.directives, :gpus_per_task)
    return backend
end

"""
    run_evaluation(backend, interface, parameter_files, path)

Run the columns of `interface` with the default parameters overridden by
`parameter_files`, saving the seasonal means to `path`: in this process, or as
a job of an `HPCBackend`.
"""
run_evaluation(
    ::ClimaCalibrate.JuliaBackend,
    interface,
    parameter_files,
    path,
) = run_columns(interface, parameter_files, path)

function run_evaluation(
    backend::ClimaCalibrate.HPCBackend,
    interface,
    parameter_files,
    path,
)
    mkpath(dirname(path))
    julia_code = """
    include("$MODEL_INTERFACE_FILE")
    run_columns($(repr(interface)), $(repr(parameter_files)), "$path")
    """
    job_script = ClimaCalibrate.make_job_script(
        backend,
        "julia --project=$(ClimaCalibrate.project_dir()) -e '$julia_code'";
        job_name = "evaluation_$(basename(dirname(path)))",
        output = joinpath(dirname(path), "model_log.txt"),
    )
    job = ClimaCalibrate.submit_job(backend, job_script)
    while !ClimaCalibrate.iscompleted(job)
        sleep(30)
    end
    return nothing
end

"""
    evaluate(backend, interface, name, parameter_files)

Return the seasonal means simulated with the default parameters overridden by
`parameter_files`, running the columns unless an evaluation called `name` is
already saved.
"""
function evaluate(backend, interface, name, parameter_files)
    path =
        joinpath(interface.output_dir, "evaluation", name, SEASONAL_MEANS_FILE)
    if !isfile(path)
        @info "Running the $name evaluation"
        run_evaluation(backend, interface, parameter_files, path)
        isfile(path) || error(
            "The $name evaluation failed, see the log in $(dirname(path))",
        )
    end
    return JLD2.load_object(path)
end

"""
    final_ekp(output_dir)

Return the `EnsembleKalmanProcess` saved after the last ensemble update in
`output_dir` (iterations are numbered from 1).
"""
function final_ekp(output_dir)
    iteration = 1
    while isfile(ClimaCalibrate.ekp_path(output_dir, iteration + 1))
        iteration += 1
    end
    return ClimaCalibrate.load_ekp_struct(output_dir, iteration)
end

"""
    main(output_dir; backend = ClimaCalibrate.JuliaBackend(), lowres_forcing = false)

Run the calibration in `output_dir` and write its results to
`output_dir/results`.
"""
function main(
    output_dir;
    backend = ClimaCalibrate.JuliaBackend(),
    lowres_forcing = false,
)
    interface = ERA5ColumnsInterface(output_dir, lowres_forcing)
    results_dir = joinpath(output_dir, "results")
    mkpath(results_dir)
    isfile(columns_path(interface)) ||
        write_columns(columns_path(interface), select_columns())
    columns = read_columns(columns_path(interface))

    # Calibrate the observed entries at columns where the default model is
    # finite
    default = evaluate(backend, interface, "default", String[])
    observed = observed_seasonal_means(columns.longlat)
    masks = Dict(
        name => vec(isfinite.(observed[name]) .& isfinite.(default[name]))
        for name in SHORT_NAMES
    )
    for name in SHORT_NAMES
        n_unobserved = count(!isfinite, observed[name])
        n_nonfinite = count(!isfinite, default[name])
        @info "$name: calibrating $(count(masks[name])) seasonal means" n_unobserved n_nonfinite
    end
    JLD2.jldsave(joinpath(output_dir, OBSERVATIONS_FILE); observed, masks)

    prior = get_calibration_prior()
    ekp = EKP.EnsembleKalmanProcess(
        make_observation(observed, masks, columns.weight),
        EKP.TransformUnscented(prior; impose_prior = true);
        rng = Random.MersenneTwister(RNG_SEED),
        verbose = true,
        scheduler = EKP.DataMisfitController(terminate_at = 100),
        failure_handler_method = EKP.SampleSuccGauss(),
    )
    ClimaCalibrate.calibrate(
        backend,
        ekp,
        interface,
        N_ITERATIONS,
        prior,
        output_dir,
    )
    ekp = final_ekp(output_dir)

    parameter_file = joinpath(results_dir, "era5_parameters.toml")
    write_parameter_file(parameter_file, interface, prior, ekp)
    calibrated = evaluate(backend, interface, "calibrated", [parameter_file])

    metrics_default = error_metrics(default, observed, masks, columns.weight)
    metrics_calibrated =
        error_metrics(calibrated, observed, masks, columns.weight)
    JLD2.jldsave(
        joinpath(results_dir, "metrics.jld2");
        metrics_default,
        metrics_calibrated,
    )
    plot_parameters(joinpath(results_dir, "parameters.png"), ekp, prior)
    plot_rmse(
        joinpath(results_dir, "rmse.png"),
        metrics_default,
        metrics_calibrated,
    )
    plot_columns(
        joinpath(results_dir, "columns.png"),
        columns.longlat,
        columns.weight,
        column_loss_change(default, calibrated, observed, masks),
    )
    write_summary(
        joinpath(results_dir, "summary.md"),
        interface,
        prior,
        parameter_file,
        metrics_default,
        metrics_calibrated,
        relpath(results_dir),
    )
    @info "Calibration results written to $results_dir"
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    lowres_forcing = haskey(ENV, "ERA5_LOWRES_FORCING")
    # Fail early if the forcing or the data that are not artifacts are missing
    lowres_forcing ||
        ClimaLand.Artifacts.find_era5_year_paths(SPINUP_START, STOP_DATE)
    isfile(SCF_FILE) ||
        error("No snow cover fraction at $SCF_FILE: set ILAMB_SCF_FILE")
    isfile("topographic_data.nc") ||
        error("No topographic_data.nc: run from the root of the repository")
    backend_name = get(ENV, "CALIBRATION_BACKEND", "julia")
    backend =
        backend_name == "caltech" ? caltech_backend() :
        backend_name == "julia" ? ClimaCalibrate.JuliaBackend() :
        error("Unknown CALIBRATION_BACKEND $backend_name: use julia or caltech")
    main(abspath(get(ARGS, 1, "era5_calibration")); backend, lowres_forcing)
end
