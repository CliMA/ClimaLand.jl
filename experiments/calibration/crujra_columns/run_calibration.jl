# Calibrate the land model of the default long run (P-model, prescribed MODIS
# LAI), forced by CRUJRA, on the columns of `columns.csv` against seasonal means
# of FLUXCOM LHF, SHF, GPP and CERES SWU, LWU, with TransformUnscented EKI.
#
# Usage:
#   julia --project=.buildkite experiments/calibration/crujra_columns/run_calibration.jl [OUTPUT_DIR]
#
# Ensemble members run one after another in this process, or concurrently on
# `CALIBRATION_N_WORKERS` local worker processes that share the device chosen by
# `CLIMACOMMS_DEVICE`. The calibrated parameters (`crujra_parameters.toml`),
# figures, and a summary are written to `OUTPUT_DIR/results`. Rerunning with the
# same `OUTPUT_DIR` resumes an interrupted calibration (workers only).

using Distributed
import ClimaCalibrate
import EnsembleKalmanProcesses as EKP
import JLD2
import Random

const MODEL_INTERFACE_FILE = joinpath(@__DIR__, "model_interface.jl")
include(MODEL_INTERFACE_FILE)
include(joinpath(@__DIR__, "analysis.jl"))

"""
    evaluate(name, override_files, output_dir)

Return the seasonal means simulated with the default parameters overridden by
`override_files`, running the columns (on a worker, if any) unless an evaluation
called `name` is already saved in `output_dir`.
"""
function evaluate(name, override_files, output_dir)
    path = joinpath(output_dir, "evaluation", name, SEASONAL_MEANS_FILE)
    if !isfile(path)
        @info "Running the $name evaluation"
        if workers() == [myid()]
            run_columns(override_files, path)
        else
            remotecall_fetch(
                run_columns,
                first(workers()),
                override_files,
                path,
            )
        end
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
    main(output_dir; n_workers = 0)

Run the calibration in `output_dir` and write its results to
`output_dir/results`.
"""
function main(output_dir; n_workers = 0)
    results_dir = joinpath(output_dir, "results")
    mkpath(results_dir)
    columns = read_columns()

    if n_workers > 0
        addprocs(n_workers; exeflags = "--project=$(Base.active_project())")
        @everywhere workers() include($MODEL_INTERFACE_FILE)
        backend = ClimaCalibrate.WorkerBackend()
    else
        backend = ClimaCalibrate.JuliaBackend()
    end

    # Calibrate the observed entries at columns where the default model is
    # finite
    default = evaluate("default", String[], output_dir)
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
        CRUJRAColumnsInterface(output_dir),
        N_ITERATIONS,
        prior,
        output_dir,
    )
    ekp = final_ekp(output_dir)

    parameter_file = joinpath(results_dir, "crujra_parameters.toml")
    write_parameter_file(parameter_file, prior, ekp)
    calibrated = evaluate("calibrated", [parameter_file], output_dir)

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
    # Fail early if the forcing is not available on this machine
    ClimaLand.Artifacts.find_crujra_year_paths(SPINUP_START, STOP_DATE)
    main(
        abspath(get(ARGS, 1, "crujra_calibration"));
        n_workers = parse(Int, get(ENV, "CALIBRATION_N_WORKERS", "0")),
    )
end
