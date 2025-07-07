using Dates
using Distributed
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaComms
import ClimaLand
import EnsembleKalmanProcesses as EKP
import JLD2

addprocs(ClimaCalibrate.SlurmManager())

include(
    joinpath(pkgdir(ClimaLand), "experiments/better_calibration/getters.jl"),
)

(; output_dir) = get_config()

include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments/better_calibration/observation_map.jl",
    ),
)

@everywhere begin
    import ClimaLand
    experiment_dir = joinpath(pkgdir(ClimaLand), "experiments")
    include(
        joinpath(experiment_dir, "better_calibration", "model_interface.jl"),
    )
end

prior = get_prior()

# Note: You should check that the ensemble size is the same as the number of
# tasks in the batch script
# For example, if you calibrating 3 parameters and are using
# EKP.TransformUnscented, then the number of tasks should be 7, since
# 3 * 2 + 1 = 7
(; n_iterations) = get_calibration_config();
ekp = get_ekp();
eki = ClimaCalibrate.calibrate(
    ClimaCalibrate.WorkerBackend,
    ekp,
    n_iterations,
    prior,
    output_dir,
)
