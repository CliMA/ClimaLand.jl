using ClimaLand
dir = pkgdir(ClimaLand)
import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL
using Distributed
addprocs(CAL.SlurmManager())

@everywhere begin
    import ClimaCalibrate:
        forward_model, parameter_path, path_to_ensemble_member
    import ClimaCalibrate as CAL

    using ClimaLand
    caldir = "calibration_output"
    dir = pkgdir(ClimaLand)

    # include(joinpath(dir,"experiments/calibration/global_bucket/climacalibrate_bucket/bucket_target_script.jl"))
    include(
        joinpath(
            dir,
            "experiments/calibration/global_land/forward_model_land.jl",
        ),
    )
end

include(
    joinpath(dir, "experiments/calibration/observations/observation_map.jl"),
)
include(
    joinpath(dir, "experiments/calibration/observations/observation_data.jl"),
)


## Priors
prior_pc = EKP.constrained_gaussian("pc", -2.5e6, 1e6, -Inf, 0);
prior_sc = EKP.constrained_gaussian("sc", -2e-6, 1e-6, 0, 1);
prior = EKP.combine_distributions([prior_pc, prior_sc]);

ensemble_size = 2
n_iterations = 2
noise = 1.0 * EKP.I # Should work, but this should be covariance of each month from observation (ERA5)

caldir = "calibration_output"

CAL.calibrate(
    CAL.WorkerBackend,
    ensemble_size,
    n_iterations,
    observations,
    noise,
    prior,
    caldir,
)
