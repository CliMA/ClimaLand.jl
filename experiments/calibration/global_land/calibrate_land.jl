using ClimaLand
dir = pkgdir(ClimaLand)
import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL
using Distributed
# addprocs(CAL.SlurmManager())

# @everywhere begin
import ClimaCalibrate:
    forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL

using ClimaLand
caldir = "calibration_output"
dir = pkgdir(ClimaLand)

model_interface = joinpath(dir,"experiments/calibration/global_land/forward_model_land.jl")
include(model_interface)
# end

include(
    joinpath(dir, "experiments/calibration/global_land/observations/observation_map.jl"),
)
include(
    joinpath(dir, "experiments/calibration/global_land/observations/observation_data.jl"),
)


## Priors
prior_pc = EKP.constrained_gaussian("pc", -2.5e6, 1e6, -Inf, 0);
prior_sc = EKP.constrained_gaussian("sc", -2e-6, 1e-6, -Inf, 0);
prior = EKP.combine_distributions([prior_pc, prior_sc]);

ensemble_size = 2
n_iterations = 2
noise = 1.0 * EKP.I # Should work, but this should be covariance of each month from observation (ERA5)

caldir = "calibration_output"

hpc_kwargs = CAL.kwargs(time = 60, ntasks = 1, gpus_per_task = 1, cpus_per_task = 4)
CAL.calibrate(
    CAL.ClimaGPUBackend,
    ensemble_size,
    n_iterations,
    observations,
    noise,
    prior,
    caldir;
    hpc_kwargs,
    model_interface,
)
