ENV["JULIA_WORKER_TIMEOUT"] = "300.0"

using ClimaLand
dir = pkgdir(ClimaLand)
import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL
using Distributed
addprocs(
         CAL.PBSManager(2),
         q = "main", # not prio
         A = "UCIT0011",
         l_select = "ngpus=1",
         l_walltime = "10:00:00",
        )

@everywhere begin
    import ClimaCalibrate:
        forward_model, parameter_path, path_to_ensemble_member
    import ClimaCalibrate as CAL

    using ClimaLand
    caldir = "calibration_output"
    dir = pkgdir(ClimaLand)

    include(joinpath(dir,"experiments/calibration/global_land/forward_model_land.jl"))
end

n_variables = 4 # lhf, shf, lwu, swu
length_one_var = 12 # 12 data point (months average) for 1 location
n_locations = 25 #
single_member_dims = (n_variables*length_one_var*n_locations,)

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

ensemble_size = 10
n_iterations = 2
noise = 1.0 * EKP.I # Should work, but this should be covariance of each month from observation (ERA5)

caldir = "calibration_output"

# hpc_kwargs = CAL.kwargs(time = 60, ntasks = 1, gpus_per_task = 1, cpus_per_task = 4)
CAL.calibrate(
    CAL.WorkerBackend,
    ensemble_size,
    n_iterations,
    observations,
    noise,
    prior,
    caldir, #;
    # hpc_kwargs,
    # model_interface,
)
