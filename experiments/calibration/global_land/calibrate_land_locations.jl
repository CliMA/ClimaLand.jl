ENV["JULIA_WORKER_TIMEOUT"] = "300.0"

using ClimaLand
dir = pkgdir(ClimaLand)
import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL
using Distributed
addprocs(
    CAL.PBSManager(8), # does that mean only 2 GPU or 2 workers run simultaneously?
    q = "main", # not prio
    A = "UCIT0011",
    l_select = "ngpus=8",
    l_walltime = "10:00:00",
)

@everywhere begin
    import ClimaCalibrate:
        forward_model, parameter_path, path_to_ensemble_member
    import ClimaCalibrate as CAL

    using ClimaLand
    caldir = "calibration_output_new"
    dir = pkgdir(ClimaLand)
    n_elements = (50, 10) # resolution for model and era5

    include(
        joinpath(
            dir,
            "experiments/calibration/global_land/forward_model.jl",
        ),
    )
end

include(joinpath(dir, "experiments/calibration/shared/stable_training_locations.jl"))

import EnsembleKalmanProcesses as EKP

## Priors
prior_pc = EKP.constrained_gaussian("pc", -2e6, 1e6, -Inf, 0);
prior_sc = EKP.constrained_gaussian("sc", 5e-6, 5e-4, 0, Inf);
prior = EKP.combine_distributions([prior_pc, prior_sc]);

ensemble_size = 16
n_iterations = 4
n_locations = length(training_locations)
nelements = (50, 10) # resolution for model and era5
l_obs = 12*4*n_locations # 12 months * 4 variables * n_locations

caldir = "calibration_output_new"

include(
    joinpath(
        dir,
        "experiments/calibration/global_land/observation_map_locations.jl",
    ),
)

include(joinpath(dir, "experiments/calibration/shared/observation_data_locations.jl"))
using LinearAlgebra
noise = Diagonal(noise_era5) # has to be a Matrix. Could also have covariances between variables.

CAL.calibrate(
    CAL.WorkerBackend,
    ensemble_size,
    n_iterations,
    observations,
    noise,
    prior,
    caldir,
)
