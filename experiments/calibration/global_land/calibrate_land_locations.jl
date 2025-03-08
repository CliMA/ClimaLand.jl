ENV["JULIA_WORKER_TIMEOUT"] = "300.0"

using ClimaLand
dir = pkgdir(ClimaLand)
import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL
using Distributed
addprocs(
    CAL.PBSManager(2), # does that mean only 2 GPU or 2 workers run simultaneously?
    q = "main", # not prio
    A = "UCIT0011",
    l_select = "ngpus=4",
    l_walltime = "10:00:00",
)

@everywhere begin
    import ClimaCalibrate:
        forward_model, parameter_path, path_to_ensemble_member
    import ClimaCalibrate as CAL

    using ClimaLand
    caldir = "calibration_output"
    dir = pkgdir(ClimaLand)
    n_elements = (50, 10) # resolution for model and era5

    include(
        joinpath(
            dir,
            "experiments/calibration/global_land/forward_model.jl",
        ),
    )
end

import EnsembleKalmanProcesses as EKP

## Priors
prior_pc = EKP.constrained_gaussian("pc", -2e6, 1e6, -Inf, 0);
prior_sc = EKP.constrained_gaussian("sc", 5e-6, 5e-4, 0, Inf);
prior = EKP.combine_distributions([prior_pc, prior_sc]);

ensemble_size = 6
n_iterations = 2
n_locations = 500 # number of random locations
nelements = (50, 10) # resolution for model and era5
l_obs = 10*4*n_locations # 10 months * 4 variables * n_locations

caldir = "calibration_output"

#include(joinpath(dir, "experiments/calibration/shared/hand_picked_locations.jl"))
include(joinpath(dir, "experiments/calibration/shared/rand_locations.jl"))
include(joinpath(dir, "experiments/calibration/shared/make_training_locations.jl"))
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
