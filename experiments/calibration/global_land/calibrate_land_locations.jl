ENV["JULIA_WORKER_TIMEOUT"] = "1000.0"

using ClimaLand
dir = pkgdir(ClimaLand)
import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL
using Distributed
addprocs(
    CAL.PBSManager(20), # n simulations in parallel
    q = "main", # not prio
    A = "UCIT0011",
    l_select = "1:ngpus=1:ncpus=4",
    l_walltime = "09:30:00",
)

@everywhere begin
    import ClimaCalibrate:
        forward_model, parameter_path, path_to_ensemble_member
    import ClimaCalibrate as CAL

    using ClimaLand
    caldir = "calibration_output_new"
    dir = pkgdir(ClimaLand)
    nelements = (50, 10) # resolution for model and era5

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
prior_sc = EKP.constrained_gaussian("sc", 5e-6, 2e-6, 0, Inf);
prior_K_sat_plant = EKP.constrained_gaussian("K_sat_plant", 3e-8, 2e-8, 0.0, Inf);
prior_a = EKP.constrained_gaussian("a", 0.00196, 0.00049, 0.0, 0.00588);
prior_h_leaf = EKP.constrained_gaussian("h_leaf", 4.0, 2.0, 0.5, 9.0);
prior_α_snow = EKP.constrained_gaussian("α_snow", 0.7, 0.08, 0.5, 1.0);
prior_α_soildry_scale = EKP.constrained_gaussian("α_soil_dry_scaler", 1.1, 0.05, 1.0, 1.35);
prior_τ_leaf_scale = EKP.constrained_gaussian("τ_leaf_scaler", 1.1, 0.05, 1.0, 1.35);
prior_α_leaf_scale = EKP.constrained_gaussian("α_leaf_scaler", 1.1, 0.05, 1.0, 1.35);
prior = EKP.combine_distributions([prior_pc,
                                   prior_sc,
                                   prior_h_leaf,
                                   prior_a,
                                   prior_K_sat_plant,
                                   prior_α_leaf_scale,
                                   prior_τ_leaf_scale,
                                   prior_α_soildry_scale,
                                   prior_α_snow]);

ensemble_size = 40 # ideally 40
n_iterations = 5 # try 5, ideally 10
n_locations = length(training_locations)
nelements = (50, 10) # resolution for model and era5
l_obs = 4*4*n_locations # 4 seasons * 4 variables * n_locations

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
