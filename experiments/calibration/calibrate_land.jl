ENV["JULIA_WORKER_TIMEOUT"] = "1000.0"

using ClimaLand
dir = pkgdir(ClimaLand)
import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL
FT = Float64
using Distributed
using Random
rng_seed = 2
rng = Random.MersenneTwister(rng_seed)

addprocs(
    CAL.PBSManager(3), # n simulations in parallel
    q = "main", # not prio
    A = "UCIT0011",
    l_select = "1:ngpus=1:ncpus=4",
    l_walltime = "11:30:00",
    l_job_priority = "premium",
)

@everywhere begin
    import ClimaCalibrate:
        forward_model, parameter_path, path_to_ensemble_member
    import ClimaCalibrate as CAL

    using ClimaLand
    caldir = "calibration_output_utki"
    dir = pkgdir(ClimaLand)
    nelements = (50, 10) # resolution for model and era5

    include(joinpath(dir, "experiments/calibration/forward_model.jl"))
end

import EnsembleKalmanProcesses as EKP

## Priors
#prior_pc = EKP.constrained_gaussian("pc", -2e6, 1e6, -Inf, 0);
#prior_sc = EKP.constrained_gaussian("sc", 5e-6, 2e-6, 0, Inf);
#prior_K_sat_plant =
#    EKP.constrained_gaussian("K_sat_plant", 8e-8, 2e-8, 0.0, Inf);
#prior_a = EKP.constrained_gaussian("a", 0.00196, 0.00049, 0.0, 0.00588);
#prior_h_leaf = EKP.constrained_gaussian("h_leaf", 4.0, 2.0, 0.5, 9.0);
prior_α_snow = EKP.constrained_gaussian("α_snow", 0.6, 0.2, 0.0, 1.0); # only applies physical bounds
#prior_α_soildry_scale =
#    EKP.constrained_gaussian("α_soil_dry_scaler", 1.15, 0.05, 1.0, 1.3); # capped at 1 from below to ensure alpha_wet < alpha_dry, from above to ensure that α_soil_scaler * α_soil_dry_scaler <~1.6, which is when dry albedos start to exceed 1
#prior_α_soil_scale =
#    EKP.constrained_gaussian("α_soil_scaler", 1, 0.15, 0.7, 1.3); # new parameter
#prior_τ_leaf_scale =
#    EKP.constrained_gaussian("τ_leaf_scaler", 1, 0.15, 0.5, 1.5);
#prior_α_leaf_scale =
#    EKP.constrained_gaussian("α_leaf_scaler", 1, 0.15, 0.5, 1.5);
prior = EKP.combine_distributions([
#    prior_pc,
#    prior_sc,
#    prior_h_leaf,
#    prior_a,
#    prior_K_sat_plant,
#    prior_α_leaf_scale,
#    prior_τ_leaf_scale,
#    prior_α_soildry_scale,
    prior_α_snow,
#    prior_α_soil_scale,
]);

# ensemble_size = 40 # ideally 40
n_iterations = 10 # try 5, ideally 10
nelements = (50, 10) # resolution for model and era5

# observationseries
include(joinpath(dir, "experiments/calibration/observationseries_era5.jl"))

n_locations = length(training_locations)
n_variables = 1 # SWU
n_time_points = 4 # 4 seasons (and not, for example, 12 months)
l_obs = n_time_points * n_variables * n_locations # 4 seasons * 3 variables * n_locations
ekp_process = EKP.Unscented(prior)
ensemble_size = ekp_process.N_ens

caldir = "calibration_output_utki"

include(joinpath(dir, "experiments/calibration/observation_map.jl"))

utki = EKP.EnsembleKalmanProcess(
    # EKP.construct_initial_ensemble(rng, prior, ensemble_size),
    observationseries,
    # EKP.Inversion();
    EKP.TransformUnscented(prior, impose_prior = true);
    verbose = true,
    # localization_method=EKP.Localizers.NoLocalization(),
    rng = rng,
    # scheduler=EKP.DefaultScheduler(1.0)
)

CAL.calibrate(CAL.WorkerBackend, utki, n_iterations, prior, caldir)
