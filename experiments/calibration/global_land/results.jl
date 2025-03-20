import EnsembleKalmanProcesses as EKP
using JLD2

## Priors
prior_pc = EKP.constrained_gaussian("pc", -2e6, 1e6, -Inf, 0);
prior_sc = EKP.constrained_gaussian("sc", 5e-6, 2e-6, 0, Inf);
prior_K_sat_plant = EKP.constrained_gaussian("K_sat_plant", 3e-8, 2e-8, 0.0, Inf);
prior_a = EKP.constrained_gaussian("a", 0.00196, 0.00049, 0.0, 0.00588);
prior_h_leaf = EKP.constrained_gaussian("h_leaf", 4.0, 2.0, 0.5, 9.0);
prior_α_snow = EKP.constrained_gaussian("α_snow", 0.6, 0.2, 0.0, 1.0); # only applies physical bounds
prior_α_soildry_scale = EKP.constrained_gaussian("α_soil_dry_scaler", 1.15, 0.1, 1.0, 1.3); # capped at 1 from below to ensure alpha_wet < alpha_dry, from above to ensure that α_soil_scaler * α_soil_dry_scaler <~1.6, which is when dry albedos start to exceed 1
prior_α_soil_scale = EKP.constrained_gaussian("α_soil_scaler", 1, 0.15, 0.7, 1.3); # new parameter
prior_τ_leaf_scale = EKP.constrained_gaussian("τ_leaf_scaler", 1, 0.15, 0.5, 1.5);
prior_α_leaf_scale = EKP.constrained_gaussian("α_leaf_scaler", 1, 0.15, 0.5, 1.5);
prior = EKP.combine_distributions([
                                   prior_pc,
                                   prior_sc,
                                   prior_h_leaf,
                                   prior_a,
                                   prior_K_sat_plant,
                                   prior_α_leaf_scale,
                                   prior_τ_leaf_scale,
                                   prior_α_soildry_scale,
                                   prior_α_snow,
                                   prior_α_soil_scale,
                                  ]);

path_to_results = "calibration_output_utki_copy/iteration_001/eki_file.jld2"

uki = JLD2.load_object(path_to_results)

errors = uki.error
normalized_errors = errors ./ errors[1] .* 100

names = ["pc", "sc", "h_leaf", "a", "K_sat_plant", "α_leaf_scaler", "τ_leaf_scaler", "α_soil_dry_scaler", "α_snow", "α_soil_scale",]

params = EKP.get_ϕ_mean.(Ref(prior), Ref(uki), collect(1:2))
initial_params = params[1]
calibrated_params = params[end]
relative_change = calibrated_params ./ initial_params .* 100

EKP.get_ϕ(prior, uki, 1)
EKP.get_ϕ(prior, uki, 2)

# Other things to do, to look at how calibration improved model - data

# Plot data vs season, and model vs season before and after calibration, for all locations, and all variables, maybe for a subset of locations for visibility (2 or 3)

# Global plots: data vs model, actual values, anomalies, RMSE, seasonality
# (use CI leaderboard figures, and ClimaLand leaderboard)




