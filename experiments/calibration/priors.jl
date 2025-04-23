#prior_pc = EKP.constrained_gaussian("pc", -2e6, 1e6, -1e7, -1e5);
#prior_sc = EKP.constrained_gaussian("sc", 5e-6, 2e-6, 1e-7, 1e-5);
#prior_K_sat_plant =
#    EKP.constrained_gaussian("K_sat_plant", 7e-8, 2e-8, 1e-9, 1e-7);
#prior_a = EKP.constrained_gaussian("a", 0.00196, 0.00049, 0.0001, 0.00588);
## prior_α_snow = EKP.constrained_gaussian("α_snow", 0.65, 0.2, 0.4, 1.0);
#prior_α_soildry_scale =
#    EKP.constrained_gaussian("α_soil_dry_scaler", 1.15, 0.05, 1.0, 1.3); # capped at 1 from below to ensure alpha_wet < alpha_dry, from above to ensure that α_soil_scaler * α_soil_dry_scaler <~1.6, which is when dry albedos start to exceed 1
#prior_α_soil_scale =
#    EKP.constrained_gaussian("α_soil_scaler", 1, 0.15, 0.7, 1.3);
#prior_τ_leaf_scale =
#    EKP.constrained_gaussian("τ_leaf_scaler", 1, 0.15, 0.5, 1.5);
#prior_α_leaf_scale =
#    EKP.constrained_gaussian("α_leaf_scaler", 1, 0.15, 0.5, 1.5);

prior_α_0 = EKP.constrained_gaussian("α_0", 0.5, 0.2, 0.2, 0.8);
prior_α_horizon = EKP.constrained_gaussian("α_horizon", 0.85, 0.1, 0.7, 1);
prior_k = EKP.constrained_gaussian("k", 10, 5, 2, 25);

prior_beta_snow = EKP.constrained_gaussian("beta_snow", 0.4, 0.2, 0.1, 0.8);
prior_x0_snow = EKP.constrained_gaussian("x0_snow", 0.4, 0.2, 0.1, 0.8);
prior_gamma_snow = EKP.constrained_gaussian("gamma_snow", 0.1, 0.05, 0.01, 0.2);
prior_beta_0 = EKP.constrained_gaussian("beta_0", 1.8, 0.1, 1.0, 2.0);
#prior_beta_min = EKP.constrained_gaussian("beta_min", 0.85, 0.1, 0.7, 1);
prior_z0_snow = EKP.constrained_gaussian("z0_snow", 0.106, 0.05, 0.01, 0.3);

prior = EKP.combine_distributions([
    #    prior_pc,
    #    prior_sc,
    #    prior_a,
    #    prior_K_sat_plant,
    #    prior_α_leaf_scale,
    #    prior_τ_leaf_scale,
    #    prior_α_soildry_scale,
    #    # prior_α_snow,
    #    prior_α_soil_scale,
    prior_α_0,
    prior_α_horizon,
    prior_k,
    prior_beta_snow,
    prior_x0_snow,
    prior_gamma_snow,
    prior_beta_0,
    #    prior_beta_min,
    prior_z0_snow,
]);
