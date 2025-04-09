prior_pc = EKP.constrained_gaussian("pc", -2e6, 1e6, -1e7, -1e5);
prior_sc = EKP.constrained_gaussian("sc", 5e-6, 2e-6, 1e-7, 1e-5);
prior_K_sat_plant =
    EKP.constrained_gaussian("K_sat_plant", 7e-8, 2e-8, 1e-9, 1e-7);
prior_a = EKP.constrained_gaussian("a", 0.00196, 0.00049, 0.0001, 0.00588);
prior_α_snow = EKP.constrained_gaussian("α_snow", 0.65, 0.2, 0.4, 1.0);
prior_α_soildry_scale =
    EKP.constrained_gaussian("α_soil_dry_scaler", 1.15, 0.05, 1.0, 1.3); # capped at 1 from below to ensure alpha_wet < alpha_dry, from above to ensure that α_soil_scaler * α_soil_dry_scaler <~1.6, which is when dry albedos start to exceed 1
prior_α_soil_scale =
    EKP.constrained_gaussian("α_soil_scaler", 1, 0.15, 0.7, 1.3);
prior_τ_leaf_scale =
    EKP.constrained_gaussian("τ_leaf_scaler", 1, 0.15, 0.5, 1.5);
prior_α_leaf_scale =
    EKP.constrained_gaussian("α_leaf_scaler", 1, 0.15, 0.5, 1.5);

prior = EKP.combine_distributions([
    prior_pc,
    prior_sc,
    prior_a,
    prior_K_sat_plant,
    prior_α_leaf_scale,
    prior_τ_leaf_scale,
    prior_α_soildry_scale,
    prior_α_snow,
    prior_α_soil_scale,
]);

