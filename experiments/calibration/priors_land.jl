prior_α_0 = EKP.constrained_gaussian("α_0", 0.7, 0.1, 0.5, 0.95);
prior_beta_snow_cover = EKP.constrained_gaussian("beta_snow_cover", 1.77, 0.2, 0.1, 2.5);
prior_z0_snow_cover = EKP.constrained_gaussian("z0_snow_cover", 0.106, 0.03, 0.0, 0.3);
prior = EKP.combine_distributions([prior_α_0, prior_beta_snow_cover, prior_z0_snow_cover]);
