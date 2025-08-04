prior_α_0 = EKP.constrained_gaussian("α_0", 0.6, 0.2, 0.2, 0.95);
prior_Δα = EKP.constrained_gaussian("Δα", 0.2, 0.1, 0.0, 1.0);
prior_k = EKP.constrained_gaussian("k", 4, 2, 0.01, 25);
#prior_α_0 = EKP.constrained_gaussian("α_0", 0.7, 0.1, 0.5, 0.95);
#prior_beta_snow_cover = EKP.constrained_gaussian("beta_snow_cover", 1.77, 0.2, 0.1, 2.5);
#prior_z0_snow_cover = EKP.constrained_gaussian("z0_snow_cover", 0.106, 0.03, 0.0, 0.3);
prior = EKP.combine_distributions([prior_α_0, prior_k, prior_Δα]);
