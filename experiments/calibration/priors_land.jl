prior_α_0 = EKP.constrained_gaussian("α_0", 0.8, 0.1, 0.5, 0.95);
prior_beta_snow = EKP.constrained_gaussian("beta_snow", 0.4, 0.2, 0.1, 1.0);
prior = EKP.combine_distributions([prior_α_0, prior_beta_snow]);
