prior_α_0 = EKP.constrained_gaussian("α_0", 0.5, 0.2, 0.2, 0.8);
prior_Δα = EKP.constrained_gaussian("Δα", 0.2, 0.1, 0.0, 1.0);
prior_k = EKP.constrained_gaussian("k", 10, 5, 2, 25);

prior = EKP.combine_distributions([prior_α_0, prior_Δα, prior_k]);
