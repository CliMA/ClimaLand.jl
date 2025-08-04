prior_ϵ_soil = EKP.constrained_gaussian("ϵ_soil", 0.97, 0.02, 0.9, 1.0);
prior_ϵ_canopy = EKP.constrained_gaussian("ϵ_canopy", 0.97, 0.02, 0.9, 1.0);
prior = EKP.combine_distributions([prior_ϵ_soil, prior_ϵ_canopy]);
