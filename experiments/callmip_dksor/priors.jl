"""
Prior parameter distributions for the DK-Sor CalLMIP calibration.

11 parameters: 5 P-model / moisture-stress, 4 DAMM heterotrophic-respiration, and
2 JULES autotrophic-respiration. Prior means are recentered on a global multi-site
posterior (also reflected in `toml/default_parameters.toml`).

Notes:
  - The CalLMIP DK-Sor objective uses NEE, LHF, and SHF only (GPP/ER are not available
    separately), so only those three fluxes constrain the calibration.
  - `pmodel_β_c4` is inert at DK-Sor (pure C3 European beech); it is kept to match the
    11-parameter set but is effectively unconstrained here.
  - The autotrophic-respiration parameters control the LAI-independent maintenance
    respiration; their lower means reduce the default over-respiration, while a higher
    `soilCO2_reference_rate` raises the (otherwise too low) heterotrophic respiration.

Include this file from every script that needs the priors to guarantee consistency.
"""

import EnsembleKalmanProcesses.ParameterDistributions as PD

function build_dk_sor_priors()
    # Priors are recentered so the calibration optimum lies within the prior mass (not
    # the tail), which gives the emulator/MCMC well-posed support. Two DAMM upper bounds
    # are widened (soilCO2_reference_rate, soilCO2_activation_energy) to admit the
    # posterior; the resulting carbon partition is physically consistent (HR/GPP ~20-24%,
    # AR/GPP ~60-67%).
    priors_vec = [
        # ── P-model + moisture stress ──────────────────────────────────────────
        PD.constrained_gaussian("pmodel_cstar",      0.3648,  0.05,  0.2,   0.7),
        PD.constrained_gaussian("pmodel_β_c3",       159.4,   40.0,  10.0,  300.0),
        PD.constrained_gaussian("pmodel_β_c4",       39.575,  5.0,   5.0,   80.0),   # inert (C3 site)
        PD.constrained_gaussian("pmodel_α",          0.9443,  0.02,  0.85,  0.999),
        PD.constrained_gaussian("moisture_stress_c", 0.3098,  0.15,  0.05,  1.0),

        # ── DAMM heterotrophic respiration ─────────────────────────────────────
        # NEE alone constrains only ecosystem respiration (AR + HR), not the AR/HR split
        # — a flat likelihood ridge. The regularized UTKI posterior stays physical
        # (HR/GPP ~20%); an unregularized sampler can drift along this ridge, so the
        # deliverable uses the UTKI posterior.
        PD.constrained_gaussian("soilCO2_reference_rate",    1.72e-6,   5.0e-7,  1.0e-9, 5.0e-6),
        PD.constrained_gaussian("soilCO2_activation_energy", 57080.0,   25000.0, 5000.0, 150000.0),
        PD.constrained_gaussian("michaelis_constant",        0.6234,    0.25,    0.001,  2.0),
        # The O2 Michaelis constant is kept near 8.7e-4: centering at the much smaller
        # global value is infeasible for the logit-normal (mean too close to the lower
        # bound). O2 limitation is secondary at this well-drained site.
        PD.constrained_gaussian("O2_michaelis_constant",     0.0008681, 0.0004,  1.0e-5, 1.0e-1),

        # ── JULES autotrophic respiration ──────────────────────────────────────
        # autotrophic_respiration_Rd_ref is the base maintenance-respiration rate
        # (R_root = Rd_ref·RAI, R_stem = Rd_ref·μs·SAI). root_leaf_nitrogen_ratio is
        # inactive in the current AR formulation, so Rd_ref is calibrated instead. As
        # with soilCO2_reference_rate, AR lies in the NEE-unidentifiable ridge, so the
        # deliverable uses the UTKI posterior.
        PD.constrained_gaussian("autotrophic_respiration_Rd_ref", 1.02e-7, 6.0e-8, 1.0e-8, 1.0e-6),
        PD.constrained_gaussian("relative_contribution_factor",   0.305, 0.15, 0.0, 1.5),
        # Energy/turbulent parameters (e.g. leaf_Cd, emissivity_bare_soil) were evaluated
        # but excluded: the NEE-dominated objective pulled SHF away from its already-good
        # default. The unmatched SHF spring peak is structural, not parametric, so the set
        # stays at 11 carbon parameters; SHF is left at its (neutral) default and leaf NIR
        # optics remain fixed at the site level.
    ]

    prior = PD.combine_distributions(priors_vec)
    return prior, priors_vec
end

# Names + bounds of the 11 parameters (must match build_dk_sor_priors above).
const DK_SOR_BOUNDS = [
    ("pmodel_cstar", 0.2, 0.7), ("pmodel_β_c3", 10.0, 300.0), ("pmodel_β_c4", 5.0, 80.0),
    ("pmodel_α", 0.85, 0.999), ("moisture_stress_c", 0.05, 1.0),
    ("soilCO2_reference_rate", 1.0e-9, 5.0e-6), ("soilCO2_activation_energy", 5000.0, 150000.0),
    ("michaelis_constant", 0.001, 2.0), ("O2_michaelis_constant", 1.0e-5, 1.0e-1),
    ("autotrophic_respiration_Rd_ref", 1.0e-8, 1.0e-6), ("relative_contribution_factor", 0.0, 1.5),
]

"""
    build_dk_sor_eki_prior(eki_mean, eki_std)

Prior for the CES *Sample* step, anchored on the *Calibrate* (EKI) result: a per-parameter
`constrained_gaussian` at the EKI final-ensemble mean ± std, reusing the calibration bounds.
CES is Calibrate–Emulate–Sample, so the Sample step should build on the calibration. Using an
uninformative prior instead lets the emulator-MCMC drift along the NEE-unidentifiable AR/HR
ridge into an unphysical regime (AR > GPP, NEE source); anchoring on the EKI posterior
regularizes those sloppy directions and keeps the CES posterior physical.
"""
function build_dk_sor_eki_prior(eki_mean::AbstractVector, eki_std::AbstractVector)
    @assert length(eki_mean) == length(DK_SOR_BOUNDS) == length(eki_std)
    priors_vec = map(enumerate(DK_SOR_BOUNDS)) do (i, (nm, lo, hi))
        σ = max(eki_std[i], 1.0e-4 * (hi - lo))   # tiny floor to avoid a degenerate prior
        PD.constrained_gaussian(nm, eki_mean[i], σ, lo, hi)
    end
    return PD.combine_distributions(priors_vec), priors_vec
end
