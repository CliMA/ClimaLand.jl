"""
Prior parameter distributions for DK-Sor CalLMIP Phase 1a calibration.

15-parameter set: P-model, moisture stress, DAMM heterotrophic respiration,
JULES autotrophic respiration, and canopy turbulent/radiative transfer.

Design decisions:
  - `pmodel_β_c4` excluded — DK-Sor is a pure C3 site (European beech).
  - `autotrophic_respiration_Q10` excluded — does not exist in current ClimaLand.
  - `soilCO2_reference_rate` and `O2_michaelis_constant` use lognormal [0, Inf]
    transforms. Bounded Gaussians near zero cause Rh ~13× too high in transform
    space (EKP initialises at bound midpoint, not prior mean).
  - DAMM means recentered from global calibration (Alexis, ar/calibrate_inversion_nee):
      soilCO2_reference_rate   6e-8   → validated vs summer Rh observations
      soilCO2_activation_energy 36897 → recentered on EKI iteration-3 result
      michaelis_constant        0.051 → tighter, recentered
  - Canopy params keep physical defaults (iter-3 means thrashed, e.g. d_coeff→0.01).

Include this file from every script that needs EKP priors to guarantee consistency.
"""

import EnsembleKalmanProcesses.ParameterDistributions as PD

function build_dk_sor_priors()
    priors_vec = [

        # ── P-model + moisture stress ──────────────────────────────────────────
        PD.constrained_gaussian("pmodel_α",          0.9758,  0.01,   0.85,   0.999),
        PD.constrained_gaussian("pmodel_β_c3",       85.06,   40.0,   10.0,   300.0),
        PD.constrained_gaussian("pmodel_cstar",      0.4127,  0.05,   0.2,    0.7),
        PD.constrained_gaussian("moisture_stress_c", 0.5675,  0.15,   0.05,   1.0),

        # ── DAMM heterotrophic respiration ─────────────────────────────────────
        # soilCO2_reference_rate: lognormal [0, Inf]; mean 6e-8 kg C m⁻³ s⁻¹
        # validated to give correct summer Rh at DK-Sor (~1.84 vs 1.80 gC/m²/d).
        PD.constrained_gaussian("soilCO2_reference_rate",  6.0e-8, 4.0e-8, 0.0, Inf),
        # Activation energy recentered on EKI result (old default 65432 too high).
        PD.constrained_gaussian("soilCO2_activation_energy", 36897.0, 20000.0, 5000.0, 150000.0),
        PD.constrained_gaussian("michaelis_constant",   0.051,   0.03,   0.001,  2.0),
        # O2_michaelis_constant: lognormal [0, Inf]
        PD.constrained_gaussian("O2_michaelis_constant", 1.9e-4, 1.5e-4, 0.0, Inf),

        # ── JULES autotrophic respiration ──────────────────────────────────────
        PD.constrained_gaussian("root_leaf_nitrogen_ratio",    2.399, 1.0, 0.1, 6.0),
        PD.constrained_gaussian("relative_contribution_factor", 1.132, 0.3, 0.0, 3.0),

        # ── Canopy aerodynamics / longwave ─────────────────────────────────────
        # canopy_K_lw and canopy_z_0m_coeff are tightly constrained by flux obs
        # (PR #1693). All others remain near physical defaults.
        PD.constrained_gaussian("leaf_Cd",           0.08,  0.03,  0.0, Inf),
        PD.constrained_gaussian("canopy_z_0m_coeff", 0.13,  0.05,  0.0, 0.5),
        PD.constrained_gaussian("canopy_z_0b_coeff", 0.013, 0.005, 0.0, 0.05),
        PD.constrained_gaussian("canopy_d_coeff",    0.67,  0.2,   0.0, 1.0),
        PD.constrained_gaussian("canopy_K_lw",       1.0,   0.2,   0.0, 2.0),
    ]

    prior = PD.combine_distributions(priors_vec)
    return prior, priors_vec
end
