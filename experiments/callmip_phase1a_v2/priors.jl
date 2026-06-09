"""
Shared prior definitions for CalLMIP Phase 1a DK-Sor single-site calibration.

16-parameter set calibrating P-model, moisture stress, DAMM heterotrophic
respiration, JULES autotrophic respiration, and canopy turbulent/radiative
transfer parameters.

Key changes vs PR #1693 (rb/callmip_phase1a):
  - `autotrophic_respiration_Q10` removed — parameter does not exist in current
    ClimaLand.jl main.
  - `soilCO2_reference_rate` and `O2_michaelis_constant` now use [0, Inf]
    lognormal transforms instead of bounded Gaussians. The old bounded Gaussians
    placed the prior mean near the lower bound in transform space, causing EKP's
    constrained_gaussian to initialise with silently wrong Rh (~13× too high).
  - `soilCO2_reference_rate` mean set to ~6e-8 (empirically validated in single-
    column sweep to give correct summer Rh magnitude at DK-Sor).

Include this file in `run_calibration.jl`, `emulate_sample.jl`, and
`run_callmip_simulations.jl` to guarantee all scripts use identical priors.
"""

import EnsembleKalmanProcesses.ParameterDistributions as PD

function build_dk_sor_priors()
    priors_vec = [

        # ── P-model + moisture stress ──────────────────────────────────────────
        # Quantum-yield efficiency (Mengoli 2022; 1 - 1/T_acclim_timescale)
        PD.constrained_gaussian("pmodel_α", 0.9758, 0.01, 0.85, 0.999),

        # Unit cost ratio β for C3/C4 photosynthesis (Stocker 2020)
        PD.constrained_gaussian("pmodel_β_c3", 68.19, 40.0, 10.0, 300.0),
        PD.constrained_gaussian("pmodel_β_c4", 53.93, 10.0, 5.0, 100.0),

        # c* in pCO2 normalisation
        PD.constrained_gaussian("pmodel_cstar", 0.4127, 0.05, 0.2, 0.7),

        # Piecewise moisture stress c parameter
        PD.constrained_gaussian("moisture_stress_c", 0.5675, 0.15, 0.05, 1.0),

        # ── DAMM heterotrophic respiration ─────────────────────────────────────
        # V_ref_sx: maximum respiration rate at T_ref = 288.15 K  (kg C m⁻³ s⁻¹)
        # LOGNORMAL [0, Inf] — old bounded-Gaussian placed mean near lower bound
        # in transform space, giving initial Rh ~13× too high.
        # Mean ~6e-8 validated via single-column sweep to produce correct DK-Sor
        # summer Rh magnitude (~1.5 gC m⁻² d⁻¹).
        PD.constrained_gaussian("soilCO2_reference_rate", 6.0e-8, 5.0e-8, 0.0, Inf),

        # Activation energy (J mol⁻¹)
        PD.constrained_gaussian(
            "soilCO2_activation_energy",
            65432.0, 25000.0, 5000.0, 150000.0,
        ),

        # Michaelis-Menten CO2 half-saturation constant (m³/m³)
        PD.constrained_gaussian(
            "michaelis_constant", 0.3669, 0.25, 0.001, 2.0,
        ),

        # O2 half-saturation constant (m³/m³)
        # LOGNORMAL [0, Inf] — same reasoning as soilCO2_reference_rate.
        PD.constrained_gaussian(
            "O2_michaelis_constant", 0.000186, 1.0e-4, 0.0, Inf,
        ),

        # ── JULES autotrophic respiration ──────────────────────────────────────
        # Ratio of root nitrogen to leaf nitrogen
        PD.constrained_gaussian(
            "root_leaf_nitrogen_ratio", 2.399, 1.0, 0.1, 6.0,
        ),

        # Growth-respiration relative contribution factor
        PD.constrained_gaussian(
            "relative_contribution_factor", 1.132, 0.3, 0.0, 3.0,
        ),

        # ── Canopy turbulent / roughness / LW ─────────────────────────────────
        PD.constrained_gaussian("leaf_Cd",            0.08, 0.03, 0.0, Inf),
        PD.constrained_gaussian("canopy_z_0m_coeff",  0.13, 0.05, 0.0, 0.5),
        PD.constrained_gaussian("canopy_z_0b_coeff",  0.013, 0.005, 0.0, 0.05),
        PD.constrained_gaussian("canopy_d_coeff",     0.67, 0.2,  0.0, 1.0),
        PD.constrained_gaussian("canopy_K_lw",        1.0,  0.2,  0.0, 2.0),
    ]

    prior = PD.combine_distributions(priors_vec)
    return prior, priors_vec
end
