"""
Shared prior definitions for DK-Sor single-site calibration.

Include this file (via `include`) in both `run_calibration.jl` and
`experiments/callmip_uq_dk_sor/emulate_sample.jl` to guarantee both scripts
use exactly the same 17-parameter priors.

The prior table below is transplanted from
`ar/calibrate_inversion_nee/experiments/calibration/configs/inversionnee_sifgpp_er_lhf_shf_lwu.jl`
as requested, while keeping the CalLMIP pipeline structure from
`rb/callmip_phase1a`.
"""

import EnsembleKalmanProcesses.ParameterDistributions as PD

function build_dk_sor_priors()
    priors_vec = [
    # P-model + moisture stress
    PD.constrained_gaussian("pmodel_cstar", 0.4126532780348918, 0.05, 0.2, 0.7),
    PD.constrained_gaussian("pmodel_β_c3", 68.19135215644542, 40.0, 10.0, 300.0),
    PD.constrained_gaussian("pmodel_β_c4", 53.92806481016041, 10.0, 5.0, 100.0),
    PD.constrained_gaussian("pmodel_α", 0.9758442675597523, 0.01, 0.85, 0.999),
    PD.constrained_gaussian("moisture_stress_c", 0.5675395275525799, 0.15, 0.05, 1.0),

    # DAMM heterotrophic respiration
    PD.constrained_gaussian("soilCO2_reference_rate", 7.885402642122577e-8, 5.0e-8, 1.0e-9, 2.0e-6),
    PD.constrained_gaussian("soilCO2_activation_energy", 65432.29096630613, 25000.0, 5000.0, 150000.0),
    PD.constrained_gaussian("michaelis_constant", 0.3668980114900562, 0.25, 0.001, 2.0),
    PD.constrained_gaussian("O2_michaelis_constant", 0.00018610344999027417, 1.0e-4, 1.0e-6, 1.0e-2),

    # JULES autotrophic respiration
    PD.constrained_gaussian("root_leaf_nitrogen_ratio", 2.3992909128669666, 1.0, 0.1, 6.0),
    PD.constrained_gaussian("relative_contribution_factor", 1.1321179749274621, 0.3, 0.0, 3.0),
    PD.constrained_gaussian("autotrophic_respiration_Q10", 2.0, 0.5, 1.0, 3.5),

    # Canopy turbulent / radiative transfer
    PD.constrained_gaussian("leaf_Cd", 0.08, 0.03, 0.0, Inf),
    PD.constrained_gaussian("canopy_z_0m_coeff", 0.13, 0.05, 0.0, 0.5),
    PD.constrained_gaussian("canopy_z_0b_coeff", 0.013, 0.005, 0.0, 0.05),
    PD.constrained_gaussian("canopy_d_coeff", 0.67, 0.2, 0.0, 1.0),
    PD.constrained_gaussian("canopy_K_lw", 1.0, 0.2, 0.0, 2.0),
    ]

    prior = PD.combine_distributions(priors_vec)
    return prior, priors_vec
end
