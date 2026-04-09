"""
Shared prior definitions for DK-Sor single-site calibration.

Include this file (via `include`) in both `run_calibration.jl` and
`experiments/callmip_uq_dk_sor/emulate_sample.jl` to guarantee both scripts
use exactly the same 16-parameter priors — avoiding the fragility of
duplicating and manually keeping two copies in sync.

Call `build_dk_sor_priors()` to get `(prior, priors_vec)` where:
  - `prior`      is the `ParameterDistribution` for the combined prior
  - `priors_vec` is the Vector of individual per-parameter distributions

Prior names MUST match ClimaParams TOML keys, since ClimaCalibrate writes
parameter TOMLs using these names and LandParameterTypes.create_toml_dict reads them.

16 parameters total:
  9 canopy / conductance parameters
  3 DAMM soil-CO₂ parameters
  3 autotrophic respiration parameters (winter NEE bias fix)
  1 canopy heat capacity
"""

import EnsembleKalmanProcesses.ParameterDistributions as PD

function build_dk_sor_priors()
    priors_vec = [
        # ── Canopy / conductance ─────────────────────────────────────────────
        PD.constrained_gaussian("moisture_stress_c",              0.5,      0.3,      0.01,     5.0),
        PD.constrained_gaussian("pmodel_cstar",                   0.43,     0.15,     0.05,     2.0),
        PD.constrained_gaussian("pmodel_β",                      20.0,      8.0,      2.0,     80.0),  # prior mean lowered: old 51 → GPP 2× too large
        PD.constrained_gaussian("leaf_Cd",                        0.1,      0.05,     0.005,    1.0),
        PD.constrained_gaussian("canopy_z_0m_coeff",              0.10,     0.04,     0.01,    0.25),  # raised: z0m ≈ 0.1h for tall forest
        PD.constrained_gaussian("canopy_z_0b_coeff",              0.001,    0.0005,   1e-5,    0.01),
        PD.constrained_gaussian("canopy_d_coeff",                 0.65,     0.12,     0.30,    0.92),  # raised: d ≈ 0.65h; default 0.007 gives d=17cm (wrong for 25m forest)
        PD.constrained_gaussian("canopy_K_lw",                    0.85,     0.25,     0.1,     2.0),
        PD.constrained_gaussian("canopy_emissivity",              0.97,     0.02,     0.9,     1.0),
        # ── DAMM soil-CO₂ ───────────────────────────────────────────────────
        PD.constrained_gaussian("soilCO2_pre_exponential_factor", 25000.0,  10000.0,  1000.0,  200000.0),
        PD.constrained_gaussian("michaelis_constant",             0.01,     0.005,    1e-4,     0.1),
        PD.constrained_gaussian("O2_michaelis_constant",          0.01,     0.005,    1e-4,     0.1),
        # ── Autotrophic respiration (winter NEE fix) ─────────────────────────
        # soilCO2_activation_energy [J/mol]: Arrhenius Ea for DAMM; default causes
        # excessive cold-temperature suppression → winter respiration too low.
        PD.constrained_gaussian("soilCO2_activation_energy",      61000.0,  15000.0,  30000.0, 100000.0),
        # root_leaf_nitrogen_ratio (μr): Ra_root = Rd * μr * RAI; calibrated to ~5 → more winter Ra
        PD.constrained_gaussian("root_leaf_nitrogen_ratio",       1.0,      0.5,      0.1,      5.0),
        # stem_leaf_nitrogen_ratio (μs): Ra_stem = Rd * μs * ... (LAI-independent)
        PD.constrained_gaussian("stem_leaf_nitrogen_ratio",       0.1,      0.05,     0.01,     0.5),
        # ── Canopy heat capacity ─────────────────────────────────────────────
        # ac_canopy [J/m²/K]: DK-Sor is dense beech forest with high wood thermal mass;
        # default 2500 J/m²/K likely too low → ac_canopy calibrated to ~9000 J/m²/K
        PD.constrained_gaussian("ac_canopy",                      2500.0,   1500.0,   500.0,   10000.0),
    ]

    prior = PD.combine_distributions(priors_vec)
    return prior, priors_vec
end
