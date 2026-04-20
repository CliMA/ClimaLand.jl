# Calibration configuration for ecosystem respiration (ER = Hr + Ra)
#
# This file is included by run_calibration.jl and defines:
#   - CALIBRATE_CONFIG (CalibrateConfig)
#   - NOISE_SCALARS (per-variable covariance scaling)
#   - get_calibration_prior() (combined prior distribution)
#
# To use this config, set the CALIBRATION_CONFIG env var when invoking
# run_calibration.sh, e.g.
#   CALIBRATION_CONFIG=er.jl bash experiments/calibration/run_calibration.sh
#
# Observational target: ILAMB FLUXCOM reco (reco_FLUXCOM_reco.nc), in g C m^-2 day^-1.
#
# Autotrophic respiration one-sided constraint:
#   root_leaf_nitrogen_ratio (μr) and relative_contribution_factor (Rel) have
#   hard upper bounds at their defaults (1.0 and 0.25). No sample drawn from
#   these priors can exceed the defaults, so Ra cannot grow above its current
#   value during calibration. EKP requires μ + σ < ub strictly, so the prior
#   means sit below the defaults (0.7 and 0.15). Under the current setup
#   (PrescribedBiomassModel with default SAI=0, see snowy_land.jl), the stem
#   term in Rpm is gated off by heaviside(SAI), so Rpm = Rd · (β + μr·RAI/LAI)
#   and reducing μr monotonically reduces Ra. This prevents the optimizer
#   from compensating a shrunk Hr with an inflated Ra. If SAI is ever turned
#   on (via the maxLAI constructor in src/standalone/Vegetation/Canopy.jl:372),
#   this guarantee holds only in expectation, not strictly.

"""Noise scalars for the covariance matrix of each observed variable.

Units follow the observational dataset:
- `er`: g C m^-2 day^-1 (ILAMB FLUXCOM reco)
"""
const NOISE_SCALARS = Dict("er" => 3.0)

"""
    get_calibration_prior()

Return the combined prior distribution for the calibration parameters.

Calibrates 6 parameters: 4 heterotrophic (DAMM) + 2 autotrophic (JULES).
Ensemble size for TransformUnscented: 6 * 2 + 1 = 13 members.

DAMM uses the centered-Arrhenius form
    Vmax = V_ref_sx * exp(-Ea_sx/R * (1/T - 1/T_ref_sx))
with T_ref_sx = 288.15 K fixed, so V_ref_sx and Ea_sx are approximately
orthogonal.
"""
function get_calibration_prior()
    priors = [
        EKP.constrained_gaussian(
            "soilCO2_reference_rate",
            2.526e-7,
            1.0e-7,
            5.0e-8,
            5.0e-7,
        ),
        EKP.constrained_gaussian(
            "soilCO2_activation_energy",
            40000.0,
            15000.0,
            20000.0,
            80000.0,
        ),
        EKP.constrained_gaussian("michaelis_constant", 0.3, 0.2, 0.01, 1.0),
        EKP.constrained_gaussian(
            "O2_michaelis_constant",
            0.005,
            0.003,
            5.0e-4,
            5.0e-2,
        ),
        # One-sided (reduce-only) priors on autotrophic respiration.
        # EKP.constrained_gaussian requires μ + σ < ub strictly (and
        # μ - σ > lb). The hard upper bounds at the defaults (1.0 and 0.25)
        # prevent any sample from exceeding them, so Ra cannot grow above
        # its current value during calibration. The means are set below
        # the defaults so the (μ, σ) constraint holds while still leaving
        # room to explore downward.
        EKP.constrained_gaussian(
            "root_leaf_nitrogen_ratio",
            0.7,
            0.2,
            0.05,
            1.0,
        ),
        EKP.constrained_gaussian(
            "relative_contribution_factor",
            0.15,
            0.05,
            0.0,
            0.25,
        ),
    ]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["er"],
    minibatch_size = 2,
    n_iterations = 10,
    # 10 yearly samples: each covers DJF, MAM, JJA (Dec 1 -> Sep 1)
    # Ending at Sep 1 avoids DJF over-representation across consecutive
    # samples that would occur with Dec-to-Dec ranges.
    # with extend = Month(3), simulation runs through Nov 30 of year+1
    sample_date_ranges = [
        ("$(year)-12-1", "$(year+1)-9-1") for year in 2001:2010
    ],
    extend = Dates.Month(3),
    # Longer spinup than GPP (3 months) so that deep soil temperature and
    # moisture profiles settle before Hr is compared to FLUXCOM reco.
    spinup = Dates.Year(1),
    nelements = (180, 360, 15),
    output_dir = OUTPUT_DIR,
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector_er.jld2",
    model_type = ClimaLand.LandModel,
)
