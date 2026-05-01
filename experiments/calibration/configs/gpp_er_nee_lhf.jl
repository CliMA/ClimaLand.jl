# Calibration configuration for GPP + ER + NEE + LHF
#
# This file is included by run_calibration.jl and defines:
#   - CALIBRATE_CONFIG (CalibrateConfig)
#   - NOISE_SCALARS (per-variable covariance scaling)
#   - get_calibration_prior() (combined prior distribution)
#
# To use this config, set the CALIBRATION_CONFIG env var when invoking
# run_calibration.sh, e.g.
#   CALIBRATION_CONFIG=gpp_er_nee_lhf.jl bash experiments/calibration/run_calibration.sh
#
# Observational targets:
#   - gpp: ILAMB FLUXCOM (gpp_FLUXCOM_gpp.nc),   g C m^-2 day^-1
#   - er:  ILAMB FLUXCOM (reco_FLUXCOM_reco.nc), g C m^-2 day^-1
#   - nee: ILAMB FLUXCOM (nee_FLUXCOM_nee.nc),   g C m^-2 day^-1
#   - lhf: ERA5,                                 W m^-2
#
# Prior structure: 5 P-model/moisture-stress parameters (from gpp_lhf.jl) +
# 4 DAMM heterotrophic-respiration parameters + 2 JULES autotrophic-
# respiration parameters (from er.jl), for 11 parameters total.
# Ensemble size for TransformUnscented: 11 * 2 + 1 = 23 members.
#
# Priors are recentered on the posterior from the previous gpp_er_nee_lhf
# calibration (now in toml/default_parameters.toml). Standard deviations are
# kept the same when feasible; bounds are expanded only for parameters whose
# posterior reached or approached the previous bound:
#   - pmodel_β_c3: lower 50 -> 10 (posterior 68.2139 with σ=40 requires lower < ~20)
#   - pmodel_β_c4: upper 40 -> 80 (posterior 39.575 hit upper bound)
#   - soilCO2_reference_rate: lower 1e-8 -> 1e-9 (drifted toward lower bound)
#   - root_leaf_nitrogen_ratio: lower 0.05 -> 0.01 (posterior 0.0626 hit lower bound)
#   - relative_contribution_factor: upper 1.0 -> 1.5 (posterior 0.844 approached upper bound)
# Standard deviations were reduced where the original σ became infeasible after
# recentering on a much smaller posterior mean (σ shrunk to preserve σ/μ ratio):
#   - O2_michaelis_constant: σ 0.004 -> 0.002 (μ dropped 0.005 -> 0.00255)
#   - root_leaf_nitrogen_ratio: σ 0.5 -> 0.04 (μ dropped 1.0 -> 0.0626)

"""Noise scalars for the covariance matrix of each observed variable.

These multiply the identity in `ScalarCovariance`, so they are variances.
Units are the square of the observational dataset units:
- `gpp`: (g C m^-2 day^-1)^2 (ILAMB FLUXCOM)
- `er`:  (g C m^-2 day^-1)^2 (ILAMB FLUXCOM reco)
- `nee`: (g C m^-2 day^-1)^2 (ILAMB FLUXCOM)
- `lhf`: (W m^-2)^2 (ERA5)
"""
const NOISE_SCALARS =
    Dict("gpp" => 2.5, "er" => 2.5, "nee" => 2.5, "lhf" => 100.0)

"""
    get_calibration_prior()

Return the combined prior distribution for the calibration parameters.

Calibrates 11 parameters: 5 P-model/moisture-stress + 4 DAMM + 2 JULES
autotrophic respiration.
"""
function get_calibration_prior()
    # Means recentered on the posterior from the previous gpp_er_nee_lhf
    # calibration (now reflected in toml/default_parameters.toml).
    # Standard deviations are unchanged. Bounds are expanded for parameters
    # whose posterior reached or approached the previous bound.
    priors = [
        # P-model + moisture stress
        EKP.constrained_gaussian("pmodel_cstar", 0.448094, 0.05, 0.2, 0.7),
        # pmodel_β_c3 lower bound expanded 50 -> 10 (μ=68.2139 with σ=40 requires lower < ~20)
        EKP.constrained_gaussian("pmodel_β_c3", 68.2139, 40.0, 10.0, 300.0),
        # pmodel_β_c4 prior posterior 39.575 hit upper bound 40 -> expand to 80
        EKP.constrained_gaussian("pmodel_β_c4", 39.575, 5.0, 5.0, 80.0),
        EKP.constrained_gaussian("pmodel_α", 0.964115, 0.02, 0.85, 0.999),
        EKP.constrained_gaussian(
            "moisture_stress_c",
            0.377838,
            0.15,
            0.05,
            1.0,
        ),
        # DAMM heterotrophic respiration
        # soilCO2_reference_rate lower bound expanded 1e-8 -> 1e-9
        EKP.constrained_gaussian(
            "soilCO2_reference_rate",
            2.17564e-7,
            2.0e-7,
            1.0e-9,
            2.0e-6,
        ),
        EKP.constrained_gaussian(
            "soilCO2_activation_energy",
            37357.0,
            25000.0,
            5000.0,
            150000.0,
        ),
        EKP.constrained_gaussian(
            "michaelis_constant",
            0.459997,
            0.25,
            0.001,
            2.0,
        ),
        # O2_michaelis_constant σ reduced 0.004 -> 0.002 (μ dropped 0.005 -> 0.00255 makes original σ infeasible)
        EKP.constrained_gaussian(
            "O2_michaelis_constant",
            0.00255006,
            0.002,
            1.0e-5,
            1.0e-1,
        ),
        # JULES autotrophic respiration
        # root_leaf_nitrogen_ratio posterior 0.0626 hit lower bound 0.05 -> expand to 0.01
        # σ reduced 0.5 -> 0.04 (μ dropped 1.0 -> 0.0626 makes original σ infeasible)
        EKP.constrained_gaussian(
            "root_leaf_nitrogen_ratio",
            0.0625785,
            0.04,
            0.01,
            3.0,
        ),
        # relative_contribution_factor posterior 0.844 approached upper bound 1.0 -> expand to 1.5
        EKP.constrained_gaussian(
            "relative_contribution_factor",
            0.844417,
            0.15,
            0.0,
            1.5,
        ),
    ]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["gpp", "er", "nee", "lhf"],
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
    # Use the longer ER spinup (1 year) so deep soil temperature and
    # moisture profiles settle before Hr is compared to FLUXCOM reco.
    spinup = Dates.Year(1),
    nelements = (180, 360, 15),
    output_dir = OUTPUT_DIR,
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector_gpp_er_nee_lhf.jld2",
    model_type = ClimaLand.LandModel,
)
