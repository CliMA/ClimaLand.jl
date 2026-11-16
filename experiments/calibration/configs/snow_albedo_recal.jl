# DRAFT (not committed to ClimaLand): calibration configuration for the snow
# season. To be copied to experiments/calibration/configs/snow_albedo_recal.jl
# on the branch rb/recal-snow-albedo (see
# papers/2026/ClimaLand-TRENDY/memos/2026-11-15_recalibration_plan.md).
#
# This file is included by run_calibration.jl and defines:
#   - CALIBRATE_CONFIG (CalibrateConfig)
#   - NOISE_SCALARS (per-variable covariance scaling)
#   - get_calibration_prior() (combined prior distribution)
#
# Use:
#   export CALIBRATION_CONFIG=snow_albedo_recal.jl
#   export CLIMALAND_MODIS_SNOW_OBS_DIR=/glade/derecho/scratch/$USER/recal_obs
#   export RECAL_TIER=core            # or "full"
#   export RECAL_N_ITERATIONS=1       # pilot; then 9 with the same output dir
#   bash experiments/calibration/run_calibration.sh /glade/derecho/scratch/$USER/recal_snow_albedo
#
# Targets (the first four are those of #1926; the last three are new):
#   - lhf, shf, lwu: ERA5, seasonal means, W m^-2
#   - sif_gpp:       GOSIF GPP (inversion_nee artifact), seasonal means, g C m^-2 day^-1
#   - modis_swu:     MODIS MCD43C3 blue-sky albedo x ERA5 SW down, MONTHLY Feb-Jun,
#                    land 40-90N without glaciers, W m^-2 (model: swu)
#   - modis_swu_ice: the same on the ice sheets (Greenland Feb-Jun, Antarctica Dec-Feb),
#                    MONTHLY, W m^-2 (model: swu)
#   - modis_scf:     MODIS snow cover (ILAMB), MONTHLY Dec-Jun, open land 30-90N,
#                    unitless (model: snowc)
#
# Code base: main + rb/snow-albedo-diffuse-zenith (zenith term on direct light
# only) + the calibration changes (target aliases and monthly targets).

"""Noise scalars for the covariance matrix of each observed variable.

These multiply the identity in `ScalarCovariance` (latitude-weighted), so they
are variances in the square of the target units. They are set so that each
target's typical squared misfit per point divided by its noise is about 2 with
the current defaults (seasonal RMSE: LHF ~14, SHF ~14, LWU ~4 W m^-2, GPP ~1.1
g C m^-2 day^-1; Feb-May SWU vs MODIS 15.4 W m^-2 and monthly open-land snow
cover 0.156 in the smoke test on main's 2-yr run), so no
target dominates the loss by its units. The ice-sheet SWU target is weighted
down (larger variance): ERA5 SW down and MODIS at low sun are less certain there.
"""
const NOISE_SCALARS = Dict(
    "lhf" => 100.0,
    "shf" => 100.0,
    "lwu" => 9.0,
    "sif_gpp" => 0.6,
    "modis_swu" => 100.0,
    "modis_swu_ice" => 225.0,
    "modis_scf" => 0.012,
)

const RECAL_TIER = get(ENV, "RECAL_TIER", "core")

"""
    get_calibration_prior()

Return the combined prior distribution for the calibration parameters.

Tier "core" (15 parameters, ensemble 31 with TransformUnscented):
  - stem reflectances (PAR, NIR), centred at 1.4x the woody-plant values;
  - snow albedo (ZenithAngleAlbedoModel: alpha_0, delta_alpha, k, beta, x0);
  - Wu & Wu snow cover (beta_0, z0), centred between the default and the
    offline MODIS fit (beta_0 1.08, z0 0.068);
  - the #1926 parameters that trade off with the stem reflectances in the
    surface energy balance: leaf_Cd, canopy_K_lw, the sub-canopy wind
    extinction and minimum reference height, the soil scalar roughness, and
    moisture_stress_c (GPP/LHF).
Tier "full" (25 parameters, ensemble 51) adds the rest of the #1926 set as it
exists on main: the P-model costs and acclimation, the canopy roughness, the
SAI coefficients (split into tree/grass on main) and the dry soil layer depth.
"""
function get_calibration_prior()
    priors = [
        # Stems: #1926 ended at 2x the woody-plant values (0.32 / 0.78) with
        # the scale at its upper bound; 1x fails SHF/LWU; ~1.4x matches MODIS.
        EKP.constrained_gaussian(
            "stem_PAR_reflectance",
            0.22,
            0.06,
            0.08,
            0.40,
        ),
        EKP.constrained_gaussian(
            "stem_NIR_reflectance",
            0.55,
            0.12,
            0.20,
            0.85,
        ),
        # Snow albedo, f(x) [(1 - fd)(alpha_0 + delta_alpha exp(-k cos θ)) + fd α_diff],
        # f(x) = min(1 - beta (ρ_snow/ρ_liq - x0), 1). Centred on the defaults.
        EKP.constrained_gaussian("alpha_0", 0.59, 0.08, 0.40, 0.85),
        EKP.constrained_gaussian("delta_alpha", 0.40, 0.10, 0.0, 0.60),
        EKP.constrained_gaussian("k", 1.96, 0.6, 0.1, 5.0),
        # beta is documented in [0, 1]; the formula holds for beta > 1 (darker
        # with density), and an upper bound at 1 would sit 0.03 from the mean.
        EKP.constrained_gaussian("beta", 0.97, 0.3, 0.0, 2.0),
        EKP.constrained_gaussian("x0", 0.2, 0.08, 0.05, 0.5),
        # Wu & Wu snow cover. At 1 deg, beta_scf = max(beta_0 + 0.5 gamma,
        # beta_min = 1); the lower bound keeps beta_scf above the floor.
        EKP.constrained_gaussian("beta_0", 1.4, 0.3, 0.96, 3.0),
        EKP.constrained_gaussian("z0", 0.085, 0.03, 0.01, 0.5),
        # #1926 parameters that trade off with the stems (recentred on main)
        EKP.constrained_gaussian("leaf_Cd", 1.17, 0.4, 0.01, 3.0),
        EKP.constrained_gaussian("canopy_K_lw", 0.999, 0.2, 0.3, 2.0),
        # #1926 posterior 0.0556 sat on the lower bound 0.05: bound lowered
        EKP.constrained_gaussian(
            "canopy_subcanopy_wind_extinction_coefficient",
            0.0556,
            0.03,
            0.005,
            2.0,
        ),
        # #1926 posterior 2.86 sat near the upper bound 3: bound raised
        EKP.constrained_gaussian(
            "canopy_subcanopy_min_reference_height",
            2.86,
            0.8,
            0.5,
            6.0,
        ),
        EKP.constrained_gaussian(
            "soil_scalar_roughness_length",
            0.00302,
            0.002,
            1.0e-5,
            0.1,
        ),
        EKP.constrained_gaussian("moisture_stress_c", 0.759, 0.1, 0.05, 1.0),
    ]
    if RECAL_TIER == "full"
        append!(
            priors,
            [
                # main (#1901) reset these to the pyrealm values; freeing them
                # again needs the agreement of whoever owns the P-model defaults
                EKP.constrained_gaussian("pmodel_cstar", 0.41, 0.05, 0.2, 0.7),
                EKP.constrained_gaussian(
                    "pmodel_β_c3",
                    146.0,
                    40.0,
                    2.0,
                    300.0,
                ),
                EKP.constrained_gaussian(
                    "pmodel_β_c4",
                    16.222,
                    5.0,
                    1.0,
                    100.0,
                ),
                EKP.constrained_gaussian(
                    "pmodel_α",
                    0.0113,
                    0.005,
                    0.002,
                    0.33,
                ),
                # #1926 posteriors at the lower bound 0; recentred off it
                EKP.constrained_gaussian(
                    "canopy_z_0m_coeff",
                    0.01,
                    0.005,
                    0.0,
                    0.5,
                ),
                EKP.constrained_gaussian(
                    "canopy_z_0b_coeff",
                    0.005,
                    0.0025,
                    0.0,
                    0.5,
                ),
                EKP.constrained_gaussian(
                    "canopy_z_0min",
                    0.0257,
                    0.02,
                    1.0e-4,
                    0.3,
                ),
                EKP.constrained_gaussian(
                    "SAI_coeff_tree",
                    0.26,
                    0.08,
                    0.0,
                    1.0,
                ),
                EKP.constrained_gaussian(
                    "SAI_coeff_grass",
                    0.05,
                    0.03,
                    0.0,
                    1.0,
                ),
                EKP.constrained_gaussian(
                    "maximum_dry_soil_layer_depth",
                    0.00296,
                    0.0015,
                    1.0e-4,
                    0.1,
                ),
            ],
        )
    end
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = [
        "lhf",
        "shf",
        "lwu",
        "sif_gpp",
        "modis_swu",
        "modis_swu_ice",
        "modis_scf",
    ],
    # As #1926: 16 yearly samples (Dec 2003 - Nov 2019) in 2-year minibatches
    # with a 1-year spin-up. Each sample covers Dec 1 -> Sep 1 (seasonal
    # targets DJF, MAM, JJA, SON with extend = 3 months; monthly targets
    # Dec-Jun). A forward run is 1 + 2 years.
    minibatch_size = 2,
    # 9 as #1926. RECAL_N_ITERATIONS=1 runs a one-iteration pilot; rerunning
    # with the same output directory and 9 continues from where it stopped.
    n_iterations = parse(Int, get(ENV, "RECAL_N_ITERATIONS", "9")),
    sample_date_ranges = [
        ("$(year)-12-1", "$(year+1)-9-1") for year in 2003:2018
    ],
    extend = Dates.Month(3),
    spinup = Dates.Year(1),
    nelements = (180, 360, 15),
    output_dir = OUTPUT_DIR,
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector_snow_albedo_recal.jld2",
    model_type = ClimaLand.LandModel,
)
