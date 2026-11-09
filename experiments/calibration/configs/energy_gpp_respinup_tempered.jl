# Calibration of the global land model against seasonal means of ERA5 energy
# fluxes (lhf, shf, swu, lwu) and GOSIF GPP, from the initial conditions of a
# new spin-up, with the noise inflated for correlated errors and a fixed UKI
# step.
#
# This file is included by run_calibration.jl and defines:
#   - CALIBRATE_CONFIG (CalibrateConfig)
#   - NOISE_SCALARS (per-variable covariance scaling)
#   - EKP_SCHEDULER (learning-rate scheduler)
#   - get_calibration_prior() (combined prior distribution)
#
# To use this config, set the CALIBRATION_CONFIG env var when invoking
# run_calibration.sh, e.g.
#   CALIBRATION_CONFIG=energy_gpp_respinup_tempered.jl bash experiments/calibration/run_calibration.sh <output_dir>

import EnsembleKalmanProcesses.ParameterDistributions as PD

"""Noise scalars for the covariance matrix of each observed variable.

These multiply the identity in `ScalarCovariance`, so they are variances. The
standard deviation of the seasonal means of each variable is its area-weighted
RMSE at the prior centre (seasonal means from December 2003 to November 2005,
from the initial conditions of the `rosetta_spunup_ic` artifact), so that each
variable counts equally in relative terms: ERA5 lhf, shf, swu, lwu (W m^-2) and
GOSIF GPP (`sif_gpp`, g C m^-2 day^-1). GPP has data in fewer cells (15138)
than ERA5 (22041), so its standard deviation is reduced by the square root of
their ratio, which makes it count as much as one energy flux.

All variances are multiplied by `NOISE_INFLATION`, the number of data per
effectively independent datum: the errors of the model are correlated in space
and time, so the ~800,000 data of a minibatch carry about a thousand
independent pieces of information. Without it, the prior has no weight in the
UKI update (`impose_prior = true`), whose steps then overshoot to the bounds.
"""
const NOISE_INFLATION = 1000
const NOISE_SCALARS = let
    seasonal_rmse = Dict(
        "lhf" => 14.5,
        "shf" => 14.2,
        "swu" => 10.6,
        "lwu" => 5.1,
        "sif_gpp" => 1.48,
    )
    Dict(
        name =>
            NOISE_INFLATION *
            (name == "sif_gpp" ? rmse * sqrt(15138 / 22041) : rmse)^2 for
        (name, rmse) in seasonal_rmse
    )
end

# Fixed UKI step: with `impose_prior = true` the iterations then converge to the
# posterior of the inflated noise
const EKP_SCHEDULER = EKP.DefaultScheduler(1)

"""
    logit_gaussian(name, median, upper, σ)

Return a prior on (0, `upper`) that is Gaussian, with standard deviation `σ`,
in the logit of the parameter over `upper`, with median `median`. Near zero it
is close to a log-normal, so that it can span orders of magnitude.
"""
function logit_gaussian(name, median, upper, σ)
    x = median / upper
    return PD.ParameterDistribution(
        Dict(
            "distribution" => PD.Parameterized(PD.Normal(log(x / (1 - x)), σ)),
            "constraint" => PD.bounded(0.0, upper),
            "name" => name,
        ),
    )
end

"""
    get_calibration_prior()

Return the combined prior distribution for the calibration parameters, centered
on the values in `toml/default_parameters.toml` except for the canopy roughness
coefficients, whose defaults are near zero: their priors span both near-zero and
physical values (0.05-0.15 of the canopy height) within two standard deviations.

`stem_reflectance_scale` multiplies the default PAR and NIR reflectances of
stems (see `model_parameters_path`). The P-model parameters are not calibrated.

Calibrates 14 parameters; the TransformUnscented ensemble has 29 members.
"""
function get_calibration_prior()
    priors = [
        # Soil moisture stress: GPP and LHF
        EKP.constrained_gaussian("moisture_stress_c", 0.759, 0.3, 0.05, 3.0),
        EKP.constrained_gaussian(
            "wilting_potential",
            -153.0,
            75.0,
            -600.0,
            -30.0,
        ),
        # Canopy turbulent exchange: LHF and SHF
        EKP.constrained_gaussian("leaf_Cd", 1.17, 0.5, 0.01, 3.0),
        logit_gaussian("canopy_z_0m_coeff", 0.003, 0.3, 2.3),
        logit_gaussian("canopy_z_0b_coeff", 0.002, 0.3, 2.3),
        EKP.constrained_gaussian("canopy_z_0min", 0.0257, 0.015, 0.001, 0.3),
        # Turbulent exchange of the ground below the canopy
        EKP.constrained_gaussian(
            "canopy_subcanopy_wind_extinction_coefficient",
            0.0556,
            0.05,
            0.005,
            2.0,
        ),
        EKP.constrained_gaussian(
            "canopy_subcanopy_min_reference_height",
            2.86,
            1.0,
            0.5,
            10.0,
        ),
        # Bare soil evaporation and heat exchange
        logit_gaussian("maximum_dry_soil_layer_depth", 0.00296, 0.1, 1.0),
        logit_gaussian("soil_scalar_roughness_length", 0.00302, 0.1, 1.0),
        # Radiation: SWU and LWU (grasses keep SAI_coeff_grass)
        EKP.constrained_gaussian("SAI_coeff_tree", 0.26, 0.13, 0.02, 1.0),
        EKP.constrained_gaussian(
            "stem_reflectance_scale",
            1.0,
            0.15,
            0.2,
            1.25,
        ),
        EKP.constrained_gaussian("canopy_K_lw", 0.999, 0.3, 0.1, 3.0),
        EKP.constrained_gaussian("alpha_0", 0.59, 0.07, 0.3, 0.8),
    ]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["lhf", "shf", "swu", "lwu", "sif_gpp"],
    minibatch_size = 2,
    n_iterations = 8,
    # 16 yearly samples, Dec 2003 to Nov 2019; each covers DJF, MAM, JJA, SON
    # (Dec 1 to Sep 1, and extend = Month(3) runs through Nov 30)
    sample_date_ranges = [
        ("$(year)-12-1", "$(year+1)-9-1") for year in 2003:2018
    ],
    extend = Dates.Month(3),
    # Start on March 1, the season of the spun-up initial conditions
    spinup = Dates.Month(9),
    nelements = (180, 360, 15),
    output_dir = OUTPUT_DIR,
    rng_seed = 42,
    obs_vec_filepath = joinpath(
        OUTPUT_DIR,
        "obs_energy_gpp_respinup_tempered.jld2",
    ),
    model_type = ClimaLand.LandModel,
    ic_path = "/glade/derecho/scratch/arenchon/claude/respinup_calib/ic/rosetta_spunup_ic.nc",
)
