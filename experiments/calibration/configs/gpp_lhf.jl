# Calibration configuration for GPP and LHF
#
# This file is included by run_calibration.jl and defines:
#   - CALIBRATE_CONFIG (CalibrateConfig)
#   - NOISE_SCALARS (per-variable covariance scaling)
#   - get_calibration_prior() (combined prior distribution)
#
# To create a new calibration experiment, copy this file and adjust the
# settings below.

"""Noise scalars for the covariance matrix of each observed variable."""
const NOISE_SCALARS = Dict("gpp" => 2.0, "lhf" => 20.0)

"""
    get_calibration_prior()

Return the combined prior distribution for the calibration parameters.

Calibrates 5 parameters: 4 P-model + 1 soil moisture stress.
Ensemble size for TransformUnscented: 5 * 2 + 1 = 11 members.

Note: ϕ0_c3/ϕ0_c4 are not calibrated because temperature_dep_yield = true
uses the quadratic coefficients (ϕa0, ϕa1, ϕa2) instead.
"""
function get_calibration_prior()
    priors = [
        EKP.constrained_gaussian("pmodel_cstar", 0.41, 0.05, 0.2, 0.7),
        EKP.constrained_gaussian("pmodel_β_c3", 146.0, 40.0, 50.0, 300.0),
        EKP.constrained_gaussian("pmodel_β_c4", 16.222, 5.0, 5.0, 40.0),
        EKP.constrained_gaussian("pmodel_α", 0.933, 0.02, 0.85, 0.999),
        EKP.constrained_gaussian("moisture_stress_c", 0.27, 0.15, 0.05, 1.0),
    ]
    return EKP.combine_distributions(priors)
end

if !TEST_CALIBRATION
    const CALIBRATE_CONFIG = CalibrateConfig(;
        short_names = ["gpp", "lhf"],
        minibatch_size = 4,
        n_iterations = 10,
        # 10 yearly samples: each covers DJF, MAM, JJA (Dec 1 -> Sep 1)
        # Ending at Sep 1 avoids DJF over-representation across consecutive
        # samples that would occur with Dec-to-Dec ranges.
        # with extend = Month(3), simulation runs through Nov 30 of year+1
        sample_date_ranges = [
            ("$(year)-12-1", "$(year+1)-9-1") for year in 2001:2010
        ],
        extend = Dates.Month(3),
        spinup = Dates.Year(1),
        nelements = (180, 360, 15),
        output_dir = "/glade/derecho/scratch/arenchon/calibration_gpp_lhf",
        rng_seed = 42,
        obs_vec_filepath = "experiments/calibration/land_observation_vector.jld2",
        model_type = ClimaLand.LandModel,
    )
else
    @info "Using calibration config for test calibration"
    const CALIBRATE_CONFIG = CalibrateConfig(;
        short_names = ["gpp", "lhf"],
        minibatch_size = 1,
        n_iterations = 1,
        sample_date_ranges = [("2007-12-1", "2007-12-1")],
        extend = Dates.Month(3),
        spinup = Dates.Month(0),
        nelements = (180, 360, 15),
        output_dir = "/glade/derecho/scratch/arenchon/calibration_gpp_lhf",
        rng_seed = 42,
        obs_vec_filepath = "experiments/calibration/land_observation_vector.jld2",
        model_type = ClimaLand.LandModel,
    )
end
