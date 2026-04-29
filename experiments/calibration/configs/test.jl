# Calibration configuration for the end-to-end test.
#
# This config is loaded by run_calibration.jl when the environment variable
# TEST_CALIBRATION is set. It calibrates a single soil parameter against lwu
# over a single sample/iteration so the test runs quickly.

"""Noise scalars for the covariance matrix of each observed variable.

These multiply the identity in `ScalarCovariance`, so they are variances.
Units are the square of the observational dataset units:
- `lwu`: (W m^-2)^2 (ERA5)
"""
const NOISE_SCALARS = Dict("lwu" => 1.0)

"""
    get_calibration_prior()

Return the combined prior for the test calibration.

Calibrates one parameter (true value is 0.96).
"""
function get_calibration_prior()
    priors =
        [EKP.constrained_gaussian("emissivity_bare_soil", 0.82, 0.12, 0.0, 2.0)]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["lwu"],
    minibatch_size = 1,
    n_iterations = 1,
    sample_date_ranges = [("2007-12-1", "2007-12-1")],
    extend = Dates.Month(3),
    spinup = Dates.Month(0),
    nelements = (180, 360, 15),
    output_dir = OUTPUT_DIR,
    rng_seed = 42,
    obs_vec_filepath = joinpath(
        "experiments",
        "calibration",
        "land_observation_vector.jld2",
    ),
    model_type = ClimaLand.LandModel,
)
