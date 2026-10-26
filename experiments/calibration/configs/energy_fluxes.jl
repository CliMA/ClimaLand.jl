# Calibration configuration for a single soil-emissivity parameter against
# energy fluxes (lwu, shf, lhf). This matches the production default that
# lived in run_calibration.jl on main prior to the config-split refactor.

"""Noise scalars for the covariance matrix of each observed variable.

These multiply the identity in `ScalarCovariance`, so they are variances.
Units are the square of the observational dataset units:
- `lwu`: (W m^-2)^2 (ERA5)
- `shf`: (W m^-2)^2 (ERA5)
- `lhf`: (W m^-2)^2 (ERA5)
"""
const NOISE_SCALARS = Dict("lwu" => 1.0, "shf" => 1.0, "lhf" => 1.0)

"""
    get_calibration_prior()

Return the combined prior distribution for the calibration parameters.

Calibrates one parameter: soil emissivity. True value is near 0.96.
"""
function get_calibration_prior()
    priors =
        [EKP.constrained_gaussian("emissivity_bare_soil", 0.82, 0.12, 0.0, 2.0)]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["lwu", "shf", "lhf"],
    minibatch_size = 1,
    n_iterations = 10,
    sample_date_ranges = [
        ("$(2000 + 2*i)-12-1", "$(2002 + 2*i)-9-1") for i in 0:9
    ], # 2000 to 2020
    extend = Dates.Month(3),
    spinup = Dates.Month(3),
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
