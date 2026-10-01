# Calibration configuration for a single soil-emissivity parameter against
# energy fluxes (lwu, shf, lhf). This matches the production default that
# lived in run_calibration.jl on main prior to the config-split refactor.

"""Noise scalars for the covariance matrix of each observed variable.

These multiply the identity in `ScalarCovariance`, so they are variances.
Units are the square of the observational dataset units:
- `scf`: unitless (MODIS)
"""
const NOISE_SCALARS = Dict("snowc" => 0.001)

"""
    get_calibration_prior()

Return the combined prior distribution for the calibration parameters.
"""
function get_calibration_prior()
    priors =
        [EKP.constrained_gaussian("beta_0", 1.8, 0.3, 0.1, 4.0),
	 EKP.constrained_gaussian("z0", 0.106, 0.03, 0.01, 0.5),
         EKP.constrained_gaussian("holding_capacity_of_water_in_snow", 0.08, 0.03, 0.0, 0.2),
         EKP.constrained_gaussian("wet_snow_hydraulic_conductivity", 5e-4, 3e-4, 0.0, 1e-2)]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["snowc",],
    minibatch_size = 1,
    n_iterations = 5,
    sample_date_ranges = [
        ("$(2001 + 2*i)-12-1", "$(2003 + 2*i)-9-1") for i in 0:5
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
