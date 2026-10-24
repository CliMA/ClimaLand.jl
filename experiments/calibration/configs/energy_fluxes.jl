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
const NOISE_SCALARS = Dict("lwu" => 10.0, "swu" => 10.0, "shf" => 15.0, "lhf" => 15.0, "gpp" =>1.5, "snowc" => 0.001)

"""
    get_calibration_prior()

Return the combined prior distribution for the calibration parameters.

Calibrates one parameter: soil emissivity. True value is near 0.96.
"""
function get_calibration_prior()
    priors = [
        EKP.constrained_gaussian("topo_scf_z0", 0.5, 0.2, 0.01, 2.0),
        EKP.constrained_gaussian("topo_scf_mu_exp", 0.3, 0.1, 0.01, 0.9),
        EKP.constrained_gaussian("topo_scf_xi_exp", 2.0, 0.3, 0.01, 4.0),
        EKP.constrained_gaussian("weibull_moisture_stress_c", 2.0, 0.5, 0.05, 6.0),
	EKP.constrained_gaussian("weibull_moisture_stress_a", 1.0, 0.2, 0.05, 3.0),
        EKP.constrained_gaussian("stem_PAR_reflectance", 0.3, 0.1, 0.01, 0.9),
        EKP.constrained_gaussian("stem_NIR_reflectance", 0.3, 0.1, 0.01, 0.9),
        EKP.constrained_gaussian("leaf_Cd", 0.07, 0.03, 0.005, 0.5),
	EKP.constrained_gaussian("canopy_subcanopy_min_reference_height", 1.0, 0.3, 0.5, 4.0),
        EKP.constrained_gaussian("canopy_subcanopy_wind_extinction_coefficient", 0.5, 0.02, 0.0, 2.0),	
        EKP.constrained_gaussian("canopy_K_lw", 1.16, 0.3, 0.1, 3.0),
        EKP.constrained_gaussian("maximum_dry_soil_layer_depth", 0.05, 0.03, 0.0, 0.2),
        EKP.constrained_gaussian("evaporation_scalar", 0.1, 0.02, 0.03, 0.8),
    ]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["lwu", "swu", "shf", "lhf","gpp","snowc"],
    minibatch_size = 1,
    n_iterations = 5,
     sample_date_ranges = [
        ("$(2001 + 2*i)-12-1", "$(2003 + 2*i)-9-1") for i in 0:4
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
