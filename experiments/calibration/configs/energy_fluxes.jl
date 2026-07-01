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
const NOISE_SCALARS = Dict("lwu" => 25.0, "shf" => 25.0, "lhf" => 25.0, "swu" => 25.0)

"""
    get_calibration_prior()

Return the combined prior distribution for the calibration parameters.

Calibrates one parameter: soil emissivity. True value is near 0.96.
"""
function get_calibration_prior()
    priors = [
        EKP.constrained_gaussian("pmodel_cstar", 0.41, 0.1, 0.2, 0.7),
        EKP.constrained_gaussian("pmodel_β_c3", 146.0, 40.0, 20.0, 400.0),
        EKP.constrained_gaussian("pmodel_β_c4", 16.222, 5.0, 1.0, 100.0),
        EKP.constrained_gaussian("pmodel_α", 0.933, 0.02, 0.85, 0.999),
        EKP.constrained_gaussian("moisture_stress_c", 0.59, 0.15, 0.05, 1.0),
        EKP.constrained_gaussian("alpha_0", 0.59, 0.15, 0.05, 1.0),
        EKP.constrained_gaussian("delta_alpha", 0.4, 0.15, 0.05, 1.0),
        EKP.constrained_gaussian("k", 1.96, 0.15, 0.05, 10.0),
        EKP.constrained_gaussian("beta", 0.97, 0.15, 0.05, 1.0),
        EKP.constrained_gaussian("leaf_Cd", 0.07, 0.03, 0.001, 0.5),
        EKP.constrained_gaussian("K_lw", 0.92, 0.2, 0.0, 2.0),
    ]
    return EKP.combine_distributions(priors)
end

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["lwu", "shf", "lhf", "swu"],
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
