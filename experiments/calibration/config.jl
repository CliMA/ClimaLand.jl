include(joinpath(pkgdir(ClimaLand), "experiments/calibration/api.jl"))

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["swu"],
    minibatch_size = 1,
    n_iterations = 10,
    sample_date_ranges = [("2007-12-1", "2008-9-1")],
    extend = Dates.Month(3),
    spinup = Dates.Month(3),
    nelements = (101, 15),
    output_dir = "experiments/calibration/land_model",
    rng_seed = 42,
)
