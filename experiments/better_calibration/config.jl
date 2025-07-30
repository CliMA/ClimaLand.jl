include(
    joinpath(pkgdir(ClimaLand), "experiments/better_calibration/api.jl")
)

const CALIBRATE_CONFIG = CalibrateConfig(;
    short_names = ["swu"],
    sample_date_ranges = [("2007-12-1", "2008-9-1")],
    extend = Dates.Month(3),
    minibatch_size = 1,
    n_iterations = 10,
    nelements = (101, 15),
    spinup = Dates.Month(3),
    output_dir = "experiments/better_calibration/test",
    rng_seed = 42,
)
