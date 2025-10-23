import Dates
import ClimaLand

"""
    struct CalibrateConfig{SPINUP <: Dates.Period, EXTEND <: Dates.Period}
        short_names::Vector{String}
        minibatch_size::Int64
        n_iterations::Int64
        sample_date_ranges::Vector{NTuple{2, DATE}}
        extend::EXTEND
        spinup::SPINUP
        nelements::Tuple{Int64, Int64}
        output_dir::String
        rng_seed::Int64
    end

A configuration struct for keeping track of multiple fields that are of interest
to a user running calibration, or that are needed in multiple places (e.g., for
ensemble members and generating observations).
"""
struct CalibrateConfig{SPINUP <: Dates.Period, EXTEND <: Dates.Period, MODEL}
    "The short names of the observations used for calibration. The short names
    should match the same names used for the diagnostics."
    short_names::Vector{String}

    "The size of the minibatch for each iteration"
    minibatch_size::Int64

    "The number of iterations to run the calibration for"
    n_iterations::Int64

    "The date ranges of the samples for calibration and used to determine the
    start and end dates of a simulation for each iteration of calibration"
    sample_date_ranges::Vector{NTuple{2, Dates.DateTime}}

    "The amount of time to run a simulation after the last date of the
    minibatch"
    extend::EXTEND

    "The amount of time to run a simulation before the first date of the
    minibatch"
    spinup::SPINUP

    "The number of horizontal and vertical elements of the model. Used for the
    simulation and determining the ocean mask"
    nelements::Tuple{Int64, Int64}

    "The directory to store the iterations and members of the calibration."
    output_dir::String

    "An integer value for ensuring calibrations are the same between multiple
    calibrations with the same settings"
    rng_seed::Int64

    "File path to JLD2 file that stores a vector of EKP.Observation"
    obs_vec_filepath::String

    """Type of model to use"""
    model_type::MODEL
end

"""
    CalibrateConfig(;
        short_names,
        sample_date_ranges,
        extend,
        spinup = Dates.Month(3),
        minibatch_size,
        n_iterations,
        nelements = (101, 15),
        output_dir = "experiments/calibration/land_model",
        rng_seed = 42,
        model_type = ClimaLand.LandModel,
    )

Initializes a CalibrateConfig, which is of interest to a user running
calibration or contains values needed in multiple places during calibration.

Keyword arguments
=====================

- `short_names`: Short names of the observations. The currently supported short
  names are `lhf`, `shf`, `lwu`, and `swu`.

- `minibatch_size`: The size of the minibatch for each iteration.

- `n_iterations`: The number of iterations to run the calibration for.

- `sample_date_ranges`: The date ranges for each sample. The dates should be the
  same as found in the time series data of the observations. Since the land
  calibration calibrates using seasonal averages, the times passed must be the
  first day of December, March, June, or September. The seasons are December to
  February (DJF), March to May (MAM), June to August (JJA), and September to
  November (SON). In addition, the start and end dates of the simulation is
  automatically determined from `sample_date_ranges`.

- `extend`: The amount of time to run the simulation after the end date
  determined by `sample_date_ranges`. For seasonal averages, `extend` should be
  `Dates.Month(3)` and for monthly averages, `extend` should be
  `Dates.Month(1)`.

- `spinup`: The amount of time to run the simulation before the start date
  determined by `sample_date_ranges`.

- `nelements`: The resolution of the model. This is also used to determine the
  mask of the observations.

- `output_dir`: The location to save the calibration at.

- `rng_seed`: An integer to ensure that calibration runs with the same settings
  are the same.

- `model_type`: A type indicating which model to use. The only supported
  models are "ClimaLand.LandModel" and "ClimaLand.Bucket.BucketModel".
"""
function CalibrateConfig(;
    short_names,
    minibatch_size,
    n_iterations,
    sample_date_ranges,
    extend,
    spinup = Dates.Month(3),
    nelements = (101, 15),
    output_dir = "experiments/calibration/land_model",
    rng_seed = 42,
    obs_vec_filepath = "experiments/calibration/land_observation_vector.jld2",
    model_type = ClimaLand.LandModel,
)
    isempty(short_names) && error("Cannot run calibration with no short names")
    isempty(sample_date_ranges) &&
        error("Cannot run calibration with no date ranges for the samples")

    sample_date_ranges = [
        (Dates.DateTime(date_pair[1]), Dates.DateTime(date_pair[2])) for
        date_pair in sample_date_ranges
    ]

    for (start_date, stop_date) in sample_date_ranges
        start_date <= stop_date || error(
            "The start date ($start_date) should be before the stop date ($stop_date)",
        )
    end
    issorted(sample_date_ranges) ||
        error("The samples in $sample_date_ranges should be sorted")

    minibatch_size > 0 ||
        error("The minibatch size ($minibatch_size) should be positive")
    n_iterations > 0 ||
        error("The number of iterations ($n_iterations) should be positive")

    num_samples = length(sample_date_ranges)
    minibatch_size > num_samples && error(
        "The minibatch size is $minibatch_size, but the number of samples is $num_samples",
    )

    remaining = num_samples % minibatch_size
    remaining == 0 || @warn(
        "Number of samples is not divisible by the minibatch size; the last $remaining samples may be missing when running the calibration"
    )

    endswith(obs_vec_filepath, ".jld2") ||
        error("The file $(basename(obs_vec_filepath)) is not a JLD2 file")

    return CalibrateConfig(
        short_names,
        minibatch_size,
        n_iterations,
        sample_date_ranges,
        extend,
        spinup,
        nelements,
        output_dir,
        rng_seed,
        obs_vec_filepath,
        model_type,
    )

end
