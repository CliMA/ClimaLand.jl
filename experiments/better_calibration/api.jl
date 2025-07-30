# TODO: Make this into a package extension
# Don't really like how we are still hijacking the data sources file and reusing that...
module LandCalibration

import Dates

export CalibrateConfig

"""
    struct CalibrateConfig{SPINUP <: Dates.Period, EXTEND <: Dates.Period}
        nelements::Tuple{Int64, Int64}
        short_names::Vector{String}
        sample_date_ranges::Vector{NTuple{2, DATE}}
        spinup::SPINUP
        extend::EXTEND
        output_dir::String
        rng_seed::Int64
        minibatch_size::Int64
        n_iterations::Int64
    end

A configuration struct for keeping track of multiple fields that are of interest
to an user running calibration.
"""
struct CalibrateConfig{
    SPINUP <: Dates.Period,
    EXTEND <: Dates.Period,
    DATE <: Dates.AbstractDateTime,
}
    "The number of horizontal and vertical elements of the model. Used to
    determine the resolution of the mask"
    nelements::Tuple{Int64, Int64}

    "The short names of the observations used for calibration"
    short_names::Vector{String}

    "The amount of time to run a simulation before the first date of the
    minibatch"
    sample_date_ranges::Vector{NTuple{2, DATE}}

    "The amount of time to run a simulation before the first date of the
    minibatch"
    spinup::SPINUP

    "The amount of time to run a simulation after the last date of the
    minibatch"
    extend::EXTEND

    "The directory to store the iterations and members of the calibration."
    output_dir::String

    "An integer value for ensuring calibration are the same between multiple
    calibration runs with the same settings"
    rng_seed::Int64

    "The size of the minibatch for each iteration"
    minibatch_size::Int64

    "The number of iterations to run the calibration for"
    n_iterations::Int64
end

"""
    CalibrateConfig(;
        short_names,
        sample_date_ranges,
        extend,
        minibatch_size,
        n_iterations,
        nelements = (101, 15),
        spinup = Dates.Month(3),
        output_dir = "experiments/better_calibration/land_model",
        rng_seed = 42,
    )

Initialize a `CalibrateConfig` which are of interest to someone running calibration or
contains values that are needed in multiple places for calibration.
"""
function CalibrateConfig(;
    short_names,
    sample_date_ranges,
    extend,
    minibatch_size,
    n_iterations,
    nelements = (101, 15),
    spinup = Dates.Month(3),
    output_dir = "experiments/better_calibration/land_model",
    rng_seed = 42,
)
    isempty(short_names) || error("Cannot run calibration with no short names")
    isempty(sample_date_ranges) || error(
        "Cannot run calibration with no dates specified for the date ranges",
    )

    sample_date_ranges = [
        (Dates.DateTime(date_pair[1]), Dates.DateTime(date_pair[2])) for
        date_pair in sample_date_ranges
    ]

    num_samples = length(sample_date_ranges)
    minibatch_size > num_samples && error(
        "The minibatch size is $minibatch_size, but the number of samples is $num_samples",
    )

    num_samples % minibatch_size == 0 || @warn(
        "Number of samples is not divisible by the minibatch size; may be missing samples when running the calibration"
    )

    return CalibrateConfig(;
        nelements = nelements,
        short_names = short_names,
        sample_date_ranges = sample_date_ranges,
        spinup = spinup,
        extend = extend,
        output_dir = output_dir,
        rng_seed = rng_seed,
        minibatch_size = minibatch_size,
        n_iterations = n_iterations,
    )

end

end
