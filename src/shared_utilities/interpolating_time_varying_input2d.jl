import ..DataHandling
import ..DataHandling:
    DataHandler, previous_time, next_time, regridded_snapshot, available_times

"""
    InterpolatingTimeVaryingInput2D

The constructor for InterpolatingTimeVaryingInput2D is not supposed to be used directly, unless you
know what you are doing. The constructor does not perform any check and does not take care of
GPU compatibility. It is responsibility of the user-facing constructor TimeVaryingInput() to do so.
"""
struct InterpolatingTimeVaryingInput2D{
    DH <: DataHandler,
    M <: AbstractInterpolationMethod,
    CC <: ClimaComms.AbstractCommsContext,
    R <: Tuple,
} <: AbstractTimeVaryingInput
    """Object that has all the information on how to deal with files, data, and so on.
       Having to deal with files, it lives on the CPU."""
    data_handler::DH

    """Interpolation method"""
    method::M

    """ClimaComms context"""
    context::CC

    """Range of times over which the interpolator is defined. range is always defined on the
    CPU. Used by the in() function."""
    range::R
end

function TimeVaryingInput(
    data_handler::DataHandler;
    method = LinearInterpolation(),
    context = ClimaComms.context(),
)
    available_times = DataHandling.available_times(data_handler)
    isempty(available_times) &&
        error("DataHandler does not contain temporal data")
    issorted(available_times) || error("Can only interpolate with sorted times")
    range = (available_times[begin], available_times[end])
    return InterpolatingTimeVaryingInput2D(data_handler, method, context, range)
end

function evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput2D,
    time,
    args...;
    kwargs...,
)
    time in itp || error("TimeVaryingInput does not cover time $time")
    evaluate!(dest, itp, time, itp.method)
    return nothing
end

function evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput2D,
    time,
    ::NearestNeighbor,
    args...;
    kwargs...,
)
    t0, t1 =
        previous_time(itp.data_handler, time), next_time(itp.data_handler, time)

    # The closest regridded_snapshot could be either the previous or the next one
    if (time - t0) <= (t1 - time)
        dest .= regridded_snapshot(itp.data_handler, t0)
    else
        dest .= regridded_snapshot(itp.data_handler, t1)
    end
    return nothing
end

function evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput2D,
    time,
    ::LinearInterpolation,
    args...;
    kwargs...,
)
    # Linear interpolation is:
    # y = y0 + (y1 - y0) * (time - t0) / (t1 - t0)
    #
    # Define coeff = (time - t0) / (t1 - t0)
    #
    # y = (1 - coeff) * y0 + coeff * y1

    t0, t1 =
        previous_time(itp.data_handler, time), next_time(itp.data_handler, time)
    coeff = (time - t0) / (t1 - t0)
    dest .=
        (1 - coeff) .* regridded_snapshot(itp.data_handler, t0) .+
        coeff .* regridded_snapshot(itp.data_handler, t1)
    return nothing
end
