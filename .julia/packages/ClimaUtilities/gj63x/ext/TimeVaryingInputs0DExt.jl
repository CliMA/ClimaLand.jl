module TimeVaryingInputs0DExt

import ClimaCore
import ClimaCore: ClimaComms
import ClimaCore.Fields: Adapt

import ClimaUtilities.Utils:
    searchsortednearest, linear_interpolation, wrap_time, isequispaced
import ClimaUtilities.TimeVaryingInputs:
    AbstractInterpolationMethod, AbstractTimeVaryingInput
import ClimaUtilities.TimeVaryingInputs:
    NearestNeighbor,
    LinearInterpolation,
    LinearPeriodFillingInterpolation,
    Throw,
    Flat,
    PeriodicCalendar
import ClimaUtilities.TimeVaryingInputs: extrapolation_bc

import ClimaUtilities.TimeVaryingInputs

import ClimaUtilities.TimeManager: ITime, date
using Dates

"""
    InterpolatingTimeVaryingInput0D

The constructor for InterpolatingTimeVaryingInput0D is not supposed to be used directly, unless you
know what you are doing.

`times` and `values` may have different types, but they must be the same length, and we
assume that they have been sorted to be monotonically increasing in time, without repeated
values for the same timestamp. `times` can have elements of type `ITime` or floats
"""
struct InterpolatingTimeVaryingInput0D{
    AA1 <: AbstractArray,
    AA2 <: AbstractArray,
    M <: AbstractInterpolationMethod,
    R <: Tuple,
} <: AbstractTimeVaryingInput
    # AA1 and AA2 could be different because of different FTs

    """Independent coordinate"""
    times::AA1

    """Variable"""
    vals::AA2

    """Interpolation method"""
    method::M

    """Range of times over which the interpolator is defined. range is always defined on the
    CPU. Used by the in() function."""
    range::R
end

"""
    in(time, itp::InterpolatingTimeVaryingInput23D)

Check if the given `time` is in the range of definition for `itp`.
"""
function Base.in(time, itp::InterpolatingTimeVaryingInput0D)
    return itp.range[1] <= time <= itp.range[2]
end

function TimeVaryingInputs.evaluate!(
    destination,
    itp::InterpolatingTimeVaryingInput0D,
    time,
    args...;
    kwargs...,
)
    if extrapolation_bc(itp.method) isa Throw
        time in itp || error("TimeVaryingInput does not cover time $time")
    end
    scalar_dest = [zero(eltype(destination))]

    TimeVaryingInputs.evaluate!(scalar_dest, itp, time, itp.method)
    fill!(destination, scalar_dest[])

    return nothing
end

function TimeVaryingInputs.TimeVaryingInput(
    times::AbstractArray,
    vals::AbstractArray;
    context = nothing,
    method::AbstractInterpolationMethod = LinearInterpolation(),
)
    ########### DEPRECATED ###############
    if !isnothing(context)
        Base.depwarn(
            "The keyword argument `context` is no longer required for TimeVaryingInputs. It will be removed.",
            :TimeVaryingInput,
        )
    end
    ########### DEPRECATED ###############

    issorted(times) || error("Can only interpolate with sorted times")
    length(times) == length(vals) ||
        error("times and vals have different lengths")

    if method isa LinearPeriodFillingInterpolation
        error(
            "LinearPeriodFillingInterpolation is not supported when the input data is 1D",
        )
    end

    if extrapolation_bc(method) isa PeriodicCalendar
        if extrapolation_bc(method) isa PeriodicCalendar{Nothing}
            if !isequispaced(
                eltype(times) <: ITime ? [float.(promote(times...))...] : times,
            )
                error(
                    "PeriodicCalendar() boundary condition cannot be used because data is defined at non uniform intervals of time",
                )
            end
        else
            # We have the period in PeriodicCalendar
            error(
                "PeriodicCalendar(period) is not supported when the input data is 1D",
            )
        end
    end

    range = (times[begin], times[end])
    return InterpolatingTimeVaryingInput0D(
        copy(times),
        copy(vals),
        method,
        range,
    )
end

function _evalulate_flat!(dest, itp::InterpolatingTimeVaryingInput0D, time)
    t_init, t_end = itp.range
    if time >= t_end
        dest .= itp.vals[end]
    else
        time <= t_init
        dest .= itp.vals[begin]
    end
    return nothing
end

function _time_range_target_time_dt(itp::InterpolatingTimeVaryingInput0D, time)
    t_init, t_end = itp.range
    dt = itp.times[begin + 1] - itp.times[begin]
    time = wrap_time(time, t_init, t_end + dt)
    return t_init, t_end, time, dt
end

function TimeVaryingInputs.evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput0D,
    time,
    ::NearestNeighbor,
)
    # Nearest neighbor interpolation: just pick the values corresponding to the entry in
    # itp.times that is closest to the given time.
    if extrapolation_bc(itp.method) isa Flat
        _evalulate_flat!(dest, itp, time)
        return nothing
    elseif extrapolation_bc(itp.method) isa PeriodicCalendar{Nothing}
        t_init, t_end, time, dt = _time_range_target_time_dt(itp, time)
        # Now time is between t_init and t_end + dt. We are doing nearest neighbor
        # interpolation here, and when time >= t_end + 0.5dt we need to use t_init instead
        # of t_end as neighbor.
        time >= t_end + 0.5dt && (time = t_init)
    end

    index = searchsortednearest(itp.times, time)

    dest .= itp.vals[index]

    return nothing
end

"""
    evaluate!(
        dest,
        itp::InterpolatingTimeVaryingInput0D,
        time,
        ::LinearInterpolation,
        )

Write to `dest` the result of a linear interpolation of `itp` on the given `time`.
"""
function TimeVaryingInputs.evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput0D,
    time,
    ::LinearInterpolation,
)
    if extrapolation_bc(itp.method) isa Flat
        _evalulate_flat!(dest, itp, time)
        return nothing
    elseif extrapolation_bc(itp.method) isa PeriodicCalendar
        t_init, t_end, time, dt = _time_range_target_time_dt(itp, time)
        # We have to handle separately the edge case where the desired time is past t_end.
        # In this case, we know that t_end <= time <= t_end + dt and we have to do linear
        # interpolation between t_init and t_end. In this case, y0 = vals[end], y1 =
        # vals[begin], t1 - t0 = dt, and time - t0 = time - t_end
        if time > t_end
            if time isa ITime
                FT = eltype(itp.vals)
                @. dest =
                    itp.vals[end] + FT(
                        (itp.vals[begin] - itp.vals[end]) / float(dt) *
                        float((time - t_end)),
                    )
            else
                @. dest =
                    itp.vals[end] +
                    (itp.vals[begin] - itp.vals[end]) / dt * (time - t_end)
            end
            return nothing
        end
    end

    indep_vars = itp.times
    indep_value = time
    dep_vars = itp.vals
    dest .= linear_interpolation(indep_vars, dep_vars, indep_value)
    return nothing
end

function TimeVaryingInputs.evaluate!(
    destination,
    itp::InterpolatingTimeVaryingInput0D{<:AbstractArray{<:Number}},
    time::ITime,
    args...;
    kwargs...,
)
    return TimeVaryingInputs.evaluate!(
        destination,
        itp,
        eltype(itp.range)(float(time)),
        args...;
        kwargs...,
    )
end

function TimeVaryingInputs.evaluate!(
    destination,
    itp::InterpolatingTimeVaryingInput0D{<:AbstractArray{<:ITime}},
    eval_date::DateTime,
    args...;
    kwargs...,
)
    t0_date = date(itp.range[1])
    diff_ms = eval_date - t0_date
    converted_time =
        ITime(diff_ms.value; period = typeof(diff_ms)(1), epoch = t0_date)
    return TimeVaryingInputs.evaluate!(
        destination,
        itp,
        converted_time,
        args...;
        kwargs...,
    )

end

function TimeVaryingInputs.evaluate!(
    destination,
    itp::InterpolatingTimeVaryingInput0D{<:AbstractArray{<:ITime}},
    time::Number,
    args...;
    kwargs...,
)
    converted_time = first(promote(ITime(time), itp.range[1]))
    return TimeVaryingInputs.evaluate!(
        destination,
        itp,
        converted_time,
        args...;
        kwargs...,
    )

end

TimeVaryingInputs.evaluate!(
    destination,
    itp::InterpolatingTimeVaryingInput0D{<:AbstractArray{<:Number}},
    time::DateTime,
    args...;
    kwargs...,
) = error(
    "Cannot evaluate InterpolatingTimeVaryingInput0D with times as numbers and inputs as DateTime",
)

end
