module TimeVaryingInputsExt

import Dates

import ClimaCore
import ClimaCore: ClimaComms
import ClimaCore: DeviceSideContext
import ClimaCore.Fields: Adapt

import ClimaUtilities.Utils:
    searchsortednearest,
    linear_interpolation,
    isequispaced,
    wrap_time,
    bounding_dates,
    beginningofperiod,
    endofperiod,
    period_to_seconds_float,
    unique_periods

import ClimaUtilities.TimeVaryingInputs
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

import ClimaUtilities.DataHandling
import ClimaUtilities.DataHandling:
    DataHandler,
    regridded_snapshot!,
    available_times,
    available_dates,
    date_to_time,
    time_to_date,
    previous_date,
    next_date


import ClimaUtilities.TimeManager: ITime, date

# Ideally, we should be able to split off the analytic part in a different
# extension, but precompilation stops working when we do so

struct AnalyticTimeVaryingInput{F <: Function} <:
       TimeVaryingInputs.AbstractTimeVaryingInput
    # func here has to be GPU-compatible (e.g., splines are not) and reasonably fast (e.g.,
    # no large allocations)
    func::F
end

# _kwargs... is needed to seamlessly support the other TimeVaryingInputs.
function TimeVaryingInputs.TimeVaryingInput(
    input::Function;
    method = nothing,
    _kwargs...,
)
    isnothing(method) ||
        @warn "Interpolation method is ignored for analytical functions"
    return AnalyticTimeVaryingInput(input)
end

function TimeVaryingInputs.evaluate!(
    dest,
    input::AnalyticTimeVaryingInput,
    time,
    args...;
    kwargs...,
)
    dest .= input.func(time, args...; kwargs...)
    return nothing
end

"""
    InterpolatingTimeVaryingInput23D

The constructor for InterpolatingTimeVaryingInput23D is not supposed to be used directly, unless you
know what you are doing. The constructor does not perform any check and does not take care of
GPU compatibility. It is responsibility of the user-facing constructor TimeVaryingInput() to do so.
"""
struct InterpolatingTimeVaryingInput23D{
    DH,
    M <: AbstractInterpolationMethod,
    CC <: ClimaComms.AbstractCommsContext,
    R <: Tuple,
    RR,
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

    """Preallocated memory for storing regridded fields"""
    preallocated_regridded_fields::RR
end

"""
    in(time, itp::InterpolatingTimeVaryingInput23D)

Check if the given `time` is in the range of definition for `itp`.
"""
function Base.in(time, itp::InterpolatingTimeVaryingInput23D)
    return itp.data_handler.available_dates[begin] <=
           time <=
           itp.data_handler.available_dates[end]
end

function Base.in(time::Number, itp::InterpolatingTimeVaryingInput23D)
    return Base.in(time_to_date(itp.data_handler, time), itp)
end


function TimeVaryingInputs.TimeVaryingInput(
    data_handler;
    method = LinearInterpolation(),
    context = ClimaComms.context(),
)
    available_times = DataHandling.available_times(data_handler)
    isempty(available_times) &&
        error("DataHandler does not contain temporal data")
    issorted(available_times) || error("Can only interpolate with sorted times")
    range = (available_times[begin], available_times[end])

    # TODO: Generalize the number of _regridded_fields depending on the interpolation
    # stencil. At the moment, we use four for LinearPeriodFilling
    _num_fields = method isa LinearPeriodFillingInterpolation ? 4 : 2
    preallocated_regridded_fields =
        ntuple(_ -> zeros(data_handler.target_space), _num_fields)

    return InterpolatingTimeVaryingInput23D(
        data_handler,
        method,
        context,
        range,
        preallocated_regridded_fields,
    )
end

function TimeVaryingInputs.TimeVaryingInput(
    file_paths,
    varnames,
    target_space;
    method = LinearInterpolation(),
    start_date::Union{Dates.DateTime, Dates.Date} = Dates.DateTime(1979, 1, 1),
    regridder_type = nothing,
    regridder_kwargs = (),
    file_reader_kwargs = (),
    compose_function = identity,
    ########### DEPRECATED ###############
    reference_date = nothing,
    t_start = nothing,
    ########### DEPRECATED ###############
)
    ########### DEPRECATED ###############
    if !isnothing(reference_date)
        start_date = reference_date
        Base.depwarn(
            "The keyword argument `reference_date` is deprecated. Use `start_date` instead.",
            :TimeVaryingInput,
        )
    end
    if !isnothing(t_start)
        Base.depwarn("`t_start` was removed will be ignored", :TimeVaryingInput)
    end
    ########### DEPRECATED ###############

    data_handler = DataHandling.DataHandler(
        file_paths,
        varnames,
        target_space;
        start_date,
        regridder_type,
        regridder_kwargs,
        file_reader_kwargs,
        compose_function,
    )
    if extrapolation_bc(method) isa PeriodicCalendar{Nothing} &&
       !isequispaced(DataHandling.available_times(data_handler))
        error(
            "PeriodicCalendar() boundary condition cannot be used because data is defined at non uniform intervals of time",
        )
    end
    context = ClimaComms.context(target_space)
    return TimeVaryingInputs.TimeVaryingInput(data_handler; method, context)
end

function TimeVaryingInputs.evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput23D,
    time,
    args...;
    kwargs...,
)
    if extrapolation_bc(itp.method) isa Throw
        time in itp || error("TimeVaryingInput does not cover time $time")
    end
    if extrapolation_bc(itp.method) isa Flat
        date_init, date_end = itp.data_handler.available_dates[begin],
        itp.data_handler.available_dates[end]
        if time >= date_end
            regridded_snapshot!(dest, itp.data_handler, date_end)
        else
            time <= date_init
            regridded_snapshot!(dest, itp.data_handler, date_init)
        end
    else
        TimeVaryingInputs.evaluate!(dest, itp, time, itp.method)
    end
    return nothing
end

function TimeVaryingInputs.evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput23D,
    time::Number,
    args...;
    kwargs...,
)
    TimeVaryingInputs.evaluate!(
        dest,
        itp,
        Dates.Millisecond(round(1_000 * time)) + itp.data_handler.start_date,
        args...,
        kwargs...,
    )
    return nothing
end

function TimeVaryingInputs.evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput23D,
    time::ITime,
    args...;
    kwargs...,
)
    TimeVaryingInputs.evaluate!(dest, itp, date(time), args..., kwargs...)
    return nothing
end

function _time_range_dt_dt_e(itp::InterpolatingTimeVaryingInput23D)
    return _time_range_dt_dt_e(itp, extrapolation_bc(itp.method))
end

function _time_range_dt_dt_e(
    itp::InterpolatingTimeVaryingInput23D,
    extrapolation_bc::PeriodicCalendar{Nothing},
)
    dt = DataHandling.dt(itp.data_handler)
    return itp.data_handler.available_dates[begin],
    itp.data_handler.available_dates[end],
    Dates.Millisecond(round(1_000 * dt)),
    Dates.Millisecond(round((1_000 * dt) / 2))
end

function _time_range_dt_dt_e(
    itp::InterpolatingTimeVaryingInput23D,
    extrapolation_bc::PeriodicCalendar,
)
    period, repeat_date = extrapolation_bc.period, extrapolation_bc.repeat_date

    date_init, date_end =
        bounding_dates(available_dates(itp.data_handler), repeat_date, period)
    # Suppose date_init date_end are 15/01/23 and 14/12/23
    # dt_e is endofperiod(14/12/23) - 14/12/23
    # dt is 15/01/23 + period - 14/12/23
    # if period = 1 Year, dt_e = 17 days (in seconds)

    t_init, t_end = date_to_time(itp.data_handler, date_init),
    date_to_time(itp.data_handler, date_end)
    # We have to add 1 Second because endofperiod(date_end, period) returns the very last
    # second before the next period
    dt_e = (endofperiod(date_end, period) + Dates.Second(1) - date_end)
    dt = (date_init + period - date_end)
    return date_init, date_end, dt, dt_e
end

"""
    _interpolation_times_periodic_calendar(time, itp::InterpolatingTimeVaryingInput23D)

Return time, t_init, t_end, dt, dt_e.

Implementation details
======================

Okay, how are we implementing PeriodicCalendar?

There are two modes, one with provided `period` and `repeat_date`, and the other without.
When it comes to implementation, we reduce the first case to the second one. So, let's start
by looking at the second case, then, we will look at how we reduce it to the first one.

In the second case, we have `t_init`, `t_end`, and a `dt`. `t_init`, `t_end` define the earliest and
latest data we are going to use and are in units of simulation time. `dt` is so that `t_init =
t_end + dt`. We are also given a `dt_e` so that, for interpolation purposes, we attribute
points that are within `t_end + dt_e` to `t_end`, and points that are beyond that to `t_init`. For
equispaced timeseries, `dt_e = 0.5dt`.

Once we have all of this, we can wrap the given time to be within `t_init` and `t_end + dt`. For
all the cases where the wrapped time is between `t_init` and `t_end`, the function can use the
standard interpolation scheme, so, the only case we have to worry about is when the wrapped
time is between `t_end` and `t_end + dt`. We handle this case manually by working explicitly
with `dt_e`.

Now, let us reduce the case where we are given dates and a period.

Let us look at an example, suppose we have data defined at these dates

16/12/22, 15/01/23 ...,  14/12/23, 13/12/24

and we want to repeat the year 2023. period will be `Dates.Year` and `repeat_date` will be
01/01/2023 (or any date in the year 2023)

First, we identify the bounding dates that correspond to the given period that has to be
repeated. We can assume that dates are sorted. In this case, this will be 15/01/23 and
14/12/23. Then, we translate this into simulation time and compute `dt` and `dt_e`. That's it!
"""
function _interpolation_times_periodic_calendar(
    time,
    itp::InterpolatingTimeVaryingInput23D,
)
    t_init, t_end, dt, dt_e = _time_range_dt_dt_e(itp)
    time = wrap_time(time, t_init, t_end + dt)
    return time, t_init, t_end, dt, dt_e
end

function TimeVaryingInputs.evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput23D,
    time,
    ::NearestNeighbor,
    args...;
    kwargs...,
)
    if extrapolation_bc(itp.method) isa PeriodicCalendar
        time, t_init, t_end, _, dt_e =
            _interpolation_times_periodic_calendar(time, itp)

        # Now time is between t_init and t_end + dt. We are doing nearest neighbor
        # interpolation here, and when time >= t_end + dt_e we need to use t_init instead of
        # t_end as neighbor.

        # TODO: It would be nice to handle this edge case directly instead of copying the
        # code
        if (time - t_end) <= dt_e
            regridded_snapshot!(dest, itp.data_handler, t_end)
        else
            regridded_snapshot!(dest, itp.data_handler, t_init)
        end
        return nothing
    end
    date0, date1 =
        previous_date(itp.data_handler, time), next_date(itp.data_handler, time)

    # The closest regridded_snapshot could be either the previous or the next one
    if (time - date0) <= (date1 - time)
        regridded_snapshot!(dest, itp.data_handler, date0)
    else
        regridded_snapshot!(dest, itp.data_handler, date1)
    end
    return nothing
end

function TimeVaryingInputs.evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput23D,
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

    field_t0, field_t1 = itp.preallocated_regridded_fields[begin:(begin + 1)]

    if extrapolation_bc(itp.method) isa PeriodicCalendar
        time, date_init, date_end, dt, _ =
            _interpolation_times_periodic_calendar(time, itp)

        # We have to handle separately the edge case where the desired time is past t_end.
        # In this case, we know that t_end <= time <= t_end + dt and we have to do linear
        # interpolation between t_init and t_end. In this case, y0 = regridded_field(t_end),
        # y1 = regridded_field(t_init), t1 - t0 = dt, and time - t0 = time - t_end

        # TODO: It would be nice to handle this edge case directly instead of copying the
        # code
        if time > date_end
            regridded_snapshot!(field_t0, itp.data_handler, date_end)
            regridded_snapshot!(field_t1, itp.data_handler, date_init)
            coeff = (time - date_end) / dt
            dest .= (1 - coeff) .* field_t0 .+ coeff .* field_t1
            return nothing
        end
    end

    # We have to consider the edge case where time is precisely the last available_time.
    # This is relevant also because it can be triggered by LinearPeriodFilling
    if time in DataHandling.available_dates(itp.data_handler)
        regridded_snapshot!(dest, itp.data_handler, time)
    else
        date0, date1 = previous_date(itp.data_handler, time),
        next_date(itp.data_handler, time)
        coeff = (time - date0) / (date1 - date0)

        regridded_snapshot!(field_t0, itp.data_handler, date0)
        regridded_snapshot!(field_t1, itp.data_handler, date1)

        dest .= (1 - coeff) .* field_t0 .+ coeff .* field_t1
    end
    return nothing
end

include("time_varying_inputs_linearperiodfilling.jl")

"""
    close(time_varying_input::TimeVaryingInputs.AbstractTimeVaryingInput)

Close files associated to the `time_varying_input`.
"""
function Base.close(
    time_varying_input::TimeVaryingInputs.AbstractTimeVaryingInput,
)
    return nothing
end

"""
    close(time_varying_input::InterpolatingTimeVaryingInput23D)

Close files associated to the `time_varying_input`.
"""
function Base.close(time_varying_input::InterpolatingTimeVaryingInput23D)
    Base.close(time_varying_input.data_handler)
    return nothing
end

end
