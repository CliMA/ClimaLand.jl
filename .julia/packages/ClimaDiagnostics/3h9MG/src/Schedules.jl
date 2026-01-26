"""
The `Schedules` module also contain a collection of predefined schedule
functions to implement the most common behaviors (e.g., run the callback every N
steps).
"""
module Schedules

import Dates

import ..seconds_to_str_short,
    ..seconds_to_str_long,
    ..time_to_date,
    ..period_to_str_short,
    ..period_to_str_long

import SciMLBase

import ClimaUtilities.TimeManager: ITime, date, counter, period, epoch

"""
    AbstractSchedule

`AbstractSchedule`s are structs that behave like functions and are used for the purpose of
defining a schedule to be used in `ScheduledDiagnostics`. They also may contain additional
information.
"""
abstract type AbstractSchedule end

"""
    short_name(schedule)

Short of name of the given `schedule`. Typically used in names of files/datasets.
"""
function short_name end

"""
    long_name(schedule)

Long of name of the given `schedule`. Typically used in attributes.
"""
function long_name end

function Base.show(io::IO, schedule::AbstractSchedule)
    # This function is used in names of files/datasets
    print(io, short_name(schedule))
end

function Base.:(==)(schedule1::T, schedule2::T) where {T <: AbstractSchedule}
    # The Schedules are the identical when the properties are the same, but for
    # Refs, we have to unpack the value because we don't want to compare
    # pointers.
    for p in propertynames(schedule1)
        if getproperty(schedule2, p) isa Base.RefValue
            getproperty(schedule1, p)[] == getproperty(schedule2, p)[] ||
                return false
        else
            getproperty(schedule1, p) == getproperty(schedule2, p) ||
                return false
        end
    end
    return true
end


"""
    DivisorSchedule

True when the iteration number is evenly divisible by a given number.

This is roughly equivalent to: "run this call back every N steps", with the difference that
no initial offset is possible.
"""
struct DivisorSchedule <: AbstractSchedule
    """Return true when the step number is divided evenly by this number (ie, step % divisor
        == 0) """
    divisor::Int
end

"""
    DivisorSchedule(integrator)

Returns true if `integrator.step` is evenly divided by the divisor.
"""
function (schedule::DivisorSchedule)(integrator)::Bool
    return rem(integrator.step, schedule.divisor) == 0
end

"""
    short_name(schedule::DivisorSchedule)

Short name of the given `schedule`. Typically used in names of files/datasets.

By default, the name of this schedule is `<divisor>it`, with `<divisor>` the value.
"""
function short_name(schedule::DivisorSchedule)
    return "$(schedule.divisor)it"
end

"""
    long_name(schedule::DivisorSchedule)

Long name of the given `schedule`. Typically used in attributes.

By default, the name of this schedule is "every <divisor> iterations" (even this
is not technically correct...).
"""
function long_name(schedule::DivisorSchedule)
    return "every $(schedule.divisor) iterations"
end

"""
    EveryStepSchedule()

Return a schedule that executes at the end of every step.
"""
function EveryStepSchedule()
    return DivisorSchedule(1)
end

"""
    EveryDtSchedule

True every time the current time is larger than the previous time this schedule was true + Dt.

Note, this function performs no checks on whether the step is aligned with `dt` or not.
"""
struct EveryDtSchedule{T} <: AbstractSchedule
    """The integrator time the last time this function returned true."""
    t_last::Base.RefValue{T}

    """The interval of time needed to elapse for the next time that this function will
    return true."""
    dt::T

    """
        EveryDtSchedule(dt; t_last = zero(dt))

    True every time the current time is larger than the last time this schedule was true + dt.

    The default value for `t_last` assumes that 0 is a relevant time in your problem (e.g.,
    the time at the beginning of the simulation). Adjust it if this is not the case.
    """
    function EveryDtSchedule(
        dt::T;
        t_last::T = zero(dt),
        ######## DEPRECATED #########
        t_start = nothing,
        ######## DEPRECATED #########
    ) where {T}
        ######## DEPRECATED #########
        if !isnothing(t_start)
            Base.depwarn(
                "`t_start` is deprecated and will be ignored",
                :EveryDtSchedule,
            )
        end
        ######## DEPRECATED #########
        new{typeof(dt)}(Ref(t_last), dt)
    end
end

"""
    EveryDtSchedule(integrator)

Returns true if `integrator.t >= last_t + dt`, where `last_t` is the last time
this function was true and `dt` is the schedule interval time.
"""
function (schedule::EveryDtSchedule)(integrator)::Bool
    next_t = schedule.t_last[] + schedule.dt
    # Dealing with floating point precision...
    if integrator.t > next_t || integrator.t ≈ next_t
        schedule.t_last[] = integrator.t
        return true
    else
        return false
    end
end

"""
    short_name(schedule::EveryDtSchedule)

Short of name of the given `schedule`. Typically used in names of files/datasets.

By default, the name of this schedule is the value converted into DDd_HHh_MMm_SSs.

Note:

This assumes that units are seconds.
"""
function short_name(schedule::EveryDtSchedule)
    return seconds_to_str_short(schedule.dt)
end

"""
    long_name(schedule::EveryDtSchedule)

Short of name of the given `schedule`. Typically used in attributes.

Note:

This assumes that units are seconds.
"""
function long_name(schedule::EveryDtSchedule)
    return seconds_to_str_long(schedule.dt)
end


"""
    EveryCalendarDtSchedule

Returns true if `dt` has passed since the last time this schedule was true. `dt`
here is a `Dates.Period` (e.g., `Dates.Month(1)`).

!!! compat "ClimaDiagnostics 0.2.4"
    This schedule was introduced in version `0.2.4`.
"""
struct EveryCalendarDtSchedule{P <: Dates.Period} <: AbstractSchedule
    """Last date this function returned true."""
    date_last::Base.RefValue{Dates.DateTime}

    """The `Dates.Period` needed to elapse for the next time that this function will return
    true."""
    dt::P

    """The `Dates.DateTime` used to convert from simulation time to date."""
    start_date::Dates.DateTime

    """
        EveryCalendarDtSchedule(dt::Dates.Period,
                                start_date::Dates.DateTime,
                                date_last::Union{Nothing, Dates.DateTime} = nothing)

    True every time the current date is larger than the last date this schedule was true + dt.

    The date is computed assuming that `integrator.t` is in seconds and using `start_date`.
    Schematically:
    ```julia
    date = start_date + Second(integrator.t)
    ```

    When `date_last` is `nothing`, `date_last` is assumed to be `start_date`.

    !!! compat "ClimaDiagnostics 0.2.4"
        This schedule was introduced in version `0.2.4`.
    """
    function EveryCalendarDtSchedule(
        dt::Union{Dates.Period, ITime};
        # TODO: When removing the deprecated arguments, remove the Nothing in this Union
        start_date::Union{Dates.Date, Dates.DateTime, ITime, Nothing} = nothing,
        date_last::Union{Dates.DateTime, ITime, Nothing} = nothing,
        ######## DEPRECATED #########
        reference_date = nothing,
        t_start = nothing,
        ######## DEPRECATED #########
    )
        ######## DEPRECATED #########
        if !isnothing(t_start)
            Base.depwarn(
                "`t_start` is deprecated and will be ignored",
                :EveryCalendarDtSchedule,
            )
        end
        if !isnothing(reference_date)
            start_date = reference_date
            Base.depwarn(
                "The keyword argument `reference_date` is deprecated. Use `start_date` instead.",
                :EveryCalendarDtSchedule,
            )
        end
        ######## DEPRECATED #########
        dt isa ITime && (dt = counter(dt) * period(dt))
        start_date isa ITime && (start_date = date(start_date))
        date_last isa ITime && (date_last = date(date_last))
        isnothing(date_last) && (date_last = start_date)
        new{typeof(dt)}(Ref(date_last), dt, start_date)
    end
end

"""
    EveryCalendarDtSchedule(dt, t::ITime)

Construct an `EveryCalendarDtSchedule` from `dt` and `t`, where `dt` is a period or an `ITime`
and `t` is an `ITime`. The start date of the schedule is `epoch(t)`, and the last date on
which the schedule returns true is `date(t)`.
"""
function EveryCalendarDtSchedule(dt, t::ITime)
    start_date = epoch(t)
    date_last = date(t)
    return EveryCalendarDtSchedule(
        dt;
        start_date = start_date,
        date_last = date_last,
    )
end

"""
    EveryCalendarDtSchedule(integrator)

Returns true if `current_date >= last_date + dt`, where `last_date` is the last time
this function was true and `dt` is the schedule interval time.

`current_date` is computed using the schedule `start_date` and the integrator time.
See constructor for more information.
"""
function (schedule::EveryCalendarDtSchedule)(integrator)::Bool
    next_date = schedule.date_last[] + schedule.dt
    start_date = schedule.start_date
    current_date = time_to_date(integrator.t, start_date)
    if current_date >= next_date
        schedule.date_last[] = current_date
        return true
    else
        return false
    end
end

"""
    short_name(schedule::EveryCalendarDtSchedule)

Short of name of the given `schedule`. Typically used in names of files/datasets.

By default, the name of this schedule is the value converted into DDd_HHh_MMm_SSs.

Note:

This assumes that units are seconds.
"""
function short_name(schedule::EveryCalendarDtSchedule)
    return period_to_str_short(schedule.dt)
end

"""
    long_name(schedule::EveryCalendarDtSchedule)

Long of name of the given `schedule`. Typically used in attributes.

This is directly the string representation of a `Dates.Period`.
"""
function long_name(schedule::EveryCalendarDtSchedule)
    return period_to_str_long(schedule.dt)
end

end
