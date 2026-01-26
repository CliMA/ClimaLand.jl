export ITime, counter, period, epoch, date, seconds

"""
    ITime ("Integer Time")

`ITime` is an integer (quantized) time.

`ITime` can be thought of as counting clock cycles (`counter`), with each tick
having a fixed duration (`period`).

Another way to think about this is that this is time with units.

This type is currently using Dates, but one of the design goals is to try to be
as agnostic as possible with respect to this so that in the future in will be
possible to use a different calendar.

When using Dates, the minimum unit of time that can be represented is 1
nanosecond. The maximum unit of time is determined by the maximum integer number
that can be represented.

Overflow occurs at `68 year * (1 Second / dt)` for Int32 and `300 gigayear * (1
Second / dt)` for Int64.

# Fields
- `counter::INT`: The number of clock cycles.
- `period::DT`: The duration of each cycle.
- `epoch::EPOCH`: An optional start date.
"""
struct ITime{INT <: Integer, DT, EPOCH <: Union{Nothing, Dates.DateTime}}
    counter::INT
    period::DT
    epoch::EPOCH

    function ITime(counter, period, epoch)
        return new{typeof(counter), typeof(period), typeof(epoch)}(
            counter,
            period,
            epoch,
        )
    end
end

"""
    ITime(counter::Integer; period::Union{Dates.FixedPeriod, Nothing}, epoch = nothing)

Construct an `ITime` from a counter, a period, and an optional start date.

If the `epoch` is provided as a `Date`, it is converted to a `DateTime`.

If period is not provided, the function assumes that `counter` describes a time in seconds,
and attempts to find a `Dates.FixedPeriod` such that `counter` seconds can be
represented as an integer multiple of that period.

If `counter` is zero, it defaults to a period of 1 second.
"""
function ITime(
    counter::Integer;
    period::Union{Dates.FixedPeriod, Nothing} = nothing,
    epoch = nothing,
)
    # Convert epoch to DateTime if it is not nothing (from, e.g., Date)
    isnothing(epoch) || (epoch = Dates.DateTime(epoch))
    if !isnothing(period)
        return ITime(counter, period, epoch)
    else
        counter == 0 && return ITime(counter, Dates.Second(1), epoch)
        potential_periods = [
            Dates.Week,
            Dates.Day,
            Dates.Hour,
            Dates.Minute,
            Dates.Second,
            Dates.Millisecond,
            Dates.Microsecond,
            Dates.Nanosecond,
        ]
        for period_i in potential_periods
            period_ns = Dates.tons(period_i(1))
            t_int = 1_000_000_000 / period_ns * counter
            if isinteger(t_int)
                return ITime(
                    typeof(counter)(div(1_000_000_000 * counter, period_ns)),
                    period_i(1),
                    epoch,
                )
            end
        end
    end
end

"""
    seconds(t::ITime)

Return the time represented by `t` in seconds, as a floating-point number.
"""
function seconds(t::ITime)
    return float(t)
end

# Accessors
"""
    counter(t::ITime)

Return the counter of the `ITime` `t`.
"""
function counter(t::ITime)
    return t.counter
end

"""
    period(t::ITime)

Return the period of the `ITime` `t`.
"""
function period(t::ITime)
    return t.period
end

"""
    epoch(t::ITime)

Return the start date of the `ITime` `t`.
"""
function epoch(t::ITime)
    return t.epoch
end

function date(t::ITime{<:Integer, <:Dates.FixedPeriod, Nothing})
    error("Time does not have epoch information")
end

"""
    date(t::ITime)

Return the date associated with `t`. If the time is fractional, round it to
millisecond.

For this to work, `t` has to have a `epoch`.
"""
function date(t::ITime)
    current_date = epoch(t)
    elapsed_period = counter(t) * period(t)

    # We only check if `elapsed_period` overflows.
    # It is reasonable to assume that the date will not overflow since the
    # maximum value for the date is 292277025-08-17T07:12:55.807
    # which is computed by
    # `DateTime(Dates.UTInstant(Millisecond(typemax(Int64))))`.
    # Meanwhile, it is possible for `elapsed_period` to overflow if the period
    # is extremely small (e.g. nanosecond) and Int64 is used for the counter.
    # Although, this will not happen for the nanosecond and Int64 case for ~294
    # years
    period_is_zero = period(t) == zero(period(t))
    no_overflow = period_is_zero || div(elapsed_period, period(t)) == counter(t)
    no_overflow && return current_date + elapsed_period
    error(
        "Overflow with counter(t) * period(t) in computing the date; try to use a bigger value for the period if possible",
    )
end

"""
    DateTime(t::ITime)

Convert an `ITime` to a `DateTime`.
"""
function Dates.DateTime(t::ITime)
    return date(t)
end

"""
    ITime(t; epoch = nothing)

Construct an `ITime` from a number `t` representing a time interval.

The function attempts to find a `Dates.FixedPeriod` such that `t` can be
represented as an integer multiple of that period.

If `t` is approximately zero, it defaults to a period of 1 second.
"""
function ITime(t; epoch = nothing)
    # If it is zero, assume seconds
    isapprox(t, 0) && return ITime(0, Dates.Second(1), epoch)

    # Promote t to Float64 to avoid loss of precision
    t = Float64(t)
    periods = [
        Dates.Week,
        Dates.Day,
        Dates.Hour,
        Dates.Minute,
        Dates.Second,
        Dates.Millisecond,
        Dates.Microsecond,
        Dates.Nanosecond,
    ]
    for period in periods
        period_ns = Dates.tons(period(1))
        t_int = 1_000_000_000 / period_ns * t
        if isinteger(t_int)
            return ITime(Int(t_int), period(1), epoch)
        end
    end
    error("Cannot represent $t as integer multiple of a Dates.FixedPeriod")
end

function Base.show(io::IO, time::ITime)
    # Hack to pretty print fractional times. We cannot just use Dates to print
    # them because they cannot be nicely converted to Periods, instead of
    # reconstruct the string from the type name and the value (obtained my
    # multiplying the counter and the number of units in the period)
    value = counter(time) * float(period(time).value)
    plural_s = abs(value) != 1 ? "s" : ""
    unit = lowercase(string(nameof(typeof(period(time))))) * plural_s

    print(io, "$value $unit ")
    # Add date, if available
    if !isnothing(epoch(time)) && counter(time) isa Integer
        print(io, "($(date(time))) ")
    end
    print(io, "[counter = $(counter(time)), period = $(period(time))")
    # Add start date, if available
    if !isnothing(epoch(time))
        print(io, ", epoch = $(epoch(time))")
    end
    print(io, "]")
end

"""
    promote(ts::ITime...)

Promote a tuple of `ITime` instances to a common type.

This function determines a common `epoch` and `period` for all the input
`ITime` instances and returns a tuple of new `ITime` instances with the common
type.  It throws an error if the start dates are different.
"""
function Base.promote(ts::ITime...)
    common_epoch = find_common_epoch(ts...)

    # Determine the common period
    common_period = reduce(gcd, (period(t) for t in ts))

    # Promote each ITime instance by computing the scaling factor needed
    return map(
        t -> ITime(
            counter(t) * typeof(t.counter)(div(period(t), common_period)),
            common_period,
            common_epoch,
        ),
        ts,
    )
end

"""
    find_common_epoch(ts::ITime...)

Find a common epoch of `ts` if one exists.
"""
function find_common_epoch(ts::ITime...)
    epochs = (epoch(t) for t in ts if !isnothing(epoch(t)))
    common_epoch = reduce(_unique_epochs, epochs, init = nothing)
    return common_epoch
end

"""
    _unique_epochs(epoch1, epoch2)

Return `epoch2` if `epoch1` and `epoch2` are the same or `epoch1` is nothing,
and an error otherwise.
"""
function _unique_epochs(epoch1, epoch2)
    if isnothing(epoch1)
        return epoch2
    elseif epoch1 == epoch2
        return epoch2
    else
        return error("Cannot find common epoch")
    end
end

"""
    Base.:*(t::ITime, a::AbstractFloat)

Multiplication between an ITime and float return an `ITime` whose counter is
`round(a * counter(t))`, the same period as `t`, and the same epoch as `t` if it
exists. The float `a` can only be between 0 and 1.

This function should only be used when subdividing time is necessary (e.g. time
stepping stages in ClimaTimeSteppers). In most cases, it is preferable to
convert `t` into a float if multiplication by a float is needed.
"""
function Base.:*(t::ITime, a::AbstractFloat)
    if a > one(a) || a < zero(a)
        return error(
            "In most cases, multiplying an ITime by a float is not desirable. Cast the ITime into a float",
        )
    else
        return ITime(round(typeof(t.counter), a * t.counter), t.period, t.epoch)
    end
end

# Same as above, but needed to make multiplication commutative
function Base.:*(a::AbstractFloat, t::ITime)
    if a > one(a) || a < zero(a)
        return error(
            "In most cases, multiplying an ITime by a float is not desirable. Cast the ITime into a float",
        )
    else
        return ITime(round(typeof(t.counter), a * t.counter), t.period, t.epoch)
    end
end

Base.:(==)(t1::ITime, t2::Number) = error("Cannot compare ITime with a Number")
Base.:(==)(t1::Number, t2::ITime) = t2 == t1

"""
    Base.:(==)(t1::ITime, t2::Dates.DateTime)

Converts `t1` to a `DateTime` and checks for equality with `t2`. `t1` must have an `epoch`
`t1==t2` ⟺ `t2==t1`
"""
Base.:(==)(t1::ITime, t2::Dates.DateTime) = date(t1) == t2
Base.:(==)(t1::Dates.DateTime, t2::ITime) = t2 == t1

"""
    Base.:(:)(start::ITime, step::ITime, stop::ITime)

Range operator. `start:step:stop` constructs a range from start to stop with a
step size equal to step.
"""
function Base.:(:)(start::ITime, step::ITime, stop::ITime)
    start, step, stop = promote(start, step, stop)
    common_epoch = find_common_epoch(start, step, stop)
    return (
        ITime(count, period = start.period, epoch = common_epoch) for
        count in (start.counter):(step.counter):(stop.counter)
    )
end

"""
    Base.mod(x::ITime, y::ITime)

Return the counter of `x` modulo counter of `y` after promote `x` and `y` to the
same period and epoch.
"""
function Base.mod(x::ITime, y::ITime)
    x, y = promote(x, y)
    reminder = mod(x.counter, y.counter)
    return ITime(reminder, period = x.period, epoch = x.epoch)
end

"""
    Base.:(%)(x::ITime, y::ITime)

Return the counter of `x` modulo counter of `y` after promote `x` and `y` to the
same period and epoch.
"""
function Base.:(%)(x::ITime, y::ITime)
    return mod(x, y)
end

"""
    Base.iszero(x::ITime)

Return `true` if the counter of `x` is zero.
"""
function Base.iszero(x::ITime)
    return iszero(x.counter)
end

"""
    Base.length(x::ITime)

Return the length of an ITime which is always one.
"""
function Base.length(x::ITime)
    return 1
end

# Pay attention to the units here! zero and one are not symmetric
"""
    Base.one(t::T) where {T <: ITime}

Return the multiplicative identity for an `ITime` which is `1`.
"""
Base.one(t::T) where {T <: ITime} = 1

"""
    Base.oneunit(t::T) where {T <: ITime}

Return `ITime(1, period(t), epoch(t))`.
"""
Base.oneunit(t::T) where {T <: ITime} =
    ITime(eltype(t.counter)(1), t.period, t.epoch)

"""
    Base.zero(t::T) where {T <: ITime}

Return the additive identity element for an `ITime` which is
`ITime(0, period(t), epoch(t))`.
"""
Base.zero(t::T) where {T <: ITime} =
    ITime(eltype(t.counter)(0), t.period, t.epoch)

"""
    float(t::ITime)

Convert an `ITime` to a floating-point number representing the time in seconds.
"""
function Base.float(t::T) where {T <: ITime}
    if VERSION >= v"1.11"
        return float(Dates.seconds(t.period)) * t.counter
    else
        return Dates.tons(t.period) / 1_000_000_000 * t.counter
    end
end

(::Type{FT})(t::ITime) where {FT <: AbstractFloat} = FT(float(t))

macro itime_unary_op(op)
    return esc(
        quote
            Base.$op(t::T) where {T <: ITime} =
                ITime($op(t.counter), t.period, t.epoch)
        end,
    )
end

macro itime_binary_op(op)
    return esc(
        quote
            function Base.$op(t1::T1, t2::T2) where {T1 <: ITime, T2 <: ITime}
                t1p, t2p = promote(t1, t2)
                ITime($op(t1p.counter, t2p.counter), t1p.period, t1p.epoch)
            end
        end,
    )
end

macro itime_binary_op_notype(op)
    return esc(
        quote
            function Base.$op(t1::T1, t2::T2) where {T1 <: ITime, T2 <: ITime}
                t1p, t2p = promote(t1, t2)
                $op(t1p.counter, t2p.counter)
            end
        end,
    )
end

@itime_unary_op abs
@itime_unary_op -

@itime_binary_op +
@itime_binary_op -

Base.isnan(t::ITime) = Base.isnan(t.counter)

@itime_binary_op_notype <=
@itime_binary_op_notype >=
@itime_binary_op_notype isless
@itime_binary_op_notype ==
@itime_binary_op_notype isequal
@itime_binary_op_notype isapprox
@itime_binary_op_notype div
@itime_binary_op_notype //
Base.:/(t1::T1, t2::T2) where {T1 <: ITime, T2 <: ITime} = t1 // t2

# Multiplication/division by numbers
Base.div(t::T1, num::Integer) where {T1 <: ITime} =
    ITime(div(t.counter, num), t.period, t.epoch)
Base.:*(num::Integer, t::T) where {T <: ITime} =
    ITime(num * t.counter, t.period, t.epoch)
Base.:*(t::T, num::Integer) where {T <: ITime} =
    ITime(num * t.counter, t.period, t.epoch)

# Behave as a scalar when broadcasted
Base.Broadcast.broadcastable(t::ITime) = Ref(t)
