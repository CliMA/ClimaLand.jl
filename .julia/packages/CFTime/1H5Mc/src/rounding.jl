"""
    dtr = round(::Type{DateTime}, dt::Union{DateTimeProlepticGregorian,DateTimeStandard,DateTimeJulian},r = RoundNearestTiesUp)

Round the date time `dt` to the nearest date time represenatable by julia's
[`DateTime`](https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.DateTime) using the rounding mode `r` (either
[`RoundNearest`](https://docs.julialang.org/en/v1/base/math/#Base.Rounding.RoundNearest) (default),
[`RoundDown`](https://docs.julialang.org/en/v1/base/math/#Base.Rounding.RoundDown), or
[`RoundUp`](https://docs.julialang.org/en/v1/base/math/#Base.Rounding.RoundUp)).
"""
function Base.round(::Type{DateTime}, dt::DateTimeProlepticGregorian, r::RoundingMode = RoundNearest)
    function round_ms(t)
        t_ms =
        if 3 + _exponent(t) > 0
            t.duration * _factor(t) * 10^(3 + _exponent(t))
        else
            t.duration * _factor(t) / 10^(-3 - _exponent(t))
        end

        t_ms_rounded =
        if t_ms isa Integer
            Int64(t_ms)
        else
            round(Int64, t_ms, r)
        end
        return Dates.Millisecond(t_ms_rounded)
    end

    origintuple = _origintuple(dt)
    origin = _origin_period(dt)

    if length(origintuple) >= 7
        # origin needs to be taken into account when rounding
        duration = round_ms(dt.instant + origin) + DATETIME_OFFSET
    else
        # better accuracy if dt.instant is a Float
        duration = round_ms(dt.instant) + origin + DATETIME_OFFSET
    end

    return DateTime(UTInstant{Millisecond}(Dates.Millisecond(duration)))
end

function Base.round(::Type{DateTime}, dt::Union{DateTimeJulian, DateTimeStandard}, r::RoundingMode = RoundNearest)
    return round(DateTime, convert(DateTimeProlepticGregorian, dt))
end

# see this discussion about changing the type parameters
# https://discourse.julialang.org/t/get-new-type-with-different-parameter/37253

for CFDateTime in (
        DateTimeStandard, DateTimeProlepticGregorian, DateTimeJulian,
        DateTimeNoLeap, DateTimeAllLeap, DateTime360Day,
    )
    @eval function Base.floor(dt::$CFDateTime, p::Period)
        origintuple = _origintuple(dt)
        origin = _origin_period(dt)
        t = dt.instant + origin
        t_mod = (t - mod(t, p)) - origin
        return $CFDateTime{typeof(t_mod), Val(origintuple)}(t_mod)
    end
end

function Base.floor(dt::AbstractCFDateTime, p::Dates.TimePeriod)
    return floor(dt, convert(Period, p))
end

function Base.floor(p::Period{T, Tfactor, Texponent}, ::Type{TDP}) where {T, Tfactor, Texponent} where {TDP <: Union{Dates.DatePeriod, Dates.TimePeriod}}
    scale =
    if TDP(1) > Dates.Second(1)
        1 // (TDP(1) ÷ Dates.Second(1))
    else
        (Dates.Second(1) ÷ TDP(1)) // 1
    end

    if _exponent(p) > 0
        d = (p.duration * _factor(p) * T(10)^_exponent(p) * scale)
    else
        d = (p.duration * _factor(p) * scale) ÷ (T(10)^(-_exponent(p)))
    end

    return TDP(floor(Int64, d))
end
