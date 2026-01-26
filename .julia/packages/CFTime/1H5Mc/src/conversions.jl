function datenum(::Type{T}, y, m, d) where {T <: AbstractCFDateTime}
    cm = _cum_month_length(T)
    # turn year equal to -1 (1 BC) into year = 0
    if (y < 0) && !_hasyear0(T)
        y = y + 1
    end

    if m < 1 || m > 12
        error("invalid month $(m)")
    end

    if d < 1 || d > (cm[m + 1] - cm[m])
        error("invalid day $(d) in $(@sprintf("%04d-%02d-%02d", y, m, d))")
    end

    return cm[end] * (y - 1) + cm[m] + (d - 1)
end

# Calendar with regular month-length

function findmonth(cm, t2)
    mo = length(cm)
    while cm[mo] > t2
        mo -= 1
    end
    return mo
end

@inline function datetuple_ymd(::Type{T}, timed_::Number) where {T <: AbstractCFDateTime}
    cm = _cum_month_length(T)
    y = fld(Int64(timed_), cm[end])
    t2 = Int64(timed_) - cm[end] * y

    # find month
    mo = findmonth(cm, t2)
    d = t2 - cm[mo]

    # day and year start at 1 (not zero)
    d = d + 1
    y = y + 1

    if (y <= 0) && !_hasyear0(T)
        y = y - 1
    end

    return (y, mo, d)
end

@inline _cum_month_length(::Type{DateTimeAllLeap}) =
    (0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366)

@inline _cum_month_length(::Type{DateTimeNoLeap}) =
    (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)

@inline _cum_month_length(::Type{DateTime360Day}) =
    (0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360)


function validargs(::Type{T}, arg...) where {T <: AbstractCFDateTime}
    return nothing
end


"""
    dt2 = reinterpret(::Type{T}, dt)

Convert a variable `dt` of type [`DateTime`](https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.DateTime), [`DateTimeStandard`](@ref), [`DateTimeJulian`](@ref),
[`DateTimeProlepticGregorian`](@ref), [`DateTimeAllLeap`](@ref), [`DateTimeNoLeap`](@ref) or
[`DateTime360Day`](@ref) into the date time type `T` using the same values for
year, month, day, minute, second, ... attosecond.
The conversion might fail if a particular date does not exist in the
target calendar.

For example, the difference of the 1 January 2000 in the Julian and
and the 1 January 2000 in the standard calendar is 13 days:


```julia
using CFTime, Dates
dt = DateTimeJulian(2000,1,1)
Dates.Day(dt - reinterpret(DateTimeStandard,dt))
# 13 days
```

"""
function reinterpret(::Type{T1}, dt::DateTime) where {T1 <: AbstractCFDateTime}
    return T1(
        Dates.year(dt), Dates.month(dt), Dates.day(dt),
        Dates.hour(dt), Dates.minute(dt), Dates.second(dt),
        Dates.millisecond(dt)
    )
end

function reinterpret(::Type{DateTime}, dt::T2) where {T2 <: AbstractCFDateTime}
    return DateTime(
        Dates.year(dt), Dates.month(dt), Dates.day(dt),
        Dates.hour(dt), Dates.minute(dt), Dates.second(dt),
        Dates.millisecond(dt)
    )
end

function reinterpret(::Type{T1}, dt::T2) where {T1 <: AbstractCFDateTime} where {T2 <: AbstractCFDateTime}
    return T1(datetuple(dt)...)
end

function parseDT(::Type{Tuple}, str)
    str = replace(str, "T" => " ")

    # remove Z time zone indicator
    # all times are assumed to be UTC anyway
    if endswith(str, "Z")
        str = str[1:(end - 1)]
    end


    negativeyear = str[1] == '-'
    if negativeyear
        str = str[2:end]
    end

    y, m, d, h, mi, s, subsec... =
    if occursin(" ", str)
        datestr, timestr = split(str, ' ')
        y, m, d = parse.(Int64, split(datestr, '-'))

        timestr, tz = if occursin("+", timestr)
            ts, tz = split(timestr, "+")
            ts, tz
        elseif occursin("-", timestr)
            ts, tz = split(timestr, "-")
            ts, string("-", tz)
        else
            timestr, "00:00"
        end
        if !all(iszero(parse.(Int, split(tz, ":"))))
            @warn "Time zones are currently not supported by CFTime. Ignoring Time zone information: $(tz)"
        end

        time_split = split(timestr, ':')
        h_str, mi_str, s_str =
        if length(time_split) == 2
            time_split[1], time_split[2], "00"
        else
            time_split
        end


        h = parse(Int64, h_str)
        mi = parse(Int64, mi_str)
        s, subsec... =
        if occursin('.', s_str)
            # seconds contain a decimal point, e.g. 00:00:00.0

            sec_str, subsec_str = split(s_str, '.', limit = 2)

            # seconds (whole part)
            s = parse(Int64, sec_str)

            # need to pad 59.99 as 59.990
            nparts = div(length(subsec_str), 3, RoundUp)

            # subsec are groups of 3 digits
            # milliseconds, microseconds...
            subsec =
                ntuple(nparts) do i
                istr = (3 * i - 2):min((3 * i), length(subsec_str))
                parse(Int64, subsec_str[istr]) * 10^(3 - length(istr))
            end

            (s, subsec...)
        else
            (parse(Int64, s_str), Int64(0))
        end

        (y, m, d, h, mi, s, subsec...)
    else
        y, m, d = parse.(Int64, split(str, '-'))
        (y, m, d, Int64(0), Int64(0), Int64(0), Int64(0))
    end

    if negativeyear
        y = -y
    end

    return (y, m, d, h, mi, s, subsec...)
end

function parseDT(::Type{DateTime}, str)
    t = parseDT(Tuple, str)
    t = t[1:min(length(t), 7)] # truncate smaller than microseconds
    return DateTime(t...)
end

function parseDT(::Type{DT}, str) where {DT <: AbstractCFDateTime}
    return DT(parseDT(Tuple, str)...)
end

# deprecated, but exported
function timeunits(::Type{DT}, units) where {DT}
    t0, factor, exponent = _timeunits(DT, units)

    exponent = exponent + 3

    if exponent >= 0
        plength = factor * Int64(10)^exponent
    else
        plength = factor // Int64(10)^(-exponent)
    end

    return t0, plength
end

function _timeunits(::Type{DT}, units) where {DT}
    tunit_mixedcase, starttime = strip.(split(units, " since "))
    tunit = lowercase(tunit_mixedcase)

    t0 = parseDT(DT, starttime)

    # make sure that plength is 64-bit on 32-bit platforms
    # plength is duration is *milliseconds*
    if (tunit == "years") || (tunit == "year")
        # SOLAR_YEAR is in ms
        return t0, SOLAR_YEAR, -3
    elseif (tunit == "months") || (tunit == "month")
        return t0, SOLAR_YEAR ÷ 12, -3
    else
        for (name, factor, exponent) in TIME_DIVISION
            if tunit == string(name, "s") || (tunit == string(name))
                return t0, factor, exponent
            end
        end
    end

    error("unknown units \"$(tunit)\"")
end


function timetype(calendar = "standard")
    DT =
    if (calendar == "standard") || (calendar == "gregorian")
        DateTimeStandard
    elseif calendar == "proleptic_gregorian"
        DateTimeProlepticGregorian
    elseif calendar == "julian"
        DateTimeJulian
    elseif (calendar == "noleap") || (calendar == "365_day")
        DateTimeNoLeap
    elseif (calendar == "all_leap") || (calendar == "366_day")
        DateTimeAllLeap
    elseif calendar == "360_day"
        DateTime360Day
    else
        error("Unsupported calendar: $(calendar)")
    end

    return DT
end

timetype(dt::Type{<:AbstractCFDateTime}) = dt
timetype(dt::Type{DateTime}) = dt

"""
    t0,plength = timeunits(units,calendar = "standard")

Parse time units (e.g. "days since 2000-01-01 00:00:00") and returns the start
time `t0` and the scaling factor `plength` in milliseconds.
"""
function timeunits(units, calendar = "standard")
    DT = timetype(calendar)
    return timeunits(DT, units)
end

# convert to Float64
_better_than_Float32(data::Float32) = Float64.(data)
_better_than_Float32(data) = data


function timedecode(::Type{DT}, data, units) where {DT <: AbstractCFDateTime}
    _convert(DTP, DDT, x) = DTP(DDT(x))
    _convert(DTP, DDT, x::Missing) = missing

    T = nonmissingtype(eltype(data))
    origintuple, factor, exponent = _timeunits(Tuple, units)
    DDT = Period{T, Val(factor), Val(exponent)}
    DTP = DT{DDT, Val(chop0(origintuple, 3))}

    return @. _convert(DTP, DDT, _better_than_Float32(data))
end


function timedecode(::Type{DateTime}, data, units)
    _convert(x, t0, plength) = t0 + Dates.Millisecond(round(Int64, plength * x))
    _convert(x::Missing, t0, plength) = missing

    t0, plength = timeunits(DateTime, units)
    return @. _convert(_better_than_Float32(data), t0, plength)
end


"""
    dt = timedecode(data,units,calendar = "standard"; prefer_datetime = true)

Decode the time information in data as given by the units `units` according to
the specified calendar.

`units` has the format `"DURATION_UNIT since TIME_ORIGIN"` where
`DURATION_UNIT` can be the year, month, $(join(getindex.(TIME_DIVISION, 1), ", ", " or "))
(in singular or plural). The string `TIME_ORIGIN` is a time instance written as
`year-month-day hour:minute:second.subseconds...`  using 24-hour clock system
where `subseconds` is the decimal fraction of the second. `day` and `hour` can
also be separated by the character `T` following ISO 8601.

Valid values for `calendar` are
`"standard"`, `"gregorian"`, `"proleptic_gregorian"`, `"julian"`, `"noleap"`, `"365_day"`,
`"all_leap"`, `"366_day"` and `"360_day"`.

If `prefer_datetime` is `true` (default), dates are
converted to the [`DateTime`](https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.DateTime) type (for the calendars
"standard", "gregorian", "proleptic_gregorian" and "julian")
unless the time unit is expressed in microseconds or smaller. Such conversion is
not possible for the other calendars.

| Calendar            | Type (prefer_datetime=true) | Type (prefer_datetime=false) |
| ------------------- | --------------------------- | ---------------------------- |
| `standard`, `gregorian` | `DateTime`                    | `DateTimeStandard`             |
| `proleptic_gregorian` | `DateTime`                    | `DateTimeProlepticGregorian`   |
| `julian`              | `DateTime`                    | `DateTimeJulian`               |
| `noleap`, `365_day`     | `DateTimeNoLeap`              | `DateTimeNoLeap`               |
| `all_leap`, `366_day`   | `DateTimeAllLeap`             | `DateTimeAllLeap`              |
| `360_day`             | `DateTime360Day`              | `DateTime360Day`               |

## Example:

```julia
using CFTime, Dates
# standard calendar
dt = CFTime.timedecode([0,1,2,3],"days since 2000-01-01 00:00:00")
# 4-element Vector{Dates.DateTime}:
#  2000-01-01T00:00:00
#  2000-01-02T00:00:00
#  2000-01-03T00:00:00
#  2000-01-04T00:00:00

dt = CFTime.timedecode([0,1,2,3],"days since 2000-01-01 00:00:00","360_day")
# 4-element Vector{DateTime360Day{CFTime.Period{Int64, Val{86400}(), Val{0}()}, Val{(2000, 1, 1)}()}}:
#  2000-01-01T00:00:00
#  2000-01-02T00:00:00
#  2000-01-03T00:00:00
#  2000-01-04T00:00:00
```

"""
function timedecode(data, units, calendar = "standard"; prefer_datetime = true)
    function datetime_convert(dt::AbstractCFDateTime{Period{T, Tfactor, Texponent}}) where
        {T, Tfactor, Texponent}

        return if unwrap(Texponent) >= -3
            # milliseconds, seconds, ...
            convert(DateTime, dt)
        else
            # do not convert microseconds or smaller
            dt
        end
    end
    datetime_convert(dt::Missing) = missing


    DT = timetype(calendar)
    dt = timedecode(DT, data, units)

    if (
            prefer_datetime &&
                (DT in [DateTimeStandard, DateTimeProlepticGregorian, DateTimeJulian])
        )

        return datetime_convert.(dt)
    else
        return dt
    end
end


"""
    data = timeencode(dt,units,calendar = "standard")

Convert a vector or array of [`DateTime`](https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.DateTime) (or [`DateTimeStandard`](@ref),
[`DateTimeProlepticGregorian`](@ref), [`DateTimeJulian`](@ref), [`DateTimeNoLeap`](@ref),
[`DateTimeAllLeap`](@ref), [`DateTime360Day`](@ref)) according to
the specified units (e.g. `"days since 2000-01-01 00:00:00"`) using the calendar
`calendar`.

`units` has the format `"DURATION_UNIT since TIME_ORIGIN"` where
`DURATION_UNIT` can be the year, month, $(join(getindex.(TIME_DIVISION, 1), ", ", " or "))
(in singular or plural). The string `TIME_ORIGIN` is a time instance written as
`year-month-day hour:minute:second.subseconds...`  using 24-hour clock system
where `subseconds` is the decimal fraction of the second. `day` and `hour` can
also be separated by the character `T` following ISO 8601.

Valid values for calendar are:
`"standard"`, `"gregorian"`, `"proleptic_gregorian"`, `"julian"`, `"noleap"`, `"365_day"`,
`"all_leap"`, `"366_day"`, `"360_day"`.


## Example:

```julia
using CFTime
dt = [DateTimeStandard(2000,1,1),DateTimeStandard(2000,1,2),DateTimeStandard(2000,1,3)]
CFTime.timeencode(dt,"days since 2000-01-01 00:00:00")
# output: [0., 1., 2.]
```

"""
function timeencode(
        data::AbstractArray{DT, N}, units,
        calendar = "standard"
    ) where {N} where {DT <: Union{DateTime, AbstractCFDateTime, Union{DateTime, AbstractCFDateTime, Missing}}}

    DT2 = timetype(calendar)
    t0, plength = timeunits(DT2, units)

    function encode(dt)
        if ismissing(dt)
            return missing
        end

        tmp =
        try
            convert.(DT2, dt)
        catch
            error("It is not possible to convert from $(DT) to $(DT2)")
        end

        return Dates.value(tmp - t0) / plength
    end
    return encode.(data)
end


function timeencode(
        data::DT, units,
        calendar = "standard"
    ) where {DT <: Union{DateTime, AbstractCFDateTime}}
    return timeencode([data], units, calendar)[1]
end


# do not transform data is not a vector of DateTime
# unused, should be removed
timeencode(data, units, calendar = "standard") = data


for CFDateTime in [
        :DateTimeStandard,
        :DateTimeJulian,
        :DateTimeProlepticGregorian,
    ]
    @eval begin
        """
            dt2 = convert(::Type{T}, dt)

        Convert a DateTime `dt` of type [`DateTimeStandard`](@ref), [`DateTimeProlepticGregorian`](@ref),
        [`DateTimeJulian`](@ref) or `DateTime` into the type `T` which can also be either
        [`DateTimeStandard`](@ref), [`DateTimeProlepticGregorian`](@ref), [`DateTimeJulian`](@ref) or [`DateTime`](https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.DateTime).

        Conversion is done such that duration (difference of DateTime types) are
        preserved. For dates on and after 1582-10-15, the year, month and days are the same for
        the types [`DateTimeStandard`](@ref), [`DateTimeProlepticGregorian`](@ref) and [`DateTime`](https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.DateTime).

        For dates before 1582-10-15, the year, month and days are the same for
        the types [`DateTimeStandard`](@ref) and [`DateTimeJulian`](@ref).
        """
        function convert(::Type{DateTime}, dt::$CFDateTime)
            origin = _origin_period(dt)
            ms = Dates.Millisecond(dt.instant + origin + DATETIME_OFFSET)
            return DateTime(UTInstant{Millisecond}(ms))
        end

        function convert(::Type{$CFDateTime}, dt::DateTime)
            T = $CFDateTime

            origin = DateTime(UTInstant{Millisecond}(Millisecond(0)))
            y, mdHMS... = (Dates.year(origin), Dates.month(origin), Dates.day(origin))

            if !_hasyear0(DateTimeProlepticGregorian) && y <= 0
                origintuple = (y - 1, mdHMS...)
            else
                origintuple = (y, mdHMS...)
            end

            p = convert(Period, dt.instant.periods)
            dt_greg = DateTimeProlepticGregorian{typeof(p), Val(origintuple)}(p)
            return convert(T, dt_greg)
        end


        # need to convert the time origin because not all origins are valid in
        # all calendars
        function convert(::Type{$CFDateTime}, dt::Union{DateTimeStandard, DateTimeProlepticGregorian, DateTimeJulian})

            T = $CFDateTime
            days, HMS... = timetuplefrac(_origin_period(dt))
            y, m, d = datetuple_ymd(T, days)

            origintuple = chop0((y, m, d, HMS...), 3)
            return T{typeof(dt.instant), Val(origintuple)}(dt.instant)
        end


    end
end


function convert(::Type{DT}, dt::AbstractCFDateTime) where {DT <: AbstractCFDateTime{T, Torigintuple}} where {T, Torigintuple}
    dt0 = DT(T(0))
    p = convert(T, dt - dt0)
    return DT(p)
end
