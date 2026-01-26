import Dates

import ClimaUtilities.TimeManager: ITime, date

"""
    seconds_to_str_short(time_seconds::Real)

Convert a time in seconds to a string representing the time in "short" units.

Examples:
========

```jldoctest
julia> import ClimaDiagnostics: seconds_to_str_short

julia> seconds_to_str_short(0.1)
"0.1s"

julia> seconds_to_str_short(60)
"1m"

julia> seconds_to_str_short(3600)
"1h"

julia> seconds_to_str_short(3601)
"1h_1s"
```
"""
function seconds_to_str_short(time_seconds::Real)
    name = ""
    days, rem_seconds = divrem(time_seconds, 24 * 60 * 60)
    hours, rem_seconds = divrem(rem_seconds, 60 * 60)
    minutes, seconds = divrem(rem_seconds, 60)

    # At this point, days, hours, minutes, seconds have to be integers.
    # Let us force them to be such so that we can have a consistent string output.
    days, hours, minutes = map(Int, (days, hours, minutes))

    round(seconds) ≈ seconds && (seconds = convert(Int, seconds))

    days > 0 && (name *= "$(days)d_")
    hours > 0 && (name *= "$(hours)h_")
    minutes > 0 && (name *= "$(minutes)m_")
    seconds > 0 && (name *= "$(seconds)s_")
    return rstrip(name, '_')
end

"""
    seconds_to_str_long(time_seconds::Real)

Convert a time in seconds to a string representing the time in "longer" units.

Examples:
========

```jldoctest
julia> import ClimaDiagnostics: seconds_to_str_long

julia> seconds_to_str_long(60)
"1 Minutes"

julia> seconds_to_str_long(3600)
"1 Hours"

julia> seconds_to_str_long(3601)
"1 Hours 1 Seconds"
```
"""
function seconds_to_str_long(time_seconds::Real)
    name = ""
    days, rem_seconds = divrem(time_seconds, 24 * 60 * 60)
    hours, rem_seconds = divrem(rem_seconds, 60 * 60)
    minutes, seconds = divrem(rem_seconds, 60)

    # At this point, days, hours, minutes, seconds have to be integers.
    # Let us force them to be such so that we can have a consistent string output.
    days, hours, minutes = map(Int, (days, hours, minutes))

    round(seconds) ≈ seconds && (seconds = convert(Int, seconds))

    days > 0 && (name *= "$(days) Days ")
    hours > 0 && (name *= "$(hours) Hours ")
    minutes > 0 && (name *= "$(minutes) Minutes ")
    seconds > 0 && (name *= "$(seconds) Seconds ")

    return rstrip(name, ' ')
end

"""
    time_to_date(time::AbstractFloat, start_date::Dates.DateTime)

Convert a `time` in seconds from `start_date` to a `Dates.DateTime`.

Examples:
========

```jldoctest
julia> import Dates: DateTime

julia> import ClimaDiagnostics: time_to_date

julia> start_date = DateTime(2024, 1, 1, 0, 0, 0);

julia> time_to_date(0.0, start_date)
2024-01-01T00:00:00

julia> time_to_date(0.5, start_date)
2024-01-01T00:00:00.500

julia> time_to_date(-1.0, start_date)
2023-12-31T23:59:59

julia> time_to_date(1.25, start_date)
2024-01-01T00:00:01.250
```
"""
function time_to_date(time::AbstractFloat, start_date::Dates.DateTime)
    # We go through milliseconds to allow fractions of a second (otherwise, Second(0.8)
    # would fail). Milliseconds is the level of resolution that one gets when taking the
    # difference between two DateTimes. In addition to this, we add a round to account for
    # floating point errors. If the floating point error is small enough, round will correct
    # it.
    time_ms = Dates.Millisecond(round(1_000 * time))
    return start_date + time_ms
end

"""
    time_to_date(time::ITime, start_date::Dates.DateTime)

Convert an `ITime` to a `Dates.DateTime`.

This function is needed to make EveryCalendarDtSchedule work when the type of
`integrator.t` is an `ITime`.
"""
function time_to_date(time::ITime, start_date::Dates.DateTime)
    return date(time)
end

"""
    period_to_str_short(time_seconds::Real)

Convert a time in seconds to a string representing the time in "short" units.

Examples:
========

# Examples
```jldoctest
julia> import ClimaDiagnostics: period_to_str_short

julia> period_to_str_short(Dates.Year(1))
"1y"

julia> period_to_str_short(Dates.Month(3))
"3M"

julia> period_to_str_short(Dates.Day(5))
"5d"

julia> period_to_str_short(Dates.Hour(2))
"2h"

julia> period_to_str_short(Dates.Minute(10))
"10m"

julia> period_to_str_short(Dates.Second(30))
"30s"
```
"""
function period_to_str_short(period::Dates.Period)
    return replace(
        string(period),
        r"year(s)?" => "y",
        r"month(s)?" => "M",
        r"day(s)?" => "d",
        r"hour(s)?" => "h",
        r"minute(s)?" => "m",
        r"second(s)?" => "s",
        " " => "",
        "," => "_",  # For Dates.CompundPeriod
    )
end

"""
    period_to_str_long(time_seconds::Real)

Convert a time in seconds to a string representing the time in "long" units.

Examples:
========

# Examples
```jldoctest
julia> import ClimaDiagnostics: period_to_str_long

julia> period_to_str_long(Dates.Year(1))
"1 Year"

julia> period_to_str_long(Dates.Month(2))
"2 Months"

julia> period_to_str_long(Dates.Day(5))
"5 Days"

julia> period_to_str_long(Dates.Hour(12))
"12 Hours"

julia> period_to_str_long(Dates.Minute(30))
"30 Minutes"

julia> period_to_str_long(Dates.Second(45))
"45 Seconds"
```
"""
function period_to_str_long(period::Dates.Period)
    period_str = replace(string(period), "," => "")
    return join(uppercasefirst.(split(period_str)), " ")
end
