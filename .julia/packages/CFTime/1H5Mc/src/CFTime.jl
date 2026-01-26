"""
`CFTime` encodes and decodes time units conforming to the Climate and Forecasting (CF) netCDF conventions.
Read more about the CF conventions [here](http://cfconventions.org/Data/cf-conventions/cf-conventions-1.12/cf-conventions.html#calendar).

For example:

```julia
using CFTime, Dates
# Decoding "360 day" calendar
dt = CFTime.timedecode([0,1,2,3],"days since 2000-01-01 00:00:00",DateTime360Day)
# Encoding
CFTime.timeencode(dt,"days since 2000-01-01 00:00:00",DateTime360Day)
```
The following types are supported [`DateTime`](https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.DateTime),
[`DateTimeStandard`](@ref), [`DateTimeJulian`](@ref), [`DateTimeProlepticGregorian`](@ref),
[`DateTimeAllLeap`](@ref), [`DateTimeNoLeap`](@ref) and [`DateTime360Day`](@ref)
"""
module CFTime

using Printf
import Dates
import Dates: UTInstant, Millisecond, DateTime, DateFormat
import Dates:
    year,
    month,
    day,
    hour,
    minute,
    second,
    millisecond,
    microsecond,
    nanosecond

import Dates:
    Year,
    Month,
    Day,
    Hour,
    Minute,
    Second,
    Millisecond,
    Microsecond,
    Nanosecond

import Dates:
    dayofyear,
    daysinmonth,
    daysinyear,
    firstdayofyear,
    len,
    monthday,
    yearmonth,
    yearmonthday,
    value

import Base: +, -, *, /, div, rem, lcm, gcd, gcdx, isless,
    abs, sign, signbit,
    string, show, convert, reinterpret, ==,
    promote_rule, floor, mod, round, typemax, typemin

include("constants.jl")
include("types.jl")
include("meeus_algorithm.jl")
include("period.jl")
include("datetime.jl")
include("accessors.jl")
include("query.jl")
include("conversions.jl")
include("ranges.jl")
include("rounding.jl")

export
    AbstractCFDateTime,
    DateTime360Day,
    DateTimeAllLeap,
    DateTimeJulian,
    DateTimeNoLeap,
    DateTimeProlepticGregorian,
    DateTimeStandard,
    datetuple,
    dayofyear,
    daysinmonth,
    daysinyear,
    firstdayofyear,
    monthday,
    timedecode,
    timeencode,
    timeunits,
    yearmonth,
    yearmonthday

function __init__()
    for CFDateTime in [
            DateTimeStandard,
            DateTimeJulian,
            DateTimeProlepticGregorian,
            DateTimeAllLeap,
            DateTimeNoLeap,
            DateTime360Day,
        ]
        Dates.CONVERSION_TRANSLATIONS[CFDateTime] = (
            Dates.Year, Dates.Month, Dates.Day,
            Dates.Hour, Dates.Minute, Dates.Second,
            Dates.Millisecond,
        )
    end
    return
end

end

#  LocalWords:  CFTime netCDF julia dt timedecode DateTime timeencode
#  LocalWords:  DateTimeStandard DateTimeJulian DateTimeAllLeap len
#  LocalWords:  DateTimeProlepticGregorian DateTimeNoLeap Printf jl
#  LocalWords:  UTInstant daysinmonth daysinyear yearmonthday isless
#  LocalWords:  yearmonth monthday dayofyear firstdayofyear accessors
#  LocalWords:  meeus AbstractCFDateTime datetuple timeunits init
#  LocalWords:  CFDateTime
