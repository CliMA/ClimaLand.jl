# CFTime.jl

In many Earth science disciplines and beyond, expressing a time instance and a duration is essential. The [CF conventions](http://cfconventions.org/Data/cf-conventions/cf-conventions-1.12/cf-conventions.html#calendar) provide a rich and flexible framework for handling time, equally applicable to observations and model data.

CFTime.jl implements the time structures standardized by the CF conventions, namely:

* Julian calendar ([`DateTimeJulian`](@ref)). A year is a leap year if it is divisible by 4. For example, 1900 and 2000 are both leap years in the Julian calendar.
* Proleptic Gregorian calendar ([`DateTimeProlepticGregorian`](@ref)). A year is a leap year if it is divisible by 4 but not 100 or if it is divisible by 400. For example, 1900 is not a leap year in the proleptic Gregorian calendar but 2000 is.
* The default mixed Gregorian/Julian calendar  ([`DateTimeStandard`](@ref)). This calendar uses the Julian calendar for time instances before 15th October 1582 and Gregorian calendar afterwards.
* A calendar without leap years ([`DateTimeNoLeap`](@ref)). All years are 365 days long.
* A calendar with only leap years ([`DateTimeAllLeap`](@ref)). All years are 366 days long.
* A calendar with every year being 360 days long (divided into 30-day months) ([`DateTime360Day`](@ref)).

The first three calendars (with different rules for leap years) can be used to express the time instances of observations or model. The remaining three calendars correspond to idealised model configurations where the duration of a year (revolution of the Earth around the Sun) is assumed to be exactly 365, 366 or 360 days.

While almost all datasets used in Earth Science use dates after the year 1582, some datasets or software systems use a time origin before this date, which makes it necessary to handle the transition from Julian to Gregorian calendar. Additionally, some dataset use microseconds and nanoseconds as time resolution, whereas Julia's [`Dates.DateTime`](https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.DateTime) has milliseconds as time resolution.

Note that time zones and leap seconds are currently not supported by `CFTime.jl`.


## Installation

Inside the Julia shell, you can download and install the package by issuing:

```julia
using Pkg
Pkg.add("CFTime")
```

CFTime is a pure Julia package and currently depends only on the modules [`Dates`](https://docs.julialang.org/en/v1/stdlib/Dates/) and [`Printf`](https://docs.julialang.org/en/v1/stdlib/Printf/), which are part of Julia’s standard library.


### Latest development version

If you want to try the latest development version, you can do this with the following commands:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/JuliaGeo/CFTime.jl", rev="master"))
```

## Types

[`DateTimeStandard`](@ref), [`DateTimeProlepticGregorian`](@ref), [`DateTimeJulian`](@ref), [`DateTimeNoLeap`](@ref), [`DateTimeAllLeap`](@ref) and [`DateTime360Day`](@ref) are the core types that define time instances according to the corresponding calendars. The main focus of this package is to instantiate and to manipulate these types.

```@docs
DateTimeStandard
DateTimeJulian
DateTimeProlepticGregorian
DateTimeAllLeap
DateTimeNoLeap
DateTime360Day
AbstractCFDateTime
```

## Time encoding and decoding

```@docs
CFTime.timedecode
CFTime.timeencode
```

## Accessor Functions

```@docs
CFTime.year(dt::AbstractCFDateTime)
CFTime.month(dt::AbstractCFDateTime)
CFTime.day(dt::AbstractCFDateTime)
CFTime.hour(dt::AbstractCFDateTime)
CFTime.minute(dt::AbstractCFDateTime)
CFTime.second(dt::AbstractCFDateTime)
CFTime.millisecond(dt::AbstractCFDateTime)
CFTime.microsecond(dt::AbstractCFDateTime)
CFTime.nanosecond(dt::AbstractCFDateTime)
CFTime.picosecond(dt::AbstractCFDateTime)
CFTime.femtosecond(dt::AbstractCFDateTime)
CFTime.attosecond(dt::AbstractCFDateTime)
```

## Query Functions


```@docs
daysinmonth
daysinyear
yearmonthday
yearmonth
monthday
firstdayofyear
dayofyear
```

## Conversion Functions


The flexibility of CFTime's datetime (related to the time origin, time resolution and type of the time counter) comes with some cost. When merging data from different sources, the resulting merged time vector may not have a concrete type, as there is no implicit conversion to a common time origin or internal unit, unlike Julia's [`DateTime`](https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.DateTime). In some cases, the user might decide to explicitly convert all times to a common time origin and internal unit for optimal performance.

The [`convert`](@ref) function can be used to convert dates between the different calendars:

```jldoctest convert
using CFTime, Dates
convert(DateTime,DateTimeJulian(2024,4,4))
# output
2024-04-17T00:00:00
```

```jldoctest convert
convert(DateTimeJulian,DateTime(2024,4,17))
# output
DateTimeJulian(2024-04-04T00:00:00)
```


```@docs
convert(::Type{DateTime}, dt::DateTimeStandard)
reinterpret
```

## Arithmetic

Adding and subtracting time periods is supported:

```jldoctest arithmetic
using CFTime, Dates
DateTimeStandard(1582,10,4) + Dates.Day(1)
# output
DateTimeStandard(1582-10-15T00:00:00)
```

1582-10-15 is the adoption of the Gregorian Calendar.

Comparison operator can be used to check if a date is before or after another date.

```jldoctest arithmetic
DateTimeStandard(2000,01,01) < DateTimeStandard(2000,01,02)
# output
true
```

## Ranges


Time ranges can be constructed using a start date, end date and a time increment:

```jldoctest ranges
using CFTime, Dates
range = DateTimeStandard(2000,1,1):Dates.Day(1):DateTimeStandard(2000,12,31)
length(range)
# output
366
```

```jldoctest ranges
step(range)
# output
1 day
```

Note that there is no default increment for range.

## Rounding

Julia's `DateTime` records the time relative to a time origin (January 1st, 1 BC or 0000-01-01 in ISO_8601) with a millisecond accuracy. Converting CFTime date time structures to
Julia's `DateTime` (using `convert(DateTime,dt)`) can trigger an inexact exception if the conversion cannot be done without loss of precision. One can use the `round` function in order to round to the nearest time represenatable by `DateTime`:

```jldoctest floor
using CFTime: DateTimeStandard
using Dates: DateTime
dt = DateTimeStandard(24*60*60*1000*1000 + 123,"microsecond since 2000-01-01")
round(DateTime,dt)
# output
2000-01-02T00:00:00
```

The functions `floor` and `ceil` are also supported. They can be used to effectively reduce the time resolution, for example:

```jldoctest floor
using CFTime: DateTimeStandard
using Dates

dt = DateTimeStandard(24*60*60,"second since 2000-01-01")

floor(dt+Second(9),Second(10)) == dt
# output
true
```

```jldoctest floor
ceil(dt+Second(9),Second(10)) == dt + Second(10)
# output
true
```


```jldoctest floor
round(dt+Second(9),Second(10)) == dt + Second(10)
# output
true
```


## Type-stable constructors

To create a type-stable date time structure, use the `DateTimeStandard` (and similar) either with the default `units` and time `origin`, a constant unit/origin or a value type of the unit and origin. For example:

```julia
using CFTime: DateTimeStandard

function foo(year,month,day,hour,minute,second)
   DateTimeStandard(year,month,day,hour,minute,second;
                    units=:second, origin=(1970,1,1))
end

# Type-stable thanks to constant propagation
@code_warntype foo(2000,1,1,0,0,0)


function foo2(year,month,day,hour,minute,second,units,origin)
   DateTimeStandard(year,month,day,hour,minute,second; units, origin)
end

# This not type-stable as the type depends on the value of units and origin
units = :second
origin = (1970,1,1)
@code_warntype foo2(2000,1,1,0,0,0,units,origin)


# But this is again type-stable
units = Val(:second)
origin = Val((1970,1,1))
@code_warntype foo2(2000,1,1,0,0,0,units,origin)
```

## Limitations and caveats

In Julia, integer overflows can occur in various operations such as:

```
using Dates
typemax(Int64)+1
# -9223372036854775808

DateTime(2000,1,1) + Dates.Day(typemax(Int64))
# 1999-12-31T00:00:00

DateTime(typemax(Int64),1,1)
# -0002-12-31T13:26:24
```

None of these calls issue a warning or an error in Julia (due to the performance impact affecting the vast majority of cases). The same approach is taken in CFTime: the user should be aware than under some (exeptional) circumstances, the default `Int64` can overflow. But in addition to Julia's `Dates`, CFTime gives us user the option to use different storage types like `Int128` or `BigInt`:

``` julia
using CFTime
DateTimeStandard(Int128,typemax(Int64),1,1)
# DateTimeStandard(9223372036854775807-01-01T00:00:00)
```

In this example, the year 9223372036854775807 (about nine quintillion years) is about 60 millions time longer than the age of the universe. It dependens on the users' application whether such long time frames are needed.

The default `Int64` time counter can also overflow when using a very small time resolution are used:


``` julia
using CFTime
DateTimeStandard(2000,1,1, units=:attosecond)
# which is the same as
DateTimeStandard(Int64,2000,1,1, units=:attosecond,origin=(1900,1,1))
# -> overflow
```

An `Int64` time counter does not allow to represent the number of attoseconds between the 1st January of the years 2000 and 1900. In fact, only 18 s around the time origin can be represent with a 64-bit integer and attoseconds as time unit.


This table sumarizes the maximal time span for `Int64` around the time origin (1900-01-01 00:00:00, per default)

| unit | time span |
|------|-----------|
| day | ±2.525e+16 years |
| hour | ±1.052e+15 years |
| minute | ±1.754e+13 years |
| second | ±2.923e+11 years |
| millisecond | ±2.923e+08 years |
| microsecond | ±2.923e+05 years |
| nanosecond | ±292.3 years |
| picosecond | ±106.8 days |
| femtosecond | ±2.562 hours |
| attosecond | ±9.223 seconds |


An `Int128` is sufficiently large for the case above:


```julia
DateTimeStandard(Int128,2000,1,1, units=:attosecond,origin=(1900,1,1))
# DateTimeStandard(2000-01-01T00:00:00)
```

Operations can typically be 2 to 4 times slower when using an `Int128` rather than an `Int64`.

It is important to note that when small durations are added to a time instances, then internal type promotion will also change the unit of the time instance (which is milliseconds since 1900-01-01):

```julia
dt = DateTimeStandard(2000,1,1) # default unit is :millisecond
dt2 = dt + Dates.Nanosecond(1) # now the unit is :nanosecond -> but no overflow
# DateTimeStandard(2000-01-01T00:00:00.000000001)

dt = DateTimeStandard(2160,1,1) # default unit is :millisecond
dt2 = dt + Dates.Nanosecond(1) # now the unit is :nanosecond -> overflow
# DateTimeStandard(1575-06-03T00:25:26.290448385)
```

Whether an overflow occurs or not cannot be statically inferred and a warning at runtime would impact performance. To avoid overflows, one should use `Int128` or `BigInt` in these cases.

Note that on a 32-bit system, the default internal time counter is still `Int64` as in Julia's `Dates` module.

## Internal API

For CFTime 0.1.4 and before all date-times are encoded using internally milliseconds since a fixed time origin and stored as an `Int64` similar to julia's `Dates.DateTime`.
However, this approach does not allow to encode time with a sub-millisecond precision allowed by the CF convention and supported by e.g. [numpy](https://numpy.org/doc/1.25/reference/arrays.datetime.html#datetime-units). While `numpy` allows attosecond precision, it can only encode a time span of ±9.2 around the date 00:00:00 UTC on 1 January 1970. In CFTime the time origin and the number containing the duration and the time precision are now encoded as two additional type parameters.

When wrapping a CFTime date-time type, it is recommended for performance reasons to make the containing structure also parametric, for example

``` julia
struct MyStuct{T1,T2}
  dt::DateTimeStandard{T1,T2}
end
```

Future version of CFTime might add other type parameters.
Internally, `T1` corresponds to a `CFTime.Period{T,Tfactor,Texponent}` structure  wrapping a number type `T` representing the duration expressed in seconds as:

```
duration * factor * 10^exponent
```

where `Tfactor` and `Texponent` are value types of `factor` and `exponent` respectively.

For example, duration 3600000 milliseconds is represented as `duration = 3600000`,
`Tfactor = Val(1)`, `Texponent = Val(-3)`, as

```
3600000 milliseconds = 3600000 * 1 * 10⁻³ seconds
```

or the duration 1 hours is `duration = 1`,  `Tfactor = Val(3600)` and `Texponent = Val(0)` since:

```
1 hour = 1 * 3600 * 10⁰ seconds
```

There is no normalization of the time duration per default as it could lead to under-/overflow.

The type parameter `T2` of `DateTimeStandard` encodes the time origin as a tuple of integers starting with the year (year,month,day,hour,minute,seconds,...attoseconds).  Only the year, month and day need to be specified; all other default to zero.
For example `T2` would be `Val((1970,1,1))` if the time origin is the 1st January 1970.

By using value types as type parametes, the time origin, time resolution... are known to the compiler.
For example, computing the difference between between two date time expressed in as the same time origin and units as a single substraction:

```julia
using BenchmarkTools
using Dates
using CFTime: DateTimeStandard

dt0 = DateTimeStandard(1,"days since 2000-01-01")
dt1 = DateTimeStandard(1000,"days since 2000-01-01")

difference_datetime(dt0,dt1) = Dates.Millisecond(dt1 - dt0).value
@btime difference_datetime($dt0,$dt1)

# output (minimum of 5 @btime trails)
# 1.689 ns (0 allocations: 0 bytes)

v0 = 1
v1 = 1000

difference_numbers(v0,v1) = (v1-v0)*(86_400_000)
@btime difference_numbers($v0,$v1)

# output (minimum of 5 @btime trails)
# 1.683 ns (0 allocations: 0 bytes)
```

The information in this section and any other information marked as internal or experimental is not part of the public API and not covered by the semantic versioning.
Future version of CFTime might add or changing the meaning of type parameters as patch-level changes. However, removing a type parameter would be considered as a breaking change.
