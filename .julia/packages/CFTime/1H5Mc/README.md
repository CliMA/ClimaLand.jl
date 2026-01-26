# CFTime

[![Build Status](https://github.com/JuliaGeo/CFTime.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaGeo/CFTime.jl/actions)
[![codecov](https://codecov.io/gh/JuliaGeo/CFTime.jl/graph/badge.svg?token=A6XMcOvIFr)](https://codecov.io/gh/JuliaGeo/CFTime.jl)
[![documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliageo.github.io/CFTime.jl/stable/)
[![documentation latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliageo.github.io/CFTime.jl/latest/)


`CFTime` encodes and decodes time units conforming to the [Climate and Forecasting (CF) conventions](https://cfconventions.org/).
`CFTime` was split out of the [NCDatasets](https://github.com/JuliaGeo/NCDatasets.jl) Julia package.

Features of CFTime include:

* Time instances as defined [Climate and Forecasting (CF) conventions](https://cfconventions.org/)
* Supporting a wide range of the time resolutions, from days down to attoseconds (for feature parity with NumPy's date time type)
* Supporting arbitrary time origins. For CFTime.jl, the time origin is part of the parametric type definition and not an additional field of the time data structure. As a consequence, a large array of date times with common time origin only need to store the time counter (64-bit integer by default) for every element, which makes this case as memory efficient as NumPy's or Julia's default date time for this common use case.
* By default, the time counter is a 64-bit integer, but other integers types (such as `Int32`, `Int128` or `BigInt`) or floating-point types can be used. Using an integer to encode a time instance is recommended for most applications, as it makes reasoning about the time resolution easier.
* Basic arithmetic such as computing the duration between two time instances
* Conversion function between CFTime types and Julia's `DateTime`.
* Regular time range based on Julia's range type. A time range is a vector of date time elements, but only the start time, the end time and the steps need to be stored in memory.

CFTime currently does not support leap seconds, which were standardized as part of CF conventions version 1.12, released in December 2024 and time zones.

## Installation

Inside the [Julia](https://julialang.org/) shell, you can download and install the package by issuing:

```julia
using Pkg
Pkg.add("CFTime")
```

After installing the package, the test suite of `CFTime` can be run using:

```julia
Pkg.test("CFTime")
```

## Example

For the [Climate and Forecasting (CF) conventions](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.12/cf-conventions.html#time-coordinate-units), time is expressed as duration since starting time. The function `CFTime.timedecode` allows to convert these
time instances as a Julia `DateTime` structure:

```julia
using CFTime, Dates

# standard calendar

dt = CFTime.timedecode([0,1,2,3],"days since 2000-01-01 00:00:00")
# 4-element Vector{Dates.DateTime}:
#  2000-01-01T00:00:00
#  2000-01-02T00:00:00
#  2000-01-03T00:00:00
#  2000-01-04T00:00:00
```


The function `CFTime.timeencode` does the inverse operation: converting a Julia `DateTime` structure to a duration since a start time:

```julia
CFTime.timeencode(dt,"days since 2000-01-01 00:00:00")
# 4-element Vector{Float64}:
#  0.0
#  1.0
#  2.0
#  3.0
```

The CF conventions also allow for different calendars, for example a calendar where every months has a duration of 30 days:

```julia
dt = CFTime.timedecode([0,1,2,3],"days since 2000-01-01 00:00:00",DateTime360Day)
# 4-element Vector{DateTime360Day{CFTime.Period{Int64, Val{86400}(), Val{0}()}, Val{(2000, 1, 1)}()}}:
#  2000-01-01T00:00:00
#  2000-01-02T00:00:00
#  2000-01-03T00:00:00
#  2000-01-04T00:00:00

CFTime.timeencode(dt,"days since 2000-01-01 00:00:00",DateTime360Day)
# 4-element Vector{Float64}:
#  0.0
#  1.0
#  2.0
#  3.0
```
You can replace in the example above the type `DateTime360Day` by the string `"360_day"` (the [name for the calendar](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.12/cf-conventions.html#calendar) according to the CF conventions).

Single time instances can also be created by calling the corresponding constructor function, e.g. `DateTimeStandard` for the standard calendar (mixed Gregorian/Julian calendar)
in a similar way than Julias `DateTime` type.
The `units` argument specifies the time resolutions (either `day`, `hour`, ... `attosecond`) for the common case where the duration is specified as an integer.
For example, the 1 January 2000 + 1 ns would be:

```julia
y,m,d = (2000,1,1)
hour,minute,sec = (0,0,0)
msec,µsec,nsec = (0,0,1)
DateTimeStandard(y,m,d,hour,minute,sec,msec,µsec,nsec; units=:nanosecond)
# DateTimeStandard(2000-01-01T00:00:00.000000001)
```

As in Julia's `DateTime`, the default time resolution is milliseconds.
The duration are encoded internally as a 64-bit signed integer. High precision integer (or floating point numbers) can also be used, for example a 128-bit signed integer:


```julia
DateTimeStandard(Int128,y,m,d,hour,minute,sec,msec,µsec,nsec; units=:nanosecond)
```

The default time origin is currently 1 January 1900 00:00:00. A different time origin can be used by setting the origin parameter:

```julia
DateTimeStandard(Int128,y,m,d,hour,minute,sec,msec,µsec,nsec; units=:nanosecond, origin=(1970,1,1))
```

The units and origin argument can be wrapped as a `Val` to ensure that these values are known at compile-time:

```julia
DateTimeStandard(Int128,y,m,d,hour,minute,sec,msec,µsec,nsec; units=Val(:nanosecond), origin=Val((1970,1,1)))
```

Several compile-time optimization have been implemented for the particular but common case where date have the same time origin and/or the same time resolution.

Arithmetic operations (`+`,`-`) and comparision operators on these types are supported, for example:


```julia
DateTimeStandard(2000,1,2) - DateTimeStandard(2000,1,1)
# 86400000 milliseconds

Dates.Day(DateTimeStandard(2000,1,2) - DateTimeStandard(2000,1,1))
# 1 day

DateTime360Day(2000,1,1) + Dates.Day(360)
# DateTime360Day(2001-01-01T00:00:00)

DateTimeStandard(2000,1,2) > DateTimeStandard(2000,1,1)
# true
```


## Parsing dates

Dates can be parsed by using `dateformat` from Julia's `Dates` module, for example:

```julia
dt = DateTimeNoLeap("21001231",dateformat"yyyymmdd");
# or
# dt = parse(DateTimeNoLeap,"21001231",dateformat"yyyymmdd")
Dates.year(dt),Dates.month(dt),Dates.day(dt)
# output (2100, 12, 31)
```

## Community guidelines 

We aim to follow the [community standards](https://julialang.org/community/standards/) of the Julia project.
 
Before creating an issue, check whether the issue has not already been reported.
When you file an issue, please include sufficient information that would _allow somebody else to reproduce the issue_, in particular:
1. Provide the code that generates the issue.
3. Make your code as simple as possible (while still showing the error and being runnable).
4. The full error message that you are seeing (in particular file names and line numbers of the stack-trace).
5. Which version of Julia and `CFTime` are you using? Please include the output of:

```julia
versioninfo()
using Pkg
Pkg.installed()["CFTime"]
```

For questions about the usage of CFTime, please create a topic in the Geo section of the Julia forum
https://discourse.julialang.org/c/domain/geo/.
Feel free to include `CC @Alexander-Barth` in your topic.

Contributions such as code or documentation are very welcomed! Here is an overview of the involved steps:

* fork the repository 
* create a new branch with a meaningful name 
* implement the new feature or fix the issue 
* add a test case and documentation for new features 
* use the [formatting tool runic](https://github.com/fredrikekre/Runic.jl?tab=readme-ov-file#quick-start)
* very that all tests pass locally 
* open the pull request with a brief description of the improvements to the code

For larger changes, it is recommended to open an issue first to better explain the purpose of the pull request.

## Trivia

[Microsoft Excel incorrectly considers 1900 as leap year](https://learn.microsoft.com/en-us/troubleshoot/microsoft-365-apps/excel/wrongly-assumes-1900-is-leap-year). For example, the difference between the dates 1-Jan-2000 and 1-Jan-1900 are off by one day. This bug was subsequently preserved for compatibility and codified in the ISO/IEC 29500:2008 standard.

## Alternatives

Julia packages:

 * [NanoDates.jl](https://github.com/JuliaTime/NanoDates.jl): Dates with nanosecond resolved days
 * [TimesDates.jl](https://github.com/JeffreySarnoff/TimesDates.jl): Nanosecond resolution for Time and Date, TimeZones
 * [AstroTime.jl](https://github.com/JuliaAstro/AstroTime.jl): Astronomical time keeping in Julia

Outside of the julia ecosystem:

* [cftime](https://unidata.github.io/cftime/) for python
* [CFtime](https://CRAN.R-project.org/package=CFtime) for R
* [cftime-rs](https://github.com/antscloud/cftime-rs) for rust with python bindings

## Acknowledgments

Thanks to Jeff Whitaker and [contributors](https://github.com/Unidata/cftime/graphs/contributors) for python's [cftime](https://github.com/Unidata/cftime) released under the MIT license which has helped the development of this package by providing reference values and a reference implementation for tests. The algorithm is based on Jean Meeus' algorithm published in Astronomical Algorithms (2nd Edition, Willmann-Bell, p. 63, 1998) adapted to years prior to 300 AC.
