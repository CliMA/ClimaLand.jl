# TimeVaryingInputs.jl
#
# This module contains structs and methods to process external data and evaluate it on the
# model. This module only concerns with evaluations in time, not in space.

# There are three possible sources of data:
# 1. Analytic functions that prescribe how a variable has to be set at a given time
# 2. 0D, single-site data, which is assumed small enough to be saved to memory
# 3. 2/3D, global data, which cannot be saved to memory in its entirety.
#
# The TimeVaryingInputs module introduces a shared interface for the three cases so that
# uses and developers do not have to worry about the details of what type of data will be
# provided. Behind the scenes, we introduce a new type AbstractTimeVaryingInput, that has
# three concrete implementation, corresponding to the three use cases described above.
# Constructors will automatically identify which of the three implementations to use based
# on the input data, and the existence of three concrete structs should be considered an
# implementation detail.
#
# The three TimeVaryingInputs are:
# - AnalyticTimeVaryingInput,
# - InterpolatingTimeVaryingInput0D,
# - InterpolatingTimeVaryingInput23D.
#
# Along side these TimeVaryingInputs, we also define InterpolationMethods that implement
# specific interpolation strategies (e.g., linear interpolation).
#
# In all cases, the TimeVaryingInputs work with simulation time (ie, seconds from the
# beginning of the reference date). It is up to the various TimeVaryingInputs to convert this
# information to an actual date (if needed).

module TimeVaryingInputs

import Dates: Year
import Dates: DateTime, DatePeriod, Date

"""
    AbstractTimeVaryingInput

Note
=====

`TimeVaryingInput`s should be considered implementation details. The exposed public interface
should only be considered
- `TimeVaryingInput(input; method, context)` for construction,
- `evaluate!(dest, input, time)` for evaluation
"""
abstract type AbstractTimeVaryingInput end

"""
    AbstractInterpolationMethod

Defines how to perform interpolation.

Not all the TimeVaryingInputs support all the interpolation methods (e.g., no interpolation
methods are supported when the given function is analytic).

`AbstractInterpolationMethod`s have to implement a `extrapolation_bc` field.
"""
abstract type AbstractInterpolationMethod end

"""
    AbstractInterpolationBoundaryMethod

Defines how to handle values outside of the data boundary.

Not all the `AbstractInterpolationMethod` support all the `AbstractInterpolationBoundaryMethod`s.
"""
abstract type AbstractInterpolationBoundaryMethod end

"""
    Throw

Throw an error when interpolating outside of range.
"""
struct Throw <: AbstractInterpolationBoundaryMethod end

"""
    PeriodicCalendar

Repeat data periodically.

`PeriodicCalendar` has two modes of operation:

First, when provided with a `period` (described a `DatePeriod`, e.g., `Dates.Month(1)` or
`Dates.Year(1)`), assume that the provided data is repeated over that calendar period. A
`date` can be passed too, indicating what data to use. Only simple periods (e.g.,
`Dates.Month(1)`) are supported. When provided a period, a `repeat_date` is required too.
This is the period of time that is repeated. For example, if `period = Dates.Month(1)` and
`repeat_date = Dates.Date(1993, 11)`, November 1993 is repeated (if available in the input
data).

!!! note

    Passing a period is not supported by all the interpolators (e.g., when reading from 1D
    files).

Second, if no period is provided, when interpolating outside of range, restart from the beginning.

For example, if the data is defined from t0 = 0 to t1 = 10, extrapolating at t=13 is
equivalent to interpolating at t=2. In practice, we identify `t1 + dt` to be `t0` again.
This is different from what you might be used to for periodic boundary conditions, where the
identification is `t1 = t0`.

This second mode of operation `PeriodicCalendar` requires data to be uniformly sampled in time.

If the data is defined on a calendar year, this second mode of operation is equivalent to
using the first mode with `period = Dates.Year` (same with other periods).
"""
struct PeriodicCalendar{
    P <: Union{Nothing, DatePeriod},
    D <: Union{Nothing, DateTime, Date},
} <: AbstractInterpolationBoundaryMethod
    period::P
    repeat_date::D

    function PeriodicCalendar(
        period::Union{Nothing, DatePeriod},
        repeat_date::Union{Nothing, DateTime, Date},
    )
        if period isa DatePeriod && period.value != 1
            error("Only simple periods are supported (e.g., Month(1))")
        end
        return new{typeof(period), typeof(repeat_date)}(period, repeat_date)
    end
end

function PeriodicCalendar()
    return PeriodicCalendar(nothing, nothing)
end

"""
    Flat

When interpolating outside of range, use the boundary value.

For example, if the data is defined from t0 = 0 to t1 = 10, extrapolating at t=13 returns
the value at t1 = 10. When interpolating at t=-3, use t0 = 0.
"""
struct Flat <: AbstractInterpolationBoundaryMethod end

"""
    TimeVaryingInput(func)
    TimeVaryingInput(times, vals; method, context)

Construct on object that knows how to evaluate the given function/data on the model times.

When passing a function
=======================

When a function `func` is passed, the function has to be GPU-compatible (e.g., no splines).

When passing single-site data
=============================

When a `times` and `vals` are passed, `times` have to be sorted and the two arrays have to
have the same length.

=======
When the input is a function, the signature of the function can be `func(time, args...;
kwargs...)`. The function will be called with the additional arguments and keyword arguments
passed to `evaluate!`. This can be used to access the state and the cache and use those to
set the output field.

For example:
```julia
CO2fromp(time, Y, p) = p.atmos.co2
input = TimeVaryingInput(CO2fromY)
evaluate!(dest, input, t, Y, p)
```
"""
function TimeVaryingInput end

"""
    evaluate!(dest, input, time, args...; kwargs...)

Evaluate the `input` at the given `time`, writing the output in-place to `dest`.

Depending on the details of `input`, this function might do I/O and communication.

Extra arguments
================

`args` and `kwargs` are used only when the `input` is a non-interpolating function, e.g.,
an analytic one. In that case, `args` and `kwargs` are passed down to the function itself.
"""
function evaluate! end

"""
    extrapolation_bc(aim::AbstractInterpolationMethod)

Return the interpolation boundary conditions associated to `aim`.
"""
function extrapolation_bc(aim::AbstractInterpolationMethod)
    return aim.extrapolation_bc
end

"""
    NearestNeighbor(extrapolation_bc::AbstractInterpolationBoundaryMethod)

Return the value corresponding to the point closest to the input time.

`extrapolation_bc` specifies how to deal with out of boundary values.
The default value is `Throw`.
"""
struct NearestNeighbor{BC <: AbstractInterpolationBoundaryMethod} <:
       AbstractInterpolationMethod
    extrapolation_bc::BC
end

function NearestNeighbor()
    return NearestNeighbor(Throw())
end

"""
    LinearInterpolation(extrapolation_bc::AbstractInterpolationBoundaryMethod)

Perform linear interpolation between the two neighboring points.

`extrapolation_bc` specifies how to deal with out of boundary values.
The default value is `Throw`.
"""
struct LinearInterpolation{BC <: AbstractInterpolationBoundaryMethod} <:
       AbstractInterpolationMethod
    extrapolation_bc::BC
end

function LinearInterpolation()
    return LinearInterpolation(Throw())
end

"""
    LinearPeriodFillingInterpolation(period::Dates.DatePeriod,
                                     extrapolation_bc::AbstractInterpolationBoundaryMethod)

Perform a period-filling linear interpolation between the two neighboring points.

This interpolation is period-filling because it ensures that there is no gap on input data
larger than the given period.

This interpolation method is best explained with an example. Let us assume data is defined
on 03/01/1985, 05/01/1985, 11/02/1985, ... 17/12/1985, 07/01/1995, 11/01/1995, 14/02/1995,
..., 19/12/1995. In this case, data is defined every 10 years (date format is DD/MM/YYYY).
You can think of `LinearPeriodFillingInterpolation` as first filling the years in between by
performing linear interpolation across years, then performs a second linear interpolation
for the specific day (in actuality, the order of operations is reversed). In this example,
suppose first we want to evaluate the input on 04/01/1985. Since we have data for 1985, we
simply perform linear interpolation between 03/01/1985 and 05/01/1985. Suppose now, we want
to evaluate on 08/01/1987. We do not have data for the year 1987, so, we perform three
linear interpolations, first to evaluate the data on 08/01/1985 and 08/01/1995 (using the
neighboring points), and then between the resulting two values.

`period` specifies the period within which we want to preserve time variations. For the
above example, `period` would be `Year(1)` as we want to preserve the seasonal cycle within
a year. For a different example where we want to preserve the diurnal cycle, `period` would
be `Day(1)`.

`extrapolation_bc` specifies how to deal with out of boundary values. The default value is
`Throw`, meaning that extrapolation is not allowed.

When is this interpolation method useful?
==========================================

This interpolation method is useful when certain variations are more important than others.
For example, suppose you have hourly data defined once a month, on the 15th. A naive
interpolation method would interpolate use data from 15 Jan at 23.00 and 15 Feb at 00.00 to
interpolate on any day/hour of the days in between 15 Jan and 15 Feb. However, this would be
remove the diurnal cycle. `LinearPeriodFillingInterpolation` solves this problem.

Implementation details
======================

How is `LinearPeriodFillingInterpolation` actually implemented?

First, let get a high level perspective. Let us continue with the example above.

There are two cases:

First, the target date falls within other dates that are boundary dates for the two
neighboring years. An example is 08/01/1986 because both 08/01/1985 and 08/01/1995 have two
neighboring dates (05/01/1985, 11/02/1985 and 07/01/1995, 11/01/1995). This is the easy
case: we interpolate 05/01/1985 and 11/02/1985 to get 08/01/1985, and 07/01/1995 and
11/01/1995 to get 08/01/1995, and interpolate the two results to get 08/01/1986.

Second, the other cases, when the target date does not have two boundary dates defined in
1995. An example is 04/01/1986 because 07/01/1995 is the earliest date available on that
year. We have to be careful with this interpolation because we do not want to mix December
1985 with January 1995. In this case, we perform another interpolation to obtain a date in
December 1994 and reduce to the first case using this date as second half of the
interpolation bracket. To find what date to interpolate to, we walk the available dates
backwards in time at find the first date that we can interpolate because it falls into the
first case. In this example, it is 17/12/1994 (which we can obtain because we have data for
17/12/1985 and we can interpolate data on 17/12/1995).
"""
struct LinearPeriodFillingInterpolation{
    BC <: AbstractInterpolationBoundaryMethod,
} <: AbstractInterpolationMethod
    period::DatePeriod
    extrapolation_bc::BC

    function LinearPeriodFillingInterpolation(
        period::DatePeriod = Year(1),
        bc::AbstractInterpolationBoundaryMethod = Throw(),
    )
        period.value == 1 ||
            error("Only simple periods are supported (e.g., Year(1))")
        period isa Year ||
            @warn("Only `period = Year(1)` has been tested, check results!")
        bc isa PeriodicCalendar && error(
            "LinearPeriodFillingInterpolation is incompatible with PeriodicCalendar",
        )
        return new{typeof(bc)}(period, bc)
    end
end

extension_fns = [
    :ClimaCore => [:TimeVaryingInput, :evaluate!],
    :NCDatasets => [:TimeVaryingInput, :evaluate!],
    :CUDA => [:TimeVaryingInput, :evaluate!],
]

"""
    is_pkg_loaded(pkg::Symbol)

Check if `pkg` is loaded or not.
"""
function is_pkg_loaded(pkg::Symbol)
    return any(k -> Symbol(k.name) == pkg, keys(Base.loaded_modules))
end

function __init__()
    # Register error hint if a package is not loaded
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(
            MethodError,
        ) do io, exc, _argtypes, _kwargs
            for (pkg, fns) in extension_fns
                if Symbol(exc.f) in fns && !is_pkg_loaded(pkg)
                    if pkg == :CUDA
                        print(io, "\nIf you are using a GPU, import CUDA.")
                    else
                        print(io, "\nImport $pkg to enable `$(exc.f)`.";)
                    end
                end
            end
            if Symbol(exc.f) == :TimeVaryingInput
                print(
                    io,
                    "\nYou might also need a regridder to use `$(exc.f)`.";
                )
            end
        end
    end
end

end
