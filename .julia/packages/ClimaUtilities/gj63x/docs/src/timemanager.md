# TimeManager

`TimeManager` defines `ITime`, the time type for CliMA simulations, alongside
with various functions to work with it.

## `ITime`

`ITime` is a type to describe times and dates. The "I" in `ITime` stands for
_integer_: internally, `ITime` uses integers to represent times and dates,
meaning that operations with `ITime` are exact and do not occur in
floating-point errors.

`ITime` can be thought a combination of three quantities, a `counter`, a
`period`, and (optionally) an `epoch`, with the `counter` counting how many
`period`s have elapsed since `epoch`. In other words, `ITime` counts clock
cycles from an epoch, with each clock cycle lasting a `period`.

Another useful mental model for `ITime` is that it is a time with some units. It
is useful to keep this abstraction in mind when working with binary operations
that involve `ITime`s because it helps with determining the output type of the
operation (more on this later).

Let us start getting familiar with `ITime` by exploring some basic operations.

### First steps with `ITime`

The first step in using `ITime` is to load it. You can load `ITime` by loading
the `TimeManager` module as in the following:

```@example example1
using ClimaUtilities.TimeManager
```

This will load `ITime` as well as some other utilities. If you only wish to load
`ITime`, you can change `using` with `import` and explicitly ask for `ITime`

```@julia example 
import ClimaUtilities.TimeManager: ITime
```

In these examples, we will stick with `using ClimaUtilites.TimeManager` because
we will call other functions as well.

By default, `ITime` assumes that the clock cycle is one second:
```@example example1
ITime(5)
```
The output that is printed shows the time represented (5 seconds), and its break
down into the integer counter (5) and the duration of each clock cycle (1
second).

This `ITime` does not come with any date information attached to it. To add it,
you can pass the `epoch` keyword
```@example example1
using Dates
ITime(5; epoch = DateTime(2012, 12, 21))
```
Now, the output also reports `epoch` and current date (5 seconds after
midnight of the 21st of December 2012).

The period can also be customized with the `period` keyword. This keyword
accepts any `Dates.FixedPeriod` (fixed periods are intervals of time that have a
fixed duration, such as `Week`, but unlike `Month`), for example
```@example example1
time1 = ITime(5; period = Hour(20), epoch = DateTime(2012, 12, 21))
```
All these quantities are accessible with functions:
```@example example1
@show "time1" counter(time1) period(time1) epoch(time1) date(time1)
```

Now that we know how to create `ITime`s, we can manipulate them. `ITime`s
support most (but not all) arithmetic operations.

For `ITime`s with an `epoch`, it is only possible to combine `ITime`s that
have the same `epoch`. `ITime`s can be composed to other times and
arithmetic operations propagate `epoch`s.

This works:
```@example example1
time1 = ITime(5; epoch = DateTime(2012, 12, 21))
time2 = ITime(15; epoch = DateTime(2012, 12, 21))
time1 + time2
```
This works too
```@example example1
time1 = ITime(5; epoch = DateTime(2012, 12, 21))
time2 = ITime(15)
time1 + time2
```
This does not work because the two `ITimes` don't have the same `epoch`
```@example example1
time1 = ITime(5; epoch = DateTime(2012, 12, 21))
time2 = ITime(15; epoch = DateTime(2025, 12, 21))
try  #hide
time1 + time2
catch err; showerror(stderr, err); end  #hide
```

When dealing with binary operations, it is convenient to think of `ITime`s as
dimensional quantities (quantities with units). In the previous examples,
addition returns a new `ITime`s. Division will behave differently
```@example example1
time1 = ITime(15)
time2 = ITime(5)
time1 / time2
```
In this case, the return value is just a number (essentially, we divided 15
seconds by 5 seconds, resulting the dimensionless factor of 3).

Similarly, adding a number to an `ITime` is not allowed (because the number
doesn't have units), but multiplying by is fine.

In the same spirit, when `ITime`s with different `periods` are combined, units
are transformed to the
```@example example1
time1 = ITime(5; period = Hour(1))
time2 = ITime(1; period = Minute(25))
time1 - time2
```
As we can see, the return value of 5 hours - 25 minutes is 275 minutes and it is
represented as 55 clock cycles of a period of 5 minutes. `Minute(5)` was picked
because it the greatest common divisor between `Hour(1)` and `Minute(25)`.

At any point, we can obtain a date or a number of seconds from an `ITime`:
```@example example1
time1 = ITime(5; period = Day(10), epoch = DateTime(2012, 12, 21))
date(time1)
```
And
```@example example1
time1 = ITime(5; period = Day(10), epoch = DateTime(2012, 12, 21))
seconds(time1)
```
In this, note that the `seconds` function always returns a Float64.

In this section, we saw that `ITime`s can be used to represent dates and times
and manipulated in a natural way.

`ITime`s support another feature, fractional times, useful for further
partitioning an interval of time.

### Dealing with times that cannot be represented

Sometimes, one need to work with fractions of a `period`. The primary
application of this is timestepping loops, which are typically divided in stages
which are a fraction of a timestep.

In these cases, we round to the nearest amount representable as an `ITime`.
Furthermore, these cases are constrained to only when an `ITime` is multiplied
by a float between zero and one. These are the only cases we consider, because
the only instance in which the problem of handling fraction of a `period` is in
the timestepping loops.

See the examples below.
```@repl example1
time = ITime(1; period = Hour(1), epoch = DateTime(2012, 12, 21))
0.3 * time # round down
0.5 * time # round down
0.7 * time # round up
1.2 * time # error because 1.2 > 1.0
-0.2 * time # error because -0.2 < 0.0
```

If multiplication between an `ITime` and a float is needed, you most likely want
to cast the `ITime` into a float. For example, in functions that compute
tendencies, there are no operations that are of the form `t + a * dt` where `a`
is a float, so it is safe to cast the `ITime` into a float as it will not lead to
a loss of accuracy in time.

In cases where you compute an expression of the form `t + a * dt` where `a` is a
float or something similar, then rounding occurs. For packages like
ClimaTimeSteppers, this means that any time stepping stages will incur a slight
inaccuracy in the time if the period is not divisible by `dt`.

### How to use `ITime`s in packages

The two key interfaces to work with `ITime`s are
[`ClimaUtilities.TimeManager.seconds`](@ref) and
[`ClimaUtilities.TimeManager.date`](@ref), which respectively return the time as
Float64 number of seconds and the associated `DateTime`.

If you want to support `AbstractFloat`s and `ITime`s, some useful functions are
`float` (which returns the number of seconds as a Float64), `zero`, `one`, and
`oneunit`. The difference between `one` and `oneunit` is that the latter returns
a `ITime` type, while the former returns a number. This is not what happens with
`zero`, which returns an `ITime`.

You might have the temptation to just sprinkle `float` everywhere to transition
your code to using `ITime`s. Resist to this temptation because it might defy the
purpose of using integer times in the first place.

For typical codes and simulation, we recommend only setting `t_start` and
`t_end` with a `epoch`: the `epoch` will be propagated naturally to
all the other times involved in the calculations.

Regarding the question on what `period` to use: if there are natural periods
(e.g., you are dealing with an hourly diagnostic variable), use it, otherwise
you can stick with the default. The `period` can always be changed by setting it
in `t_start` and `t_end`.

We provide a constructor from floating point numbers to assist you in
transitioning your package to using `ITime`s. This constructor guesses and uses
the largest period that can be used to properly represent the data.
```@example example1
ITime(60.0)
```

Beside the rounding which can lead to different results compared to using
floats, the other change can come from using `float` to convert an ITime to a
floating point type. At the moment, `float(t::ITime)` returns a `Float64` which
can cause changes if the model is ran with Float32 instead of Float64.

!!! note "Compatibility with ClimaTimeSteppers"
    Only IMEXAlgorithms and SSPKnoth in ClimaTimeSteppers are compatible with
    ITime.

### Common problems

#### How do I make several `ITime`s have the same type?

One can use `promote` to make all the variables have the same type. See the
example below of making `t0` and `tf` have the same type.

```@repl example1
t0 = ITime(0.25; epoch = DateTime(2012, 12, 21))
tf = ITime(10; period = Hour(1), epoch = DateTime(2012, 12, 21))
typeof(t0) == typeof(tf)
t0_same_type, tf_same_type = promote(t0, tf)
typeof(t0_same_type) == typeof(tf_same_type)
```

Similarly, one can also make `ITime`s in an array have the same type with a bit
more effort.

```@repl example1
itime_array = [t0, tf]
eltype(itime_array)
same_type_itime_arr = [promote(itime_array...)...]
eltype(same_type_itime_arr)
```

#### How do I multiply by a number by an ITime?

One can use `float` to cast the `ITime` into a floating point number
representing the number of seconds since the epoch. Casting to a floating
point number is okay as long as the floating point number is not used to keep
track of time of the simulation.

```@repl example1
t = ITime(1, period = Minute(1), epoch = DateTime(2010))
0.1 * float(t) # 6 seconds
```

For more information, see the section [Dealing with times that cannot be represented](@ref).


## Developer notes

### Why not use `Dates` directly?

Periods in Julia's `Dates` are also implemented as integer counters attached to
a type that encode its units. So why not using them directly?

There are a few reasons why `ITime` is preferred over Julia's `Dates`:
- `Dates` do not support fractional intervals, which are needed for the
  timestepping loop;
- `Dates` only support the Gregorian calendar. `ITime` provides an abstraction
  layer that will allow us to support other calendars without changing packages;
- `Dates` only allow the given periods, but it often natural to pick the
  simulation timestep as `period`;
- Julia's `Dates` are not necessarily faster or better and, being part of the
  standard library, means that it is hard to improve them.

### Why not support `Rational`?

Handling the case when the counter is rational introduces more complexity than
necessary. For instance, it is unclear how to handle the case when a float is
not perfectly representable as a rational number. This case comes up when
examining the time stepping stages of ARS343 in ClimaTimeSteppers. In
particular, the value ``\gamma`` is the middle root of the polynomial ``6x^3 -
18x^2 + 9x - 1 = 0`` and irrational. One could approximate ``\gamma`` as a
rational number, but large integers for the numerator and denominator are needed
to approximate ``\gamma`` to high accuracy. For example, ``\gamma`` approximated
as a rational number with tolerance within machine epsilon of `Float64` is
`19126397//43881317`. This could lead to overflowing in either the numerator or
denominator as ``\gamma`` propagates through the code. Hence, rational numbers
are not considered for this reason.

## TimeManager API

```@docs
ClimaUtilities.TimeManager.ITime
ClimaUtilities.TimeManager.ITime(counter::Union{Integer, Rational}; period::Dates.FixedPeriod = Dates.Second(1), epoch = nothing)
ClimaUtilities.TimeManager.ITime(t; epoch = nothing)
ClimaUtilities.TimeManager.seconds
ClimaUtilities.TimeManager.counter
ClimaUtilities.TimeManager.period
ClimaUtilities.TimeManager.epoch
ClimaUtilities.TimeManager.date
Base.promote(ts::ClimaUtilities.TimeManager.ITime...)
Base.:(:)(start::ClimaUtilities.TimeManager.ITime, step::ClimaUtilities.TimeManager.ITime, stop::ClimaUtilities.TimeManager.ITime)
Base.mod(x::ClimaUtilities.TimeManager.ITime, y::ClimaUtilities.TimeManager.ITime)
Base.iszero(x::ClimaUtilities.TimeManager.ITime)
Base.length(x::ClimaUtilities.TimeManager.ITime)
Base.float(t::T) where {T <: ClimaUtilities.TimeManager.ITime}
Base.one(t::T) where {T <: ClimaUtilities.TimeManager.ITime}
Base.oneunit(t::T) where {T <: ClimaUtilities.TimeManager.ITime}
Base.zero(t::T) where {T <: ClimaUtilities.TimeManager.ITime}
Base.:*(t::ClimaUtilities.TimeManager.ITime, a::AbstractFloat)
```
