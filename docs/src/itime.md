## ITime

`ITime`, or _integer time_, is a time type used by CliMA simulations to keep
track of simulation time. For more information, refer to the
[TimeManager section](https://clima.github.io/ClimaUtilities.jl/stable/timemanager/)
in ClimaUtilities and the
[ITime section](https://clima.github.io/ClimaAtmos.jl/stable/itime/) in ClimaAtmos.

### How do I use ITime?

There are three fields to an `ITime` which are `counter`, `period`, and `epoch`.
An `ITime` represents the amount of time, computed as `counter * period`, that
has passed since the `epoch`. See the example below of constructing an `ITime`.

```@repl itime
using ClimaUtilities.TimeManager, Dates # ITime is from ClimaUtilities
t = ITime(3, period = Minute(1), epoch = DateTime(2008))
counter(t)
period(t)
epoch(t)
```

The `ITime` `t` represents that three minutes have elapsed since 2008. The
`epoch` is chosen to be the start date of the simulation. Note that the counter
must be an integer. This means no floating point error can occur!

However, what are the differences between the following constructions of
`ITime`?

```@repl itime
dt1 = ITime(60000, period = Millisecond(1), epoch = DateTime(2010))
dt2 = ITime(60, period = Second(1), epoch = DateTime(2010))
dt3 = ITime(1, period = Minute(1), epoch = DateTime(2010))
```

The `ITime`s `dt1`, `dt2`, and `dt3` represent the same quantity which is 1
minute since 2010. However, the periods are different between the `ITime`s. As a
result, `t1`, `t2`, and `t3` can only represent times in terms of milliseconds,
seconds, and minutes respectively.

Furthermore, `dt1`, `dt2`, and `dt3` are all different types.
```@repl itime
typeof(dt1)
typeof(dt2)
typeof(dt3)
```

This raises two concerns to keep in mind when working with `ITime`s. The period
in `ITime` represents how accurate time is kept tracked of. In the example
above, `dt3` can only keep track of time that is accurate up to a minute.  This
can result in a loss of precision in a simulation if time needs to be resolved
with a higher precision than a minute. For example, in the time stepping stages,
the result of `t + α * dt` may be rounded. However, one cannot use nanoseconds
for the period because the counter will be large. This runs the risk of integer
overflow in the counter. As such, one needs to be careful about what period is
chosen for the `ITime`.

!!! note "Maximum representable date"
    With a period of `Dates.Second(1)`, an epoch of `2008`, and a counter of
    type `Int32`, the maximum representable date is `2076-01-19T03:14:07`. With
    `Int64` instead for the type of the counter, the maximum representable date is
    `292277025-08-17T07:12:55.807`.

    With a period of `Dates.Millisecond(1)`, an epoch of `2008`, and a counter
    of type `Int32`, the maximum representable date is
    `2008-01-25T20:31:23.647`. With `Int64` instead for the type of the counter,
    the maximum representable date is `292277025-08-17T07:12:55.807`.

!!! note "No period provided"
    If no period is provided, then the constructor for `ITime` will assume the
    provided value is in seconds and choose a reasonable value for the period.
    See
    [TimeManager section](https://clima.github.io/ClimaUtilities.jl/stable/timemanager/)
    in ClimaUtilities for more information.

Next, the different types of `dt1`, `dt2`, and `dt3` can lead to problems as
most functions and structs expect a single type for time. This can be solved
by using `promote` which makes all the `ITime`s have the same type. See the
example below.

```@repl itime
dt1, dt2, dt3 = promote(dt1, dt2, dt3)
typeof(dt1) == typeof(dt2)
typeof(dt2) == typeof(dt3)
```

Finally, other useful functions to know about `ITime` are `date` and `float`.
The function `date` returns the date of `t` and the function `float` convert `t`
into a floating point number. See the examples below.

```@repl itime
t = ITime(60, period = Second(1), epoch = DateTime(2008))
date(t)
float(t)
```
