# [What do I have to do to use `ClimaDiagnostics`?](@id user_guide_header)

In this page, we describe the low level interface that `ClimaDiagnostics` offers
to work with diagnostics. Most packages implement addition interface to
streamline computing and outputting diagnostics, so you should first refer to
the manual of your package of interest. Come back here if you want to go beyond
what the package developers offer and unlock the full power of
`ClimaDiagnostics`.

There are two fundamental objects in `ClimaDiagnostics`: the
`DiagnosticVariable`, and the `ScheduledDiagnostic`.

#### `DiagnosticVariable`s

A `DiagnosticVariable` is a recipe on how to compute something alongside with
some metadata.

For example, a `DiagnosticVariable` might be the air temperature. In pseudocode,
some of the information we might want to include attach to the air temperature
are:
```yaml
short_name: "ta"
long_name: "Air Temperature"
units: "K"
how_to_compute: state.ta
...
```

Conceptually, a `DiagnosticVariable` is a variable we know how to compute from the state.
We attach more information to it for documentation and to reference to it with its
short name. `DiagnosticVariables` can exist irrespective of the existence of an
actual simulation that is being run. Science packages are encouraged to define
their set of pre-made `DiagnosticVariables`, for example, `ClimaAtmos` comes with
several diagnostics already defined (in the `ALL_DIAGNOSTICS` dictionary).

Let us see how we would define a `DiagnosticVariable`
```julia
import ClimaDiagnostics: DiagnosticVariable

function compute_ta(state, cache, time)
    return state.ta
end

var = DiagnosticVariable(;
    short_name = "ta",
    long_name = "Air Temperature",
    standard_name = "air_temperature",
    comments = "Measured assuming that the air is in quantum equilibrium with the metaverse",
    units = "K",
    compute = compute_ta
)
```

`compute_ta` is the key function here. It determines how the variable should be
computed from the `state`, `cache`, and `time` of the simulation. Typically,
these are packaged within an `integrator` object (e.g., `state = integrator.u`
or `integrator.Y`).

!!! compat "ClimaDiagnostics 0.2.13"

    Support for `compute` was introduced in version `0.2.13`. Prior to this
    version, the in-place `compute!` had to be provided. In this case, `compute`
    has to take an extra argument, `out`. `out` is an area of memory managed by
    `ClimaDiagnostics` that is used to reduce the number of allocations needed
    when working with diagnostics. The first time the diagnostic is called, an
    area of memory is allocated and filled with the value (this is when `out` is
    `nothing`). All the subsequent times, the same space is overwritten, leading
    to much better performance. You should follow this pattern in all your
    diagnostics. This is left to developer to implement, so `compute_ta` would
    look like

    ```julia
    function compute_ta!(out, state, cache, time)
        if isnothing(out)
            return state.ta
        else
            out .= state.ta
        end
    end
    ```

    In general, we do not recommend implementing `compute!`, unless required for
    backward compatibility.

When the expression is anything more complicated than just returning a `Field`,
it is best to return an unevaluated expression represented by a
`Base.Broadcast.Broadcasted` object (such as the ones produced with
[LazyBroadcast.jl](https://github.com/CliMA/LazyBroadcast.jl)). Consider the
following example where we want to shift the temperature to Celsius:
```julia
function compute_ta(state, cache, time)
    return state.ta .- 273.15
end
```

This `compute` function is inefficient because it allocates an entire `Field`
before returning it. Instead, we can return just a recipe on how the diagnostic
should be computed: Using `LazyBroadcast.jl`, the snippet above can be rewritten
as
```julia
import LazyBroadcast: lazy

function compute_ta(state, cache, time)
    return lazy.(state.ta .- 273.15)
end
```
The return value of `compute_ta` is a `Base.Broadcast.Broadcasted` object and
`ClimaDiagnostics` knows how to handle it efficiently, avoiding the intermediate
allocations.

A `DiagnosticVariable` defines what a variable is and how to compute it, but
does not specify when to compute/output it. For that, we need
`ScheduledDiagnostic`s.

#### `ScheduledDiagnostic`s

A `ScheduledDiagnostic` is a `DiagnosticVariable` with attached a schedule on
when it should be computed and output, as well as what reductions should be
performed and how the file should be written.

Continuing our example on `ta`. Suppose we want to compute the average of the
air temperature over a month. We would package this in a `ScheduledDiagnostic`
that knows that we want to compute the air temperature, and we want it averaged
over a month.

Let us examine what is in a `ScheduledDiagnostic` in more details:
- `variable`, the `DiagnosticVariable` we want to compute.
- two `schedule` functions that determine when the variable should be computed
  and output (`compute_schedule_func` and `output_schedule_func`). We have two
  separate entries one for compute and one for output because we might want to
  control them separately. For example, we might want to take the average of
  something every 10 steps, and output it the average every 100 iterations.
  `schedule` functions are powerful, so there is an entire section dedicated to
  them below. `compute_schedule_func` and `output_schedule_func` are likely
  going to be the same unless there are temporal reductions.
- an `output_writer`, an object that knows what to do with the output.
  Examples of writers might be the `DictWriter`, which saves the output to a
  dictionary, or the `NetCDFWriter`, which saves the output to NetCDF files. A
  more complete description of the available writers is in [Saving the
  diagnostics](@ref) page.

- `output_short_name` and `output_long_name`, two strings that specify the names
  that should be used for the output. Typically, `output_short_name` is used for
  file/key names, `output_long_name` is used for descriptive attributes. If none
  is provided, one is automatically generated by the [`output_short_name`](@ref
  ClimaDiagnostics.ScheduledDiagnostics.output_short_name) and
  [`output_long_name`](@ref
  ClimaDiagnostics.ScheduledDiagnostics.output_long_name) functions.
- `reduction_time_func`, a function that implements a temporal reduction.
  Discussed later. This is what you need to implement operations like arithmetic
  averages. A `pre_output_hook!` function can also be passed to do some basic
  normalization operations.

Note that we can have multiple `ScheduledDiagnostic`s for the same
`DiagnosticVariable` (e.g., daily and monthly average temperatures).

##### [`Schedules`](@id schedules_header)

`ScheduledDiagnostic`s contain two arguments `compute_schedule_func` and
`output_schedule_func` which dictate when the variable should be computed and
when it should be output. These objects have to be functions that take a single
argument (the integrator) and return a boolean value.

For example, if we want to call a callback every even step, we could pass
```julia
function compute_every_even(integrator)
    return mod(integrator.step, 2) == 0
end
```
Schedules can be arbitrary. For example, we might want to compute something if
the value of the variable `var` is greater than 100 anywhere. The relevant
schedule for this would be
```julia
function compute_if_larger_than100(integrator)
    return maximum(integrator.u.var) > 100
end
```

Strictly speaking, schedules do not have to be functions, but callable objects. For
example, the `compute_every_even` schedule we defined earlier could be written
for a more general divisor
```julia
struct EveryDivisor
    divisor::Int
end

function (schedule::EveryDivisor)(integrator)
    return mod(integrator.step, schedule.divisor) == 0
end

compute_every_even = EveryDivisor(2)
```
This gives schedules great flexibility because it allows them to contain a state
that can be changed.

`ClimaDiagnostics` define an `AbstractSchedule` type to implement generic
schedules following the pattern just illustrated. One of the main roles of
`AbstractSchedule`s is to have meaningful names that can be used in
files/datasets/error messages, and so on. For this reason, `Schedule`s in
`ClimaDiagnostics` define methods for `short_name` and `long_name`.

If you define your own schedule, you are encouraged to define those methods too.

Let us see a complete example of a new schedule that returns true when a
variable is greater than a threshold.
```julia
import ClimaDiagnostics

struct ExceedThresholdSchedule <: ClimaDiagnostics.AbstractSchedule
    var::Symbol
    threshold::Float64
end

function (schedule::ExceedThresholdSchedule)(integrator)
    return maximum(getproperty(integrator.u, schedule.var)) > schedule.threshold
end

function ClimaDiagnostics.Callback.short_name(schedule::ExceedThresholdSchedule)
    return "$(schedule.var)_more_than_$(schedule.threshold)"
end

function ClimaDiagnostics.Callback.long_name(schedule::ExceedThresholdSchedule)
    return "when max($(schedule.var)) >= $(schedule.threshold)"
end
```
Names are not too important, but they should be meaningful to you.

`ClimaDiagnostics` comes with some predefined schedules for common operations,
such out every N timesteps, or every calendar period. Refer to the
[`Schedules`](@ref schedules_header) section below for more information on what
is already implemented.

!!! note

    `Schedule`s store some information about the last time they were called, so
    different `Schedule`s have to be used and created for different purposes.
    You can use the `deepcopy` function to quickly create a new `Schedule`.

##### Temporal reductions

It is often useful to compute aggregate data (e.g., monthly averages). In
`ClimaDiagnostics`, this is implemented with through temporal reductions.

Let us assume we want to compute the maximum of the air temperature within a
month. To achieve this, we simply pass the `max` function to
`reduction_time_func` and choose our window in the `output_schedule_func`.

The only temporal reductions allowed are ones defined by associative operations,
that is, functions `f` so that `f(a, b, c, d, ...) = f(a, f(b, f(c, f(d,
...))))` (such as the sum). The reason for this restriction comes from the fact
that we do not store all the intermediate values (which would lead to large
consumption of memory). Instead, we accumulate intermediate results. So, the
only statistics that can be computed are the ones that can be computed by adding
one element at the time.

More specifically, when a `ScheduledDiagnostic` is created with a
`reduction_time_func`, `ClimaDiagnostics` allocates an extra area of space
`accumulated` for the accumulated value. Every time `compute_schedule_func` is
true, the `DiagnosticVariable` is computed and saved to `out`. Then,
`accumulated` is updated with the return value of
`reduction_time_func(accumulated, out)`. When `output_schedule_func` is true,
the accumulated value is written with the `writer` and the state reset to the
neutral state.

To allow for greater flexibility, `ClimaDiagnostics` also provides the option to
evaluate a function before the output is saved. This is the `pre_output_hook!`
function that can be provided when defining a `ScheduledDiagnostic`. The
signature for `pre_output_hook!` has to be `pre_output_hook!(accumulated_value,
counter)`, where `counter` is the number of times the diagnostic was called.
Given this, the arithmetic average is obtained with a `+` time reduction and a
`pre_output_hook! = (acc, counter) -> acc .= acc ./ counter`. Given that
averages are very common operations, `ClimaDiagnostics` directly provides the
`pre_output_hook`. So, to define an average, you can directly import and use
`ClimaDiagnostics.average_pre_output_hook!`.

The following is a sketch of what happens at the end of each step for each
`ScheduledDiagnostic`:
```
if compute_schedule_func is true:
    out = compute!
    if reduction_time_func is not nothing:
        accumulated_value = reduction_time_func(accumulated_value, out)
        counter += 1
if output_schedule_func is true:
    pre_output_hook(accumulated_value, counter)
    interpolate(accumulated_value)
    dump(accumulated_value)
    reset(accumulated_value)
    reset(counter)
```

