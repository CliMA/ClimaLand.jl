# [Some notes about the internals of `ClimaDiagnostics`](@id internals_header)

There are multiple moving parts to this package. In this page, we provide some
notes about the internal design. This page also aims at _clarifying the whys_,
i.e., explaining why things are the way they are. Learning about that might help
you extend this package further.

## Schedules

Diagnostics are computed with a callback that is called at the end of each step.
The centerpiece of is the `orchestrate_diagnostics` function, a master callback
that is unconditionally executed at the end of _each_ step.
`orchestrate_diagnostics` loops over each registered diagnostic, computing and
output those for which their trigger condition is met. Conditions are specified
with schedule functions, function that take the `integrator` as single argument
and return a boolean value that determines whether the callback should be
executed or note. Most often, `schedules_func` are not simple functions, but
callable objects (subtypes of `AbstractSchedule`). There are two reasons for
this:
1. Most realistic schedules need to hold additional data (e.g., the last time the function was called)
2. We want to attach names to use in the output

Most of the details regarding schedules are described in the [user guide](@ref user_guide_header).
An internal detail that is not there is related to names. We define a method for
`show` for `AbstractSchedules`. This method calls the `short_name` function.
```julia
Base.show(io::IO, schedule::AbstractSchedule) =  print(io, short_name(schedule))
```
This allows us to set names of `ScheduledDiagnostics` with `"$schedule_func"` in
both the case `schedule_func` is a normal function, or an `AbstractSchedule`.

## Accumulation

Several diagnostics require performing reductions, such as taking the maximum or
the average. Since it is not feasible to store all the lists of all the
intermediate values, we aggregate the results in specific storage areas (e.g.,
we take `max(max(max(max(t1, t2), t3), t4), t5)` instead of `max(t1, t2, t3, t4,
t5)` In this, it is convenient to preallocate the space where we want to
accumulate the intermediate.

Accumulation is accomplished by the `accumulate!` function. All this function
does is applying the binary `reduction_time_func` to the previous accumulated
value and the newly computed one and store the output to the accumulator.

After an accumulated variable is output, the accumulator is reset to its natural
state. This is achieved with the `reset_accumulator!` function. However, we have
to fill the space with something that does not affect the reduction. This, by
definition, is the identity of the operation. The identity of the operation `+`
is `0` because `x + 0 = x` for every `x`.

We have to know the identity for every operation we want to support. Of course,
users are welcome to define their own by adding new methods to
identity_of_reduction.

For instance, to define the identity of the reduction `-`, one would write
```julia
function ClimaDiagnostics.Diagnostics.identity_of_reduction(::typeof(-))
    return 0
end
```
(Or add this to the `reduction_identities.jl` file.)

### On the design of the `DiagnosticsHandler`

There are two possible choices for accumulation of variables: each scheduled
diagnostic can carry its accumulator and counters, or all the accumulators and
counters are managed by a single central handler. `ClimaDiagnostics` implements
this second approach. The author of this package has not decided whether this is
a good idea or not. On one side, this allows us to have a concretely typed and
well defined `DiagnosticsHandler` struct. On the other side, it forces us to
initialize all the diagnostics at the very beginning of the simulation (this can
also be a positive side, because it allows us to compile all the diagnostics at
once). It might be worth exploring the alternative design where the
`ScheduledDiagnostics` get their storage space the first time they are called.

Given this restriction, the main entry point for `ClimaDiagnostics` is the
`IntegratorWithDiagnostics` function. This function is a little dissatisfying
because it creates a new integrator obtained by copying all the fields of the
old one and adding the diagnostics (with
[`Accessors`](https://github.com/JuliaObjects/Accessors.jl)).

## Orchestrate diagnostics

One of the design goals for `orchestrate_diagnostics` is to keep all the
broadcasted expression in the same function scope. This opens a path to optimize
the number of GPU kernel launches. 
