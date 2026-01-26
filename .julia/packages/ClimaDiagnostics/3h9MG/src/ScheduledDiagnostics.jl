module ScheduledDiagnostics

import ..AbstractWriter
import ..Schedules: EveryStepSchedule
import ..DiagnosticVariables:
    DiagnosticVariable, descriptive_short_name, descriptive_long_name

"""
    ScheduledDiagnostic

Conceptually, a `ScheduledDiagnostics` is a `DiagnosticVariable` we want to compute in a given
simulation. For example, it could be the temperature averaged over a day. We can have
multiple `ScheduledDiagnostics` for the same `DiagnosticVariable` (e.g., daily and monthly
average temperatures).
"""
struct ScheduledDiagnostic{
    T1,
    T2,
    OW <: AbstractWriter,
    F1,
    PO,
    DV <: DiagnosticVariable,
}
    """The `DiagnosticVariable` that has to be computed and output"""
    variable::DV

    """A boolean function that determines when this diagnostic should be output. It has to
    take one argument, the integrator. Most typically, only `integrator.t` or
    `integrator.step` are used. Could be a Callback.AbstractSchedule."""
    output_schedule_func::T1

    """Struct that controls out to save the computed diagnostic variable to disk.
    `output_writer` has to implement a method `write_field!` that takes three arguments: the
    value that has to be output, the `ScheduledDiagnostic`, and the integrator. Internally,
    the integrator contains extra information (such as the current timestep). It is
    responsibility of the `output_writer` to properly use the provided information for
    meaningful output."""
    output_writer::OW

    """If not `nothing`, this `ScheduledDiagnostic` receives an area of scratch space `acc`
    where to accumulate partial results. Then, as directed by the `compute_schedule_func`,
    `reduction_time_func` is computed between the previously stored value in `acc` and the
    new value. This implements a running reduction. For example, if `reduction_time_func =
    max`, the space `acc` will hold the running maxima of the diagnostic. To implement
    operations like the arithmetic average, the `reduction_time_func` has to be chosen as
    `sum`, and a `pre_output_hook!` that renormalizes `acc` by the number of samples has to
    be provided. For custom reductions, it is necessary to also specify the identity of
    operation by defining a new method to `identity_of_reduction`."""
    reduction_time_func::F1

    """A boolean function that determines when this diagnostic should be computed. It has to
    take one argument, the integrator. Most typically, only `integrator.t` or
    `integrator.step` are used. Could be a Callback.AbstractSchedule."""
    compute_schedule_func::T2

    # Design note: pre_output_hook!
    #
    # One of our key requirements is to be able to compute arithmetic averages.
    # Unfortunately, computing an arithmetic average requires keeping track of how many
    # elements we are summing up. pre_output_hook! was introduced so that we can think of an
    # average as a sum coupled with division, and perform the division (by the number of
    # elements) before output. pre_output_hook! could be used for other operations, but we
    # decided to keep it simple and target directly the most important use case for us.
    #
    # This choice restricts what reductions can be performed. For example, it is not
    # possible to have a geometric average. If more complex reduction are needed, this
    # mechanism has to be changed.

    """Function that has to be run before saving to disk for reductions (mostly used to
    implement averages). The function `pre_output_hook!` is called with two arguments: the value
    accumulated during the reduction, and the number of times the diagnostic was computed from
    the last time it was output. `pre_output_hook!` should mutate the accumulator in place. The
    return value of `pre_output_hook!` is discarded. An example of `pre_output_hook!` to compute
    the arithmetic average is `pre_output_hook!(acc, N) = @. acc = acc / N`."""
    pre_output_hook!::PO

    """Short name used to output this ScheduledDiagnostic for file names or datasets."""
    output_short_name::String

    """Descriptive name used to output this ScheduledDiagnostic in metadata or attributes."""
    output_long_name::String
end

"""
     ScheduledDiagnostic(; variable::DiagnosticVariable,
                           output_schedule_func,
                           output_writer,
                           reduction_time_func = nothing,
                           reduction_space_func = nothing,
                           compute_schedule_func = isa_reduction ? 1 : output_every,
                           pre_output_hook! = nothing,
                           output_short_name = descriptive_short_name(self),
                           output_short_name = descriptive_long_name(self))


 A `DiagnosticVariable` that has to be computed and output during a simulation with a cadence
 defined by the number of iterations, with an optional reduction applied to it (e.g., compute
 the maximum temperature over the course of every 10 timesteps). This object is turned into
 two callbacks (one for computing and the other for output) and executed by the integrator.

 Keyword arguments
 =================

 - `variable`: The `DiagnosticVariable` that has to be computed and output.

 - `output_every`: Save the results to disk every `output_every` iterations. If `output_every`
                   is non-positive, only output at the first time step.

 - `output_writer`: Function that controls out to save the computed diagnostic variable to
                    disk. `output_writer` has to take three arguments: the value that has to
                    be output, the `ScheduledDiagnostic`, and the integrator. Internally, the
                    integrator contains extra information (such as the current timestep). It
                    is responsibility of the `output_writer` to properly use the provided
                    information for meaningful output.

 - `reduction_time_func`: If not `nothing`, this `ScheduledDiagnostic` receives an area of
                          scratch space `acc` where to accumulate partial results. Then, at
                          every `compute_every`, `reduction_time_func` is computed between
                          the previously stored value in `acc` and the new value. This
                          implements a running reduction. For example, if
                          `reduction_time_func = max`, the space `acc` will hold the running
                          maxima of the diagnostic. To implement operations like the
                          arithmetic average, the `reduction_time_func` has to be chosen as
                          `+`, and a `pre_output_hook!` that renormalize `acc` by the number
                          of samples has to be provided. For custom reductions, it is
                          necessary to also specify the identity of operation by defining a
                          new method to `identity_of_reduction`.

 - `compute_every`: Run the computations every `compute_every` iterations. This is not
                    particularly useful for point-wise diagnostics, where we enforce that
                    `compute_every` = `output_every`. For time reductions, `compute_every` is
                    set to 1 (compute at every timestep) by default. `compute_every` has to
                    evenly divide `output_every`.

 - `pre_output_hook!`: Function that has to be run before saving to disk for reductions
                       (mostly used to implement averages). The function `pre_output_hook!`
                       is called with two arguments: the value accumulated during the
                       reduction, and the number of times the diagnostic was computed from
                       the last time it was output. `pre_output_hook!` should mutate the
                       accumulator in place. The return value of `pre_output_hook!` is
                       discarded. An example of `pre_output_hook!` to compute the arithmetic
                       average is `pre_output_hook!(acc, N) = @. acc = acc / N`.

- `output_short_name`: A descriptive name for this particular diagnostic. If none is
                       provided, one will be generated mixing the short name of the
                       variable, the reduction, and the period of the reduction.
                       Normally, it has to be unique. In `ClimaAtmos`, we follow the CMIP
                       conventions for this.

- `output_long_name`: A descriptive name for this particular diagnostic. If none is
                      provided, one will be generated mixing the short name of the
                      variable, the reduction, and the period of the reduction.

 """
function ScheduledDiagnostic(;
    variable::DiagnosticVariable,
    output_writer,
    reduction_time_func = nothing,
    compute_schedule_func = EveryStepSchedule(),
    output_schedule_func = isnothing(reduction_time_func) ?
                           deepcopy(compute_schedule_func) :
                           EveryStepSchedule(),
    pre_output_hook! = (accum, count) -> nothing,
    output_short_name = descriptive_short_name(
        variable,
        output_schedule_func,
        reduction_time_func,
        pre_output_hook!,
    ),
    output_long_name = descriptive_long_name(
        variable,
        output_schedule_func,
        reduction_time_func,
        pre_output_hook!,
    ),
)
    # pre_output_hook! has to be a function, but it is much more intuitive to specify
    # `nothing` when we want nothing to happen. Here, we convert the nothing keyword
    # into a function that does nothing
    if isnothing(pre_output_hook!)
        pre_output_hook! = (accum, count) -> nothing
    end

    T = typeof(variable)
    T1 = typeof(output_schedule_func)
    T2 = typeof(compute_schedule_func)
    OW = typeof(output_writer)
    F1 = typeof(reduction_time_func)
    PO = typeof(pre_output_hook!)

    ScheduledDiagnostic{T1, T2, OW, F1, PO, T}(
        variable,
        output_schedule_func,
        output_writer,
        reduction_time_func,
        compute_schedule_func,
        pre_output_hook!,
        output_short_name,
        output_long_name,
    )
end

"""
    output_short_name(sd::ScheduledDiagnostic)

Return the short name to use for output of the `ScheduledDiagnostic` `sd`.
"""
function output_short_name(sd::ScheduledDiagnostic)
    return sd.output_short_name
end

"""
    output_long_name(sd::ScheduledDiagnostic)

Return the long name to use for output of the `ScheduledDiagnostic` `sd`.
"""
function output_long_name(sd::ScheduledDiagnostic)
    return sd.output_long_name
end

function Base.:(==)(sd1::T, sd2::T) where {T <: ScheduledDiagnostic}
    # We provide == because we don't want to compare with === because we have
    # RefValues
    return all(
        getproperty(sd1, p) == getproperty(sd2, p) for p in propertynames(sd1)
    )
end


end
