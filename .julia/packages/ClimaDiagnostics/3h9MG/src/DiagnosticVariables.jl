module DiagnosticVariables

import ..Schedules: AbstractSchedule, long_name

"""
    DiagnosticVariable(;
        compute!,
        short_name = "",
        long_name = "",
        standard_name = "",
        units = "",
        comments = ""
    )

A recipe to compute a diagnostic variable from the state, along with some useful metadata.

The primary use for `DiagnosticVariable`s is to be embedded in a `ScheduledDiagnostic` to
compute diagnostics while the simulation is running.

The metadata is used exclusively by the `output_writer` in the `ScheduledDiagnostic`. It is
responsibility of the `output_writer` to follow the conventions about the meaning of the
metadata and their use.

In `ClimaAtmos`, we roughly follow the naming conventions listed in this file:
https://airtable.com/appYNLuWqAgzLbhSq/shrKcLEdssxb8Yvcp/tblL7dJkC3vl5zQLb

!!! compat "ClimaDiagnostics 0.2.13"

    Support for `compute` was introduced in version `0.2.13`. Prior to this version, the
    in-place `compute!` had to be provided.

Keyword arguments
=================

- `compute!`: Function that computes the diagnostic variable from the `state`, `cache`, and
              `time`. In addition to these three arguments, `compute!` has to take four
              arguments, the destination where to write the result of the computation.

- `compute`: Function that computes the diagnostic variable from the `state`, `cache`, and
             `time` and returns either a `Field` or an un-evaluated expression (e.g., with
             `LazyBroadcast.lazy`).

- `short_name`: Name used to identify the variable in the output files and in the file
                names. Short but descriptive. `ClimaAtmos` follows the CMIP conventions and
                the diagnostics are identified by the short name.

- `long_name`: Name used to describe the variable in the output files.

- `standard_name`: Standard name, as in
                   http://cfconventions.org/Data/cf-standard-names/71/build/cf-standard-name-table.html

- `units`: Physical units of the variable.

- `comments`: More verbose explanation of what the variable is, or comments related to how
              it is defined or computed.
"""
struct DiagnosticVariable{
    F1 <: Union{Function, Nothing},
    F2 <: Union{Function, Nothing},
}
    compute!::F1
    compute::F2
    short_name::String
    long_name::String
    standard_name::String
    units::String
    comments::String
end

function DiagnosticVariable(;
    compute = nothing,
    compute! = nothing,
    short_name = "",
    long_name = "",
    standard_name = "",
    units = "",
    comments = "",
)
    only_one_compute_provided = xor(isnothing(compute), isnothing(compute!))

    if !only_one_compute_provided
        error(
            "One and only one between `compute` and `compute!` has to be provided",
        )
    end

    return DiagnosticVariable(
        compute!,
        compute,
        short_name,
        long_name,
        standard_name,
        units,
        comments,
    )
end

"""
    short_name(dv::DiagnosticVariable)

Return the short name associated to the given `DiagnosticVariable`.
"""
function short_name(dv::DiagnosticVariable)
    return dv.short_name
end

"""
    long_name(dv::DiagnosticVariable)

Return the long name associated to the given `DiagnosticVariable`.
"""
function long_name(dv::DiagnosticVariable)
    return dv.long_name
end

"""
    average_pre_output_hook!

Function to use as `pre_output_hook!` for a `ScheduledDiagnostic` to compute an arithmetic average.
"""
function average_pre_output_hook!(accum, counter)
    @. accum = accum / counter
    return nothing
end

"""
    descriptive_short_name(variable::DiagnosticVariable,
                           output_schedule_func,
                           reduction_time_func,
                           pre_output_hook!)

Return a compact, unique-ish, identifier generated from the given information. This function
is useful for filenames and error messages.
"""
function descriptive_short_name(
    variable::DiagnosticVariable,
    output_schedule_func,
    reduction_time_func,
    pre_output_hook!;
)
    var = "$(variable.short_name)"
    isa_reduction = !isnothing(reduction_time_func)

    if isa_reduction
        red = "$(reduction_time_func)"

        # Let's check if we are computing the average. Note that this might slip under the
        # radar if the user passes their own pre_output_hook!.
        if reduction_time_func == (+) &&
           nameof(pre_output_hook!) == :average_pre_output_hook!
            red = "average"
        end
        suffix = "$red"
    else
        suffix = "inst"
    end
    return "$(var)_$(output_schedule_func)_$(suffix)"
end

"""
    descriptive_long_name(variable::DiagnosticVariable,
                          output_every,
                          reduction_time_func,
                          pre_output_hook!)

Return a verbose description of the given output variable.

This function is useful for attributes in output files.
"""
function descriptive_long_name(
    variable::DiagnosticVariable,
    output_schedule_func,
    reduction_time_func,
    pre_output_hook!;
)
    var = "$(variable.long_name)"
    isa_reduction = !isnothing(reduction_time_func)

    if isa_reduction
        red = "$(reduction_time_func)"

        # Let's check if we are computing the average. Note that this might slip under the
        # radar if the user passes their own pre_output_hook!.
        if reduction_time_func == (+) &&
           pre_output_hook! == average_pre_output_hook!
            red = "average"
        end

        if output_schedule_func isa AbstractSchedule
            suffix = "$(red) within $(long_name(output_schedule_func))"
        else
            suffix = red
        end
    else
        suffix = "Instantaneous"
    end
    return "$(var), $(suffix)"
end
end
