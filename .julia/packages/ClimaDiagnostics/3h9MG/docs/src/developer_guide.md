# I am a developer, how do I add `ClimaDiagnostics.jl` to my package?

This page provides additional documentation on abstractions to use
`ClimaDiagnostics`. Before reading this page, make sure you are familiar with
the terminology. You need to know what a [`DiagnosticVariable`](@ref ClimaDiagnostics.DiagnosticVariables.DiagnosticVariable) and a
[`ScheduledDiagnostic`](@ref ClimaDiagnostics.ScheduledDiagnostics.ScheduledDiagnostic) are.

There are two components needed to add support for `ClimaDiagnostics.jl` in your package.

1. A way to convert users' intentions to a list of [`ScheduledDiagnostic`](@ref
   ClimaDiagnostics.ScheduledDiagnostics.ScheduledDiagnostic)
2. A call to [`IntegratorWithDiagnostics`](@ref
   ClimaDiagnostics.IntegratorWithDiagnostics), if you have an integrator (if
   you do not have an integrator read the [What if I do not have an
   integrator?](@ref) section)

## Step 2

Let us assume that `scheduled_diagnostics` is the list of `ScheduledDiagnostic`s
obtained from step 1. (more on this later), and `integrator` a `SciML`
integrator.

All we need to do to add diagnostics is
```julia
import ClimaDiagnostics: IntegratorWithDiagnostics

integrator = IntegratorWithDiagnostics(integrator, scheduled_diagnostics)
```
Creating an `IntegratorWithDiagnostics` results in calling all the diagnostics
once. Therefore, the compile and runtime of this function can be significant if
you have a large number of diagnostics.

You can learn about what is happening under the hook in the [Internals](@ref internals_header)
page.

This is pretty much all that you need to know about step 2.

### What if I do not have an integrator?

`ClimaDiagnostics` does not assume that you are working with an integrator
object, so you can use it in your project even when you hold your simulation
fields in other data structures. 

`DiagnosticVariable`s in `ClimaDiagnostics` come with a `compute!` function with
a specific signature `compute!(out, state, cache, time)`. This is required
because `ClimaDiagnostics` calls them by calling `compute!(out, integrator.u,
integrator.p, integrator.t)`. If you do not have an integrator, all you have to
do is create a wrapper that exposes this interface.

For example, you could have `fields`, `parameters`, and `time`, and you could
wrap them in a `my_integrator = (; u = fields, p = parameters, t = time)` and
use the objects in the `compute!` function as you like. 

The integrator-equivalent object is also used in [`Schedules`](@ref
schedules_header). For time-dependent schedules, `integrator.t` is used. For
step-dependent schedules (e.g.,
[`ClimaDiagnostics.Schedules.EveryStepSchedule`](@ref)), the `step` field is
used. If you are using these schedules, do also add `step` to your
`my_integrator`.

If you do not have an integrator, it is likely that you do not have a callback
system either. In this case, you would have to manually call the
[`ClimaDiagnostics.orchestrate_diagnostics`](@ref) at the end of each of your
steps. `orchestrate_diagnostics` is the heart of `ClimaDiagnostics`: it runs all
the calculations are saves the output. Usually, it is added as a callback called
at the end of every step. If this is not possible, you should manually call it.

In this case, your interface would probably require you to manually construct
the [`ClimaDiagnostics.DiagnosticsHandler`](@ref)

!!! summary
    
    - Wrap your objects in `my_integrator = (; u = fields, p = parameters, t =
      time, step = step)`
    - Manually construct a [`ClimaDiagnostics.DiagnosticsHandler`](@ref) from
      your [`ClimaDiagnostics.ScheduledDiagnostic`](@ref)
    - Manually call a [`ClimaDiagnostics.orchestrate_diagnostics`](@ref) at the end of each step

## Step 1

Step 1 in the recipe to bring `ClimaDiagnostics` to your package strongly
depends on you.

In this section, I will present a tower of interfaces that you can put in place
to make it more convenient for your users. Keep in mind that each layer trades
more convenience for less flexibility. So, as you set up your interfaces, I
recommend you keep them exposed so that your users can access lower-level
functions if they need to.

### Level 0: do nothing

At the zero-th level, you let your users work directly with `ClimaDiagnostics`.
This means that they will have to define their own `DiagnosticVariable`s and
`ScheduledDiagnostic`s. This also requires that your simulation is executed as a
julia script.

It is a good idea for your users to be aware of this possibility because it
brings enormous power. `ScheduledDiagnostic`s can be triggered on arbitrary
conditions, and your users could be creative with that. For example, users might
want to compute and output a variable `var1` when they find that the maximum of
variable `var2` is greater than a threshold (e.g., for debugging).

Let us see the simplest example to accomplish this
```julia
import ClimaDiagnostics: DiagnosticVariable, ScheduledDiagnostic
import ClimaDiagnostics.Writers: DictWriter

myvar = DiagnosticVariable(; compute = (u, p, t) -> u.var1)

myschedule = (integrator) -> maximum(integrator.u.var2) > 10.0

diag = ScheduledDiagnostic(variable = myvar,
                           compute_schedule_func = myschedule,
                           output_schedule_func = myschedule,
                           output_writer = DictWriter())
```

Now we can go to step 2 and 3 in the previous list and pass `[diag]` to the
`DiagnosticsHandler`.

Point your users to the documentation of this package for them to learn how to
use it in its full power.

### Level 1: provide a database of `DiagnosticVariable`s

As a package developer, you know that there is a large collection of variables
that several users will be interested in. For example, if you are running an
atmospheric simulation, your users will want to be able to look at the air
temperature. For this reason, it is a very good (and user-friendly) idea to
provide a collection of `DiagnosticVariable`s ready to be used. In this section,
I sketch how you could go about and implement this.

Your `DiagnosticVariable`s database can be represented as a dictionary
`ALL_DIAGNOSTICS` indexed over the short name of the variable. Then, you could
provide adders and accessors.

This might look like the following:
```julia
module Diagnostics
import ClimaDiagnostics: DiagnosticVariable

const ALL_DIAGNOSTICS = Dict{String, DiagnosticVariable}()

"""

    add_diagnostic_variable!(; short_name,
                               long_name,
                               standard_name,
                               units,
                               description,
                               compute)


Add a new variable to the `ALL_DIAGNOSTICS` dictionary (this function mutates the state of
`ALL_DIAGNOSTICS`).

If possible, please follow the naming scheme outline in
https://airtable.com/appYNLuWqAgzLbhSq/shrKcLEdssxb8Yvcp/tblL7dJkC3vl5zQLb

Keyword arguments
=================

- `short_name`: Name used to identify the variable in the output files and in the file
                names. Short but descriptive. Diagnostics are identified by the short name.

- `long_name`: Name used to identify the variable in the output files.

- `standard_name`: Standard name, as in
                   http://cfconventions.org/Data/cf-standard-names/71/build/cf-standard-name-table.html

- `units`: Physical units of the variable.

- `comments`: More verbose explanation of what the variable is, or comments related to how
              it is defined or computed.

- `compute`: Function that computes the diagnostic variable from the state, cache, and time. The function
             should return a `Field` or a `Base.Broadcast.Broadcasted` expression. It should not allocate
             new `Field`: if you find yourself using a dot, that is a good indication you should be using
             `lazy`.
"""
function add_diagnostic_variable!(;
    short_name,
    long_name,
    standard_name = "",
    units,
    comments = "",
    compute,
)
    haskey(ALL_DIAGNOSTICS, short_name) && @warn(
        "overwriting diagnostic `$short_name` entry containing fields\n" *
        "$(map(
            field -> "$(getfield(ALL_DIAGNOSTICS[short_name], field))",
            # We cannot really compare functions...
            filter(field -> !(field in (:compute!, :compute)), fieldnames(DiagnosticVariable)),
        ))"
    )

    ALL_DIAGNOSTICS[short_name] = DiagnosticVariable(;
        short_name,
        long_name,
        standard_name,
        units,
        comments,
        compute,
    )

"""
    get_diagnostic_variable!(short_name)

Return a `DiagnosticVariable` from its `short_name`, if it exists.
"""
function get_diagnostic_variable(short_name)
    haskey(ALL_DIAGNOSTICS, short_name) ||
        error("diagnostic $short_name does not exist")

    return ALL_DIAGNOSTICS[short_name]
end

end
```
Of course, you should have the fields and comments that are relevant to your package.

Next, as a developer, you will use `add_diagnostic_variable!` to populate your
database. You can also expose your users to this function so that they can
extend their personal database in their simulations.

A simple example of a new variable might look like
```julia
###
# Density (3d)
###
add_diagnostic_variable!(
    short_name = "rhoa",
    long_name = "Air Density",
    standard_name = "air_density",
    units = "kg m^-3",
    compute = (state, cache, time) -> state.c.ρ,
)
```

When writing compute functions, make them lazy with
[LazyBroadcast.jl](https://github.com/CliMA/LazyBroadcast.jl) to improve
performance and avoid intermediate allocations. To do that, add `LazyBroadcast`
to your dependencies and import `lazy`. A slight variation of the previous
example would look like

```julia
###
# Density (3d)
###
add_diagnostic_variable!(
    short_name = "rhoa",
    long_name = "Air Density",
    standard_name = "air_density",
    units = "kg m^-3",
    compute = (state, cache, time) -> lazy.(1000 .* state.c.ρ),
)
```
Where we added the `1000` to simulate a more complex expression. If you didn't have
`lazy`, the diagnostic would allocate an intermediate `Field`, severly hurting performance.

It is a good idea to put safeguards in place to ensure that your users will not
be allowed to call diagnostics that do not make sense for the simulation they
are running. If your package has a notion of `Model` that is stored in `p`, you
can dispatch over that and return an error. A simple example might be
```julia
###
# Specific Humidity
###
compute_hus(state, cache, time) =
    compute_hus(state, cache, time, cache.atmos.moisture_model)

compute_hus(state, cache, time) =
    compute_hus!(state, cache, time, cache.model.moisture_model)
compute_hus(_, _, _, model::T) where {T} =
    error("Cannot compute hus with $model")

function compute_hus(
    state,
    cache,
    time,
    moisture_model::T,
) where {T <: Union{EquilMoistModel, NonEquilMoistModel}}
    return lazy.(state.c.ρq_tot ./ state.c.ρ)
end

add_diagnostic_variable!(
    short_name = "hus",
    long_name = "Specific Humidity",
    standard_name = "specific_humidity",
    units = "kg kg^-1",
    comments = "Mass of all water phases per mass of air",
    compute = compute_hus,
)
```
This relies on dispatching over `moisture_model`. If `model` is not in
`Union{EquilMoistModel, NonEquilMoistModel}`, the code returns an informative error.

If you provide a database, users can create their `ScheduledDiagnostic`s
directly from the `DiagnosticVariable`s you provided.

For instance to output the specific humidity every 5 iterations:
```julia
import ClimaDiagnostics: ScheduledDiagnostic
import ClimaDiagnostics.Callbacks: DivisorSchedule
import ClimaDiagnostics.Writers: DictWriter

diag = ScheduledDiagnostic(variable = get_diagnostic_variable!("hus"),
                           output_schedule_func = DivisorSchedule(5),
                           output_writer = DictWriter())
```

Alongside with providing the `DiagnosticVariable`s, you can also provide
convenience functions for standard operations.

For example, you could provide
```julia
using ClimaDiagnostics.Callbacks: EveryStepSchedule, EveryDtSchedule

function monthly_average(FT, short_name; output_writer)
    period = 30 * 24 * 60 * 60 * one(FT)
    return ScheduledDiagnostic(
            variable = get_diagnostic_variable(short_name),
            compute_schedule_func = EveryStepSchedule(),
            output_schedule_func = EveryDtSchedule(period),
            reduction_time_func = (+),
            output_writer = output_writer,
            pre_output_hook! = average_pre_output_hook!,
        )
end
```
Allowing users to just call `monthly_average(Float32, "hus", writer)`.

> Note: `ClimaDiagnostics` will probably provided these schedules natively at
> some point in the future.

### Level 2: Provide higher-level interfaces (e.g., YAML)

Finally, you can set in place that parses user input (e.g., from command line or
text files) into `ScheduledDiagnostics` using the short names in your database.
Of course, this interface will be limited to what you expose.

For example, a simple parser that allow users to specify `ScheduledDiagnostics`
by their short name, accumulation/output period, and their writer might look
like the following:
```julia
import ClimaDiagnostics: average_pre_output_hook!, HDF5Writer, NetCDFWriter, ScheduledDiagnostic

function parse_yaml(parsed_args, source_space)
    # We either get the diagnostics section in the YAML file, or we return an empty list
    # (which will result in an empty list being created by the map below)
    yaml_diagnostics = get(parsed_args, "diagnostics", [])

    # ALLOWED_REDUCTIONS is the collection of reductions we support. The keys are the
    # strings that have to be provided in the YAML file. The values are tuples with the
    # function that has to be passed to reduction_time_func and the one that has to passed
    # to pre_output_hook!

    # We make "nothing" a string so that we can accept also the word "nothing", in addition
    # to the absence of the value
    #
    # NOTE: Everything has to be lowercase in ALLOWED_REDUCTIONS (so that we can match
    # "max" and "Max")
    ALLOWED_REDUCTIONS = Dict(
        "nothing" => (nothing, nothing), # nothing is: just dump the variable
        "max" => (max, nothing),
        "min" => (min, nothing),
        "average" => ((+), average_pre_output_hook!),
    )

    output_dir = parsed_args.output_dir

    hdf5_writer = HDF5Writer(output_dir)
    netcdf_writer = CAD.NetCDFWriter(
        source_space,
        output_dir,
    )
    writers = (hdf5_writer, netcdf_writer)

    # The default writer is HDF5
    ALLOWED_WRITERS = Dict(
        "nothing" => netcdf_writer,
        "h5" => hdf5_writer,
        "hdf5" => hdf5_writer,
        "nc" => netcdf_writer,
        "netcdf" => netcdf_writer,
    )

    diagnostics_ragged = map(yaml_diagnostics) do yaml_diag
        short_names = yaml_diag["short_name"]
        output_name = get(yaml_diag, "output_name", nothing)

        map(short_names) do short_name
            # Return "nothing" if "reduction_time" is not in the YAML block
            #
            # We also normalize everything to lowercase, so that can accept "max" but
            # also "Max"
            reduction_time_yaml =
                lowercase(get(yaml_diag, "reduction_time", "nothing"))

            if !haskey(ALLOWED_REDUCTIONS, reduction_time_yaml)
                error("reduction $reduction_time_yaml not implemented")
            else
                reduction_time_func, pre_output_hook! =
                    ALLOWED_REDUCTIONS[reduction_time_yaml]
            end

            writer_ext = lowercase(get(yaml_diag, "writer", "nothing"))

            if !haskey(ALLOWED_WRITERS, writer_ext)
                error("writer $writer_ext not implemented")
            else
                writer = ALLOWED_WRITERS[writer_ext]
            end

            haskey(yaml_diag, "period") ||
                error("period keyword required for diagnostics")

            period_seconds = FT(time_to_seconds(yaml_diag["period"]))

            if isnothing(reduction_time_func)
                compute_every = CAD.EveryDtSchedule(period_seconds)
            else
                compute_every = CAD.EveryStepSchedule()
            end

            ScheduledDiagnostic(
                variable = get_diagnostic_variable(short_name),
                output_schedule_func = CAD.EveryDtSchedule(period_seconds),
                compute_schedule_func = compute_every,
                reduction_time_func = reduction_time_func,
                pre_output_hook! = pre_output_hook!,
                output_writer = writer,
            )
        end
    end

    # Flatten the array of arrays of diagnostics
    diagnostics = vcat(diagnostics_ragged...)
end
```
This will be controlled by YAML blocks like
```yaml
diagnostics:
    - short_name: ["ta", "va"]
      period: 60s
      writer: nc
    - short_name: ["ua"]
      period: 1200s
      reduction_time: "average"
```
It is typically a good idea to add the default diagnostics to the set of
YAML-specified ones.

