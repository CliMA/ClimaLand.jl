# `SpaceVaringInputs` and `TimeVaryingInputs`

Most models require external inputs to work. Examples of inputs are an analytic
function that prescribes the sea-surface temperature in time, or a file that
describes the types of plants on the surface of the globe. The
`SpaceVaringInputs` and `TimeVaryingInputs` modules provide a unified
infrastructure to handle all these cases.

## [`TimeVaryingInputs`](@id timevaryinginput)

> This extension is loaded when loading `ClimaCore` is loaded. In addition to
> this, if NetCDF files are used, `NCDatasets` has to be loaded too. Finally, a
> `Regridder` is needed (which might require importing additional packages).

A `TimeVaryingInput` is an object that knows how to fill a `ClimaCore` `Field`
at a given simulation time `t`. `TimeVaryingInputs` can be constructed in a
variety of ways, from using analytic functions, to NetCDF data. They expose one
interface, `evaluate!(dest_field, tv, time)`, which can be used by model
developers to update their `Field`s.

This example shows that `TimeVaryingInput` can take different types of inputs
and be used with a single interface (`evaluate!`). In all of this,
`TimeVaryingInput`s internally handle all the complexity related to reading
files (using [`FileReaders`](@ref file_reader_module)), dealing with parallelism and GPUs,
regridding onto the computational domains (using [`Regridders`](@ref regridder_module) and
[`DataHandling`](@ref datahandling_module)), and so on.

`TimeVaryingInputs` support:
- analytic functions of time;
- pairs of 1D arrays (e.g., for `PointSpaces` or constant fields);
- 2/3D NetCDF files (including composing multiple variables from one or more files into one variable);
- linear interpolation in time (default), nearest neighbors, and "period filling";
- boundary conditions and repeating periodic data.

It is possible to pass down keyword arguments to underlying constructors in the
`Regridder` with the `regridder_kwargs` and `file_reader_kwargs`. These have to
be a named tuple or a dictionary that maps `Symbol`s to values.

### NetCDF file inputs
2D or 3D NetCDF files can be provided as inputs using `TimeVaryingInputs`. This
could be a single variable provided in a single file, multiple variables provided
in a single file, or multiple variables each coming from a unique file.
When using multiple variables, a composing function must be provided as well,
which will be used to combine the input variables into one data variable that
is ultimately stored in the `TimeVaryingInput`. In this case, the order of
variables provided in `varnames` determines the order of the arguments
passed to the composing function.

Note that if a non-identity pre-processing function is provided as part of
`file_reader_kwargs`, it will be applied to each input variable before they
are composed.
All input variables to be composed together must have the same spatial and
temporal dimensions.

Composing multiple input variables is currently only supported with the
[`InterpolationsRegridder`](@ref interp_regridder), not with [`TempestRegridder`](@ref tempest_regridder). The regridding
is applied after the pre-processing and composing.

Composing multiple input variables in one `Input` is also possible with
a `SpaceVaryingInput`, and everything mentioned here applies in that case.

#### Example: NetCDF file input with multiple input variables

Suppose that the input NetCDF file `era5_example.nc` contains two variables `u`
and `v`, and we care about their sum `u + v` but not their individual values.
We can provide a pointwise composing function to perform the sum, along with
the [`InterpolationsRegridder`](@ref interp_regridder) to produce the data we want, `u + v`.
The `preprocess_func` passed in `file_reader_kwargs` will be applied to `u`
and to `v` individually, before the composing function is applied. The regridding
is applied after the composing function. `u` and `v` could also come from separate
NetCDF files, but they must still have the same spatial and temporal dimensions.

```julia
# Define the pointwise composing function we want, a simple sum in this case
compose_function = (x, y) -> x + y
# Define pre-processing function to convert units of input
unit_conversion_func = (data) -> 1000 * data

timevaryinginput = TimeVaryingInputs.TimeVaryingInput("era5_example.nc",
                                        ["u", "v"],
                                        target_space,
                                        start_date = Dates.DateTime(2000, 1, 1),
                                        regridder_type = :InterpolationsRegridder,
                                        file_reader_kwargs = (; preprocess_func = unit_conversion_func),
                                        compose_function)
```

The same arguments (excluding `start_date`) could be passed to a
`SpaceVaryingInput` to compose multiple input variables with that type.

#### Example: Data split across multiple NetCDF files

Often, large datasets come chunked, meaning that the data is split across
multiple files with each file containing only a subset of the time interval.
`TimeVaryingInput`s know to combine data across multiple files as it were
provided in a single file. To do use this feature, just pass the list of file
paths. While it is not required for the files to be in order, it is good
practice to pass them in ascending order by time.

For example:
```julia
timevaryinginput = TimeVaryingInputs.TimeVaryingInput(["era5_1980.nc", "era5_1981.nc"],
                                                       "u",
                                                       target_space,
                                                       start_date = Dates.DateTime(1980, 1, 1),
                                                       regridder_type = :InterpolationsRegridder
                                                       )
```

This capability is only available for the `InterpolationsRegridder`.

Read more about this feature in the page about [`DataHandler`](@ref datahandling_module).

### Extrapolation boundary conditions

`TimeVaryingInput`s can have multiple boundary conditions for extrapolation. By
default, the `Throw` condition is used, meaning that interpolating onto a point
that is outside the range of definition of the data is not allowed. Other
boundary conditions are allowed. With the `Flat` boundary condition, when
interpolating outside of the range of definition, return the value of the of
closest boundary is used instead.

To set these boundary conditions, construct the relevant method passing the
argument. For example, to combine `NearestNeighbor` with `Flat`:
```julia
import ClimaUtilities: TimeVaryingInputs

method = TimeVaryingInputs.NearestNeighbor(TimeVaryingInputs.Flat())
```

A boundary condition that is often useful is `PeriodicCalendar`, which repeats
the data over and over.

In general `PeriodicCalendar` takes two inputs: the `period` and `repeat_date`.
The repeat period is a `Dates.DatePeriod` (e.g., `Dates.Year(1)`) that defines
the duration of the period that has to be repeated. The `repeat_date` defines
what date range needs to be repeated. For example, if `period = Dates.Month(1)`
and `repeat_date = Dates.Date(1993, 11)`, November 1993 will be repeated.

The two inputs are not required. When they are not provided, `ClimaUtilities`
will assume that the input data constitutes one period and use that. For
example, if the data is defined from `t0` to `t1` (e.g., 1 and 5), interpolating
over `t > t1` (e.g., 7) is equivalent to interpolating to `t*` where `t*` is the
modulus of `t` and the range (3 in this case). In this case, `PeriodicCalendar`
requires the data to be uniformly spaced in time. To enable this boundary
condition, pass `LinearInterpolation(PeriodicCalendar())` to the
`TimeVaryingInput` (or `NearestNeighbor(PeriodicCalendar())`).

!!! note

    This `PeriodicCalendar` is different from what you might be used to,
    where the identification is `t1 = t0`. Here, we identify `t1 + dt = t0`.
    This is so that we can use it to repeat calendar data.

### `LinearPeriodFillingInterpolation`

Often, data is not available at the frequency we would like it to be. For
example, we might have hourly data for a given quantity but only on the 15th of
the month. Performing linear interpolation with data with this type of gap is
typically not accurate. Consider the example of a quantity with a diurnal cycle
but measured only once a month. If we were to blindly perform linear
interpolation, we would find that the diurnal cycle is completely removed for
every day of the month but the 15th. This is because we would interpolate the
last point for the day of a given month, with the first for the following.

`LinearPeriodFillingInterpolation` is an interpolation method that solves this
problem by preserving periodic structures. This is accomplished by performing
linear interpolation across corresponding periods (in the case of the day,
across corresponding hours of different days). For more information, please
refer to the docstring.

### Example

Let `target_space` be the computational domain (a `ClimaCore` `Space`) and
`cesm_albedo.nc` a NetCDF file containing albedo data as a function of time in a
variable named `alb`.

```julia
import ClimaUtilities: TimeVaryingInputs
import ClimaCore
import NCDatasets
import ClimaCoreTempestRemap
# Loading ClimaCore, NCDatasets, ClimaCoreTempestRemap loads the extensions we need

function evolve_model(albedo_tv, albedo_field)
    new_t = t + dt
    # First, we update the albedo to the new time
    evaluate!(albedo_field, albedo_tv, new_t)
    # Now we can do all the operations we want we albedo_filed
    # rhs = ...
end

# Let us prepare an empty Field that will contain the albedo
albedo_field = zero(target_space)

# If the albedo is an analytic function of time
albedo_tv_an = TimeVaryingInput((t) -> 0.5)

# If the albedo comes from data

# start_date is the calendar date at the beginning of our simulation
start_date = Dates.DateTime(2000, 1, 1)
albedo_tv = TimeVaryingInputs.TimeVaryingInput("cesem_albedo.nc", "alb", target_space;
                                               start_date, regridder_kwargs = (; regrid_dir = "/tmp"))
# When using data from files, the data is automatically interpolated to the correct
# time

# In either cases, we can always call evolve_model(albedo_tv, albedo_field), so
# model developers do not have to worry about anything :)
```

As seen in this example, `Inputs` can take keyword arguments and pass them down
to other constructors. This often used to preprocess files that are being read
(most commonly to change units). For example, if we want to multiply the albedo
by a factor of 100, we would change `albedo_tv` with
```julia
albedo_tv = TimeVaryingInputs.TimeVaryingInput("cesem_albedo.nc", "alb", target_space;
                                               start_date, regridder_kwargs = (; regrid_dir = "/tmp"),
                                               file_reader_kwargs = (; preprocess_func = (x) -> 100x))
```

!!! note

    In this example we used the [`TempestRegridder`](@ref tempest_regridder).
    This is not the best
    choice in most cases because the [`TempestRegridder`](@ref tempest_regridder) is slower, and
    not well-compatible with MPI and GPUs (`ClimaUtilities` implements
    workarounds for this, so the code would still work).
    [`InterpolationsRegridder`](@ref interp_regridder) should be preferred, unless there is a
    strict requirement of conservation: while [`TempestRegridder`](@ref tempest_regridder) is
    guaranteed to conserve various properties, [`InterpolationsRegridder`](@ref interp_regridder)
    is not.

## [`SpaceVaryingInputs`](@id spacevaryinginput)

> This extension is loaded when loading `ClimaCore` is loaded. In addition to
> this, if NetCDF files are used, `NCDatasets` has to be loaded too. Finally, a
> `Regridder` is needed (which might require importing additional packages).

`SpaceVaryingInput`s uses the same building blocks as `TimeVaryingInput`
(chiefly the [`DataHandling`](@ref datahandling_module) datahandling_module) to construct a `Field` from
different sources.

`SpaceVaryingInputs` support:
- analytic functions of coordinates;
- pairs of 1D arrays (for columns);
- 2/3D NetCDF files (including composing multiple variables from one or more files into one variable).

In some ways, a `SpaceVaryingInput` can be thought as an alternative constructor
for a `ClimaCore` `Field`.

It is possible to pass down keyword arguments to underlying constructors in the
`Regridder` with the `regridder_kwargs` and `file_reader_kwargs`. These have to
be a named tuple or a dictionary that maps `Symbol`s to values.

`SpaceVaryingInputs` support reading individual input variables from NetCDF files,
as well as composing multiple input variables into one `SpaceVaryingInput`.
See the [`TimeVaryingInput`](@ref timevaryinginput) "NetCDF file inputs" section for more
information about this feature.

### Example

Let `target_space` be a `ClimaCore` `Space` where we want the `Field` to be
defined on and `cesm_albedo.nc` a NetCDF file containing albedo data as a time
in a variable named `alb`.

```julia
import ClimaUtilities: SpaceVaryingInputs
import ClimaCore
import NCDatasets
import ClimaCoreTempestRemap
# Loading ClimaCore, NCDatasets, ClimaCoreTempestRemap loads the extensions we need

# Albedo as an analytic function of lat and lon
albedo_latlon_fun = (coord) -> 0.5 * coord.long * coord.lat

albedo = SpaceVaryingInputs.SpaceVaryingInput(albedo_latlon_fun, target_space)

albedo_from_file = SpaceVaryingInputs.SpaceVaryingInput("cesm_albedo.nc", "alb", target_space, regridder_kwargs = (; regrid_dir = "/tmp"))
```

## API

```@docs
ClimaUtilities.SpaceVaryingInputs.SpaceVaryingInput
ClimaUtilities.TimeVaryingInputs.AbstractInterpolationMethod
ClimaUtilities.TimeVaryingInputs.AbstractInterpolationBoundaryMethod
ClimaUtilities.TimeVaryingInputs.NearestNeighbor
ClimaUtilities.TimeVaryingInputs.LinearInterpolation
ClimaUtilities.TimeVaryingInputs.Throw
ClimaUtilities.TimeVaryingInputs.PeriodicCalendar
ClimaUtilities.TimeVaryingInputs.Flat
ClimaUtilities.TimeVaryingInputs.evaluate!
ClimaUtilities.TimeVaryingInputs.extrapolation_bc
Base.in
Base.close(::ClimaUtilities.TimeVaryingInputs.AbstractTimeVaryingInput)

```
