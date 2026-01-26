ClimaUtilities.jl Release Notes
===============================

main
------

#### Allow regridding 2D data onto 2D spaces with LatLongZ coordinates. PR[#176](https://github.com/CliMA/ClimaUtilities.jl/pull/176)

Some 2D spaces use `LatLongPoint` coordinates and others use `LatLongZPoint`,
depending on how they were constructed. Both cases are now supported, so we can
read in 2D data onto a 2D space with either coordinate type.

#### Interpolation method for InterpolationsRegridder

`InterpolationsRegridder` now accepts the new keyword argument
`interpolation_method` for choosing the mode of the gridded interpolation. The
default is `Interpolations.Linear()`.

```julia
import ClimaUtilities.Regridders
import ClimaCore, Interpolations

linear_reg = Regridders.InterpolationsRegridder(
    target_space,
    interpolation_method = Intp.Linear(),
)
constant_reg = Regridders.InterpolationsRegridder(
    target_space,
    interpolation_method = Intp.Constant(),
)
```

v0.1.23
------

#### Support for Box Interpolations Regridder. PR[#151](https://github.com/CliMA/ClimaUtilities.jl/pull/151)

`Regridders.IterpolationsRegridder` now supports regridding on
`ClimaCore.Geometry.XYZPoint` objects which allows for interpolation
onto boxes and single column simulations.

v0.1.22
------

#### New integer time type, `TimeManager.ITime`. PR [#149](https://github.com/CliMA/ClimaUtilities.jl/pull/149)

`ClimaUtilities` now comes with a new type that represent times and date,
`ITime`. `ITime` stands for "integer time" and is a time type that does not
incur in floating point errors. `ITime`s also encode dates and support
arithmetic operations and fractions (to be used in timestepping loops).

```julia-repl
julia> using ClimaUtilities.TimeManager, Dates;
julia> time1 = ITime(5; period = Hour(20), start_date = DateTime(2012, 12, 21))
100 hours (2012-12-25T04:00:00) [counter = 5, period = 20 hours, start_date = 2012-12-21T00:00:00]
julia> seconds(time1)
360000.0
julia> date(time1)
2012-12-25T04:00:00
```

The [documentation](https://clima.github.io/ClimaUtilities.jl/dev/timemanager/)
provides further information about this new type and has a section dedicated to
helping developers port their codes to `ITime`s.

v0.1.21
------

#### Simplified `TimeVaryingInput0D`, removed `context` argument. PR [#148](https://github.com/CliMA/ClimaUtilities.jl/pull/148)

The need for a dedicated CUDA extension was by leveraging `ClimaCore` functions.
As a result, the code for `TimeVaryingInput0D` could be significantly
simplified, while attaining greater performance at the same time.

As a result, the keyword argument `context` is no longer required in
constructing this type of `TimeVaryingInput`s. In the future, the argument will
be removed.

v0.1.20
------

#### New module `OnlineLogging`. PR [#137](https://github.com/CliMA/ClimaUtilities.jl/pull/137)

A new module, `OnlineLogging`, provides tools to report messages while a
simulation is running. Currently, the only feature implemented is
`WallTimeInfo`, to report on timing information.

With `WallTimeInfo`, one can add reports like the following to their simulations
```
┌ Info: Progress
│   simulation_time = "2 seconds"
│   n_steps_completed = 20
│   wall_time_per_step = "10 milliseconds, 100 microseconds"
│   wall_time_total = "1 second, 10 milliseconds"
│   wall_time_remaining = "808 milliseconds, 35 microseconds"
│   wall_time_spent = "202 milliseconds, 8 microseconds"
│   percent_complete = "20.0%"
│   estimated_sypd = "0.027"
│   date_now = 2024-12-03T16:01:13.660
└   estimated_finish_date = 2024-12-03T16:01:14.660
```

Find more information in the
[documentation](https://clima.github.io/ClimaUtilities.jl/dev/onlinelogging.html).

#### Support reading time data across multiple files. PRs [#127](https://github.com/CliMA/ClimaUtilities.jl/pull/127), [#132](https://github.com/CliMA/ClimaUtilities.jl/pull/132)

`NCFileReader`s can now read multiple files at the same time. The files have to
contain temporal data for the given variable and they are aggregated along the
time dimension. To use this feature, just pass a vector of file paths to the
constructor.

This capability is also available to `DataHandler`s and `TimeVaryingInput`s. To
use this feature, just pass the list of files that contain your variable of
interested, for example
```julia
timevaryinginput = TimeVaryingInputs.TimeVaryingInput(["era5_1980.nc", "era5_1981.nc"],
                                                      "u",
                                                      target_space,
                                                      start_date = Dates.DateTime(1980, 1, 1),
                                                      regridder_type = :InterpolationsRegridder)
```
You can also compose variables
```julia
timevaryinginput = TimeVaryingInputs.TimeVaryingInput(["era5_1980.nc", "era5_1981.nc", "era5_1982.nc"],
                                                      ["u", "v"],
                                                      target_space,
                                                      start_date = Dates.DateTime(1980, 1, 1),
                                                      regridder_type = :InterpolationsRegridder,
                                                      compose_function = (x, y) -> sqrt(x^2 + y^2))
```

When you compose variables, pay attention that `TimeVaryingInput` implements
some heuristics to disambiguate the case where the passed list of files is split
along the time or the variable dimension. You can always pass a list of lists to
be explicit in your intentions. Read the
[documentation](https://clima.github.io/ClimaUtilities.jl/dev/datahandling.html#Heuristics-to-do-what-you-mean)
to learn more about this.

This capability is only available for the `InterpolationsRegridder`.

#### Reduced default size of cache in `DataHandler`. PR [#133](https://github.com/CliMA/ClimaUtilities.jl/pull/133)

The default cache size for regridded fields in `DataHandler` was reduced from
128 to 2, reducing the memory footprint. You can pass the `cache_max_size`
keyword argument to control this value.

#### Improved error handling

When using a function that depends on loading another package, the error now
tells the user which package should be loaded if the package has not been loaded
already. Furthermore, if there is an error with `@clima_artifact`, the error
message tells the user that the artifact might be undownloadable and recommends
downloading it themselves.

v0.1.19
------

#### `ClimaComms` is now a required dependency. PR [#128](https://github.com/CliMA/ClimaUtilities.jl/pull/128)

`ClimaComms` used to be an optional dependency and was turned into a required
one. The reason for this change is to improve robustness in MPI settings.

The new version should also further reduce the "Path not properly synced" error.

v0.1.18
------

### Bug fixes

- Fixed `@clima_artifact` for Julia 1.12. PR [#123](https://github.com/CliMA/ClimaUtilities.jl/pull/123)
- Increased sleep time across attempts in checking for synced filesystems. This
  should give the distributed filesystems more time to sync and reduce the
  occurrence of the "Path not properly synced" error. PR
  [#125](https://github.com/CliMA/ClimaUtilities.jl/pull/125)

v0.1.17
------

#### Changed signature for `OutputPathGenerator.detect_restart_file`. PR [#122](https://github.com/CliMA/ClimaUtilities.jl/pull/122)

The signature for `OutputPathGenerator.detect_restart_file` was changed. This
change is motivated by the fact that `OutputPathGenerator.detect_restart_file`
only makes sense with the `ActiveLinkStyle`. This led to a cumbersome usage in
downstream packages where `ActiveLinkStyle` had to be imported just to call
`detect_restart_file`.

The old signature was deprecated.

### Bug fixes

#### Fixed parsing of `date`s in NetCDF files. PR [#122](https://github.com/CliMA/ClimaUtilities.jl/pull/122)

This release fixes an issue with reading NetCDF files that have `date` as a time
dimension. Only the YYYYMMDD format is currently supported.

v0.1.16
------

### Bug fixes

- `Utils.isequispaced` is now more efficient: it fails fast and does not allocate
as much. More redundant allocations due to `Utils.isequispaced` were fixed.

- Reduced allocations in regridding. New method `read!`.
Existing `read` now returns a copy of the data when returning from the cache.
PR [#119](https://github.com/CliMA/ClimaUtilities.jl/pull/119)

### Features

#### Added `OutputPathGenerator.detect_restart_file`. PR [#120](https://github.com/CliMA/ClimaUtilities.jl/pull/120)

A new function, `OutputPathGenerator.detect_restart_file`, has been added to
help automatically locate the most recent restart file generated by a
simulation. This function is particularly useful when restarting simulations
from a previous checkpoint. It works with the `ActiveLinkStyle` and can be
customized to handle different restart file naming conventions and sorting
criteria.

v0.1.15
-------

### Features

#### Added support for composing input variables. PR [#105](https://github.com/CliMA/ClimaUtilities.jl/pull/105/)

This allows a list of `varnames` (and possibly `file_paths`) to be passed to
`TimeVaryingInput` or `SpaceVaryingInput`, along with a `compose_function` to
compose them, as so:

```julia
# Define the pointwise composing function we want, a simple sum in this case
compose_function = (x, y) -> x + y
# Define pre-processing function to convert units of input
unit_conversion_func = (data) -> 1000 * data

data_handler = TimeVaryingInputs.TimeVaryingInput("era5_example.nc",
                                        ["u", "v"],
                                        target_space,
                                        start_date = Dates.DateTime(2000, 1, 1),
                                        regridder_type = :InterpolationsRegridder,
                                        file_reader_kwargs = (; preprocess_func = unit_conversion_func),
                                        compose_function)
```

See the `TimeVaryingInput` or `DataHandler` docs "NetCDF file input" sections
for more details.

### Deprecations

`reference_date` was renamed to `start_date` and `t_start` was dropped from the
constructors for `DataHandler` and `TimeVaryingInput`.

These changes are due to the fact that these arguments should not be needed.

v0.1.14
-------

- Reduce floating point errors in times for DataHandling. PR
  [#101](https://github.com/CliMA/ClimaUtilities.jl/pull/101)

v0.1.13
-------

- Fix `show` for `LRUCache` (and `iterate` with a state). PR
  [#96](https://github.com/CliMA/ClimaUtilities.jl/pull/96)

v0.1.12
-------

- Add support for interpolating while "period filling". PR
  [#85](https://github.com/CliMA/ClimaUtilities.jl/pull/85)
- Add support for boundary conditions in interpolation. PR
  [#84](https://github.com/CliMA/ClimaUtilities.jl/pull/84)
- Increased allocations in regridding. `read!` method removed. PR
  [#84](https://github.com/CliMA/ClimaUtilities.jl/pull/84)

v0.1.11
------

- Reduced allocations in regridding. New method `read!`. PR
  [#83](https://github.com/CliMA/ClimaUtilities.jl/pull/83)

v0.1.10
------

- Reduced allocations in regridding. New method `regridded_snapshot!`. PR
  [#72](https://github.com/CliMA/ClimaUtilities.jl/pull/72)

v0.1.9
------

- Extensions are internally reorganized, removing precompilation errors. PR
  [#69](https://github.com/CliMA/ClimaUtilities.jl/pull/69)

v0.1.8
------

- `generate_output_path(..., ::ActiveLinkStyle)` now returns the folder instead
  of the link. Links are still being created and managed. PR
  [#63](https://github.com/CliMA/ClimaUtilities.jl/pull/63)

v0.1.7
------

- Fix compatibility with ClimaComms 0.6. PR [#54](https://github.com/CliMA/ClimaUtilities.jl/pull/54)

v0.1.6
-------
- `OutputPathGenerator` now tries to create an active link when one is not available but some data is already there [#50](https://github.com/CliMA/ClimaUtilities.jl/pull/50)
- Fix compatibility with ClimaCore 0.14. PR [#50](https://github.com/CliMA/ClimaUtilities.jl/pull/50)

v0.1.5
-------
- Support passing down regridder and file reader arguments from higher level constructors. PR [#40](https://github.com/CliMA/ClimaUtilities.jl/pull/40)

v0.1.4
-------
- Fix and test MPI compatibility. PRs [#33](https://github.com/CliMA/ClimaUtilities.jl/pull/33), [#37](https://github.com/CliMA/ClimaUtilities.jl/pull/37)
- Select default regridder type if multiple are available. PR [#32](https://github.com/CliMA/ClimaUtilities.jl/pull/32)

v0.1.3
-------
- Add `DataStructures` module containing `LRUCache` object. PR [#35](https://github.com/CliMA/ClimaUtilities.jl/pull/35)
- Add `OutputPathGenerator`. PR [#28](https://github.com/CliMA/ClimaLand.jl/pull/28)

[badge-💥breaking]: https://img.shields.io/badge/💥BREAKING-red.svg
