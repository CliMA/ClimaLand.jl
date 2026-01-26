# NEWS

v0.2.14
-------

## Bug fixes

Fixed the default target coordinates for the `NetCDFWriter` for boxes with only
one point in the horizontal space.

v0.2.13
-------

## Features

### Add time and date bounds for netCDF files
In a netCDF file produced by ClimaDiagnostics, there are now the dimensions
`time_bnds` and `date_bnds`. Each time or date value is a representative of the
corresponding time or date bound. For example, if the time is `10.0` and
corresponding time bound is `[0.0, 10.0]`, then the time of `10.0` represents
the interval `[0.0, 10.0]`. If one knows that the data represents a time
average, then the time of `10.0` is the time average over the interval
`[0.0, 10.0]`.

### NetCDF writer now defaults to a reasonable number of points

`ClimaDiagnostics.Writers.NetCDF` now has a new default argument that depends on
the input Space. With this new default, obtained by calling the
`ClimaDiagnostics.Writers.default_num_points(space)` function, the output
diagnostics will be sampled with approximately the same resolution as the given
`space`.

### Support for `lazy`

Starting version `0.2.13`, `ClimaDiagnostics` supports diagnostic variables
specified with un-evaluated expressions (as provided by
[LazyBroadcast.jl](https://github.com/CliMA/LazyBroadcast.jl)).

Instead of
```julia
function compute_ta!(out, state, cache, time)
    if isnothing(out)
        return state.ta .- 273.15
    else
        out .= state.ta .- 273.15
    end
end
```
You can now write
```julia
import LazyBroadcast: lazy

function compute_ta(state, cache, time)
    return lazy.(state.ta .- 273.15)
end
```
Or, for `Field`s
```julia
function compute_ta(state, cache, time)
    return state.ta
end
```

v0.2.12
-------

## Features

### Support with ITime

This release adds support for working with `ITime`s. In particular,
`EveryCalendarDtSchedule` is compatible with `ITime` and a new constructor is
provided to make an `EveryCalendarDtSchedule` using `ITime`s. Lastly, there are
small changes to the writers to support `ITime`s.

v0.2.12
-------
## Bug fixes

- `NetCDFWriter` now correctly writes purely vertical and point spaces.

v0.2.11
-------
## Bug fixes

- Times in `DictWriter` are now correctly sorted.

v0.2.10
-------
## Bug fixes

Fixed broken `start_date` feature.

v0.2.9
-------

## Features

### Add a `start_date` attribute to NetCDFWriter PR [#94](https://github.com/CliMA/ClimaDiagnostics.jl/pull/94).

Prior to this version, users had to go to the simulation to find the start date.
It can now be saved as an attribute, making it easily accessible.
To do so, users need to pass the kwarg `start_date` when calling `NetCDFWriter`.

## Bug fixes

### Acquiring ownership with `compute!` PR [#88](https://github.com/CliMA/ClimaDiagnostics.jl/pull/88).

Prior to this version, `ClimaDiagnostics` would directly store use the output
returned by `compute!` functions the first time they are called. This leads to
problems when the output is a reference to an existing object since multiple
diagnostics would modify the same object. Now, `ClimaDiagnostics` makes a copy
of the return object so that it is no longer necessary to do so in the
`compute!` function.

### Correctly de-duplicate `ScheduledDiagnostics` [#93](https://github.com/CliMA/ClimaDiagnostics.jl/pull/93).

This version fixes a bug where `ScheduledDiagnostics` were not correctly
de-duplicated because `==` was not implemented correctly.

v0.2.8
-------

## Bug fixes

- `IntegratorWithDiagnostics` advertised a feature that was not implemented:
  `IntegratorWithDiagnostics` claimed that passing `state_name` and `cache_name`
  would allow users to customize the name of the state and cache inside the
  integrator. Now, this is implemented.

v0.2.7
-------

## Bug fixes

- `scheduled_diagnostics` are now internally saved as vectors instead of tuples.
  This has significant compile-time/inference benefits.

v0.2.6
-------

## Features

### More matadata in NetCDF files

Release `0.2.6` improves compatibility with CF conventions by adding
- standard and long name for the `time`, `longitude`, and `latitude` dimensions

## Bug fixes

- The default constructor for `ScheduleDiagnostic`s no longer uses reference of
  `Schedule`s but create a new copy.

### Deprecations

`reference_date` was renamed to `start_date` and `t_start` was dropped from the
constructors for the schedules. These changes are due to the fact that these
arguments should not be needed.

v0.2.5
-------

## Features

### Add support for box spaces with LatLong points in `NetCDFWriter`.

The `NetCDFWriter` can now work with regional boxes with `LatLong` points. Due
to incompatibility in `ClimaCore`, only `LatLong` points are supported (and not
`LongLat` points). This means that the box has to be created with latitude on
the x axis and longitude on the y axis.

## Bug fixes

- Ensure that `DictWriter` can only be constructed with dictionary-like objects.

v0.2.4
-------

- Add `EveryCalendarDtSchedule` for schedules with calendar periods.


v0.2.3
-------

- Detect and ignore duplicated diagnostics.

v0.2.2
-------

- Fix support for caches that are not NamedTuples in `NetCDFWriter`.

v0.2.1
-------

- Fix support for purely horizontal spaces in `NetCDFWriter`.

v0.2.0
-------

- ![][badge-💥breaking] The `NetCDFWriter` now outputs points on the vertical levels by default.
- ![][badge-💥breaking] `disable_vertical_interpolation` is removed in favor of `z_sampling_method`.

[badge-💥breaking]: https://img.shields.io/badge/💥BREAKING-red.svg
