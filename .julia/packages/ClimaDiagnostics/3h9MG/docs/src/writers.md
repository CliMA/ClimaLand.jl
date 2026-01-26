# Saving the diagnostics

Writers are needed to save the computed diagnostics.

`ClimaDiagnostics` comes with three writers:
- `NetCDFWriter`, to interpolate and save to NetCDF files;
- `DictWriter`, to save `Field`s to dictionaries in memory;
- `HDF5Writer`, to save `Field`s to HDF5 files.

(There is an additional `DummyWriter` that does nothing. It is mostly used
internally for testing and debugging.)

Users are welcome to implement their own writers. A writer has to be a subtype
of `AbstractWriter`and has to implement the `interpolate_field!` and
`write_field!` methods. `interpolate_field!` can return `nothing` is no
interpolation is needed.

## `NetCDFWriter`

The `NetCDFWriter` resamples the input `Field` to a rectangular grid and saves
the output to a NetCDF file.

The `NetCDFWriter` relies on the `Remappers` module in `ClimaCore` to
interpolate onto the rectangular grid. Horizontally, this interpolation is a
Lagrange interpolation, vertically, it is a linear. This interpolation is not
conservative. Also note that, the order of vertical interpolation drops to zero
in the first and last vertical elements of each column.

To create a `NetCDFWriter`, you need to specify the source `ClimaCore` `Space`
and the output directory where the files should be saved. By default, the
`NetCDFWriter` appends to existing files and create new ones if they do not
exist. The `NetCDFWriter` does not overwrite existing data and will error out if
existing data is inconsistent with the new one.

Optionally (recommended), you can pass an optional argument `start_date`, which
will be saved as an attribute of your NetCDF file, easily accessible.

`NetCDFWriter`s take as one of the inputs the desired number of points along
each of the dimensions. For the horizontal dimensions, points are sampled
linearly. For the vertical dimension, the behavior can be customized by passing
the `z_sampling_method` variable. When `z_sampling_method =
ClimaDiagnostics.Writers.LevelMethod()`, points evaluated on the grid levels
(and the provided number of points ignored), when `z_sampling_method =
ClimaDiagnostics.Writers.FakePressureLevelsMethod()`, points are sampled
uniformly in simplified hydrostatic atmospheric model.

The output in the `NetCDFWriter` roughly follows the CF conventions.

Each `ScheduledDiagnostic` is output to a different file with name determined by
calling the `output_short_name` on the `ScheduledDiagnostic`. Typically, these
files have names like `ta_1d_max.nc`, `ha_20s_inst.nc`, et cetera. The files
define their dimensions (`lon`, `lat`, `z`, ...). Time is always the first
dimension is any dataset.

Do not forget to close your writers to avoid file corruption!

To help reducing data loss, `NetCDFWriter` can force __syncing__, i.e. flushing
the values to disk. Usually, NetCDF buffers writes to disk (because they are
expensive), meaning values are not immediately written but are saved to disk in
batch. This can result in data loss, and it is often useful to force NetCDF to
write to disk (this is especially the case when working with GPUs). To do so,
you can pass the `sync_schedule` function to the constructor of `NetCDFWriter`.
When not `nothing`, `sync_schedule` is a callable that takes one argument (the
`integrator`) and returns a bool. When the bool is true, the files that were
modified since the last `sync` will be `sync`ed. For example, to force sync
every 1000 steps, you can pass the
`ClimaDiagnostics.Schedules.DivisorSchedule(1000)` schedule. By default, on
GPUs, we call `sync` at the end of every time step for those files that need to
be synced.

Variables are saved as datasets with attributes, where the attributes include
`long_name`, `standard_name`, `units`...

!!! note
    The `NetCDFWriter` cannot save raw `ClimaCore.Fields`, only fields that are
    resampled onto a Cartesian grids are supported. If you need such capability,
    consider using the [`ClimaDiagnostics.Writers.HDF5Writer`](@ref).

```@docs
ClimaDiagnostics.Writers.NetCDFWriter(space, output_dir; num_points, sync_schedule, z_sampling_method)
ClimaDiagnostics.Writers.interpolate_field!(writer::ClimaDiagnostics.Writers.NetCDFWriter, field, diagnostic, u, p, t)
ClimaDiagnostics.Writers.write_field!(writer::ClimaDiagnostics.Writers.NetCDFWriter, field, diagnostic, u, p, t)
ClimaDiagnostics.Writers.sync(writer::ClimaDiagnostics.Writers.NetCDFWriter)
Base.close(writer::ClimaDiagnostics.Writers.NetCDFWriter)
```

Sampling methods for the vertical direction:
```@docs
ClimaDiagnostics.Writers.AbstractZSamplingMethod
ClimaDiagnostics.Writers.LevelsMethod
ClimaDiagnostics.Writers.FakePressureLevelsMethod
```


## `DictWriter`

The `DictWriter` is a in-memory writer that is particularly useful for
interactive work and debugging.

```@docs
ClimaDiagnostics.Writers.DictWriter()
ClimaDiagnostics.Writers.write_field!(writer::ClimaDiagnostics.Writers.DictWriter, field, diagnostic, u, p, t)
```

## `HDF5Writer`

 The `HDF5Writer` writes the `Field` directly to an HDF5 file in such a way that
it can be later read and imported using the `InputOutput` module in `ClimaCore`.

The `HDF5Writer` writes one file per variable per timestep. The name of the file
is determined by the `output_short_name` field of the `ScheduledDiagnostic` that
is being output.

> Note: The `HDF5Writer` in `ClimaDiagnostics` is currently the least developed
> one. If you need this writer, we can expand it.

```@docs; canonical=false
ClimaDiagnostics.Writers.HDF5Writer
ClimaDiagnostics.Writers.write_field!(writer::ClimaDiagnostics.Writers.HDF5Writer, field, diagnostic, u, p, t)
Base.close(writer::ClimaDiagnostics.Writers.HDF5Writer)
```
