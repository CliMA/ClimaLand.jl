# ClimaUtilities.jl

`ClimaUtilities.jl` provides a toolkit of functions that cover needs that are
shared across repositories within the [`CliMA`](https://github.com/CliMA)
project.

`ClimaUtilities.jl` is designed to have the minimum possible number of direct
dependencies. Instead, everything is implemented in [Julia
extensions](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions))
and modules are conditionally loaded when key packages are imported.

`ClimaUtilities.jl` also aims to provide an abstraction to commonly required
features (e.g., regridding), so that improvements can be made under the hood
without affecting users.

## Features

### `ClimaArtifacts`

`ClimaArtifacts` provides a MPI-safe way to lazily download artifacts and to tag
artifacts that are being accessed in a given simulation.

### `SpaceVaryingInputs` and `TimeVaryingInputs`

`SpaceVaryingInputs` and `TimeVaryingInputs` provide functions to seamlessly map
input to a `Field`. The input could be a function, a 1D array, a NetCDF file,
and it could be static or time varying. `SpaceVaryingInputs` and
`TimeVaryingInputs` objects can be `evaluate!`d to set the value of `Field`
(potentially using the `Regridders`). `TimeVaryingInputs` implement linear and
nearest neighbor interpolations.

### `FileReaders`

The `FileReaders` module provides a way to efficiently read data from files.
Efficiently might mean chunked/threaded/cached/something else. Currently, this
is mostly and interface barrier to provide a path for future improvements.

### `Regridders`

`ClimaUtilities` comes with two modules to map rectangular grids two (extruded)
finite spectral elements, `InterpolationsRegridder` and [`TempestRegridder`](@ref tempest_regridder).
These modules are primarily used to ingest data and resample it onto the
computational grid.

#### `InterpolationsRegridder`

`InterpolationsRegridder` uses
[Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) to perform
non-conservative linear interpolation onto lat-long(-z) grids.
`InterpolationsRegridder` is fully compatible with MPI/GPUs.

#### `TempestRegridder`

`TempestRegridder` uses
[TempestRemap](https://github.com/ClimateGlobalChange/tempestremap) to perform
conservative interpolation onto lat-long grids.

!!! note

    At the moment, `TempestRegridder` *does not* support MPI/GPUs and can
    only perform interpolation onto lat-long grids (not on z).

### `OnlineLogging`

`OnlineLogging` provides tools to produce informative messages while a
simulation is running (e.g., current/average timing information).

### `OutputPathGenerator`

`OutputPathGenerator` handles the directory structure for the output of a
simulation. If you are a package developer, use this module to set up the output
path for your simulation.

### `DataHandling`

The `DataHandling` module bundles a `Regridder` and a `FileReader` together to
serve regridded fields at a given time upon request. The main interface for
`DataHandling` is `regridded_snapshot(data_handler, date)`, a function that
returns a `Field` with data read from file for the given `date`. The
`DataHandler` maintains an least-recently-used (LRU) cache of regridded fields
to amortize the cost of (expensive) regridding operations.

### `DataStructures`

The `DataStructures` module implements helpful data structures to be used by
other ClimaUtilities.jl modules or external packages. Currently it contains an
LRU cache that is used in `DataHandlingExt` and `NCFileReaderExt`.
