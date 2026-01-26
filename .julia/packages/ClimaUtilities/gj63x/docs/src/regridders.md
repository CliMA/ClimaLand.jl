# [`Regridders`](@id regridder_module)

Simulations often need to import external data directly onto the computational
grid. The `Regridders` module implements different schemes to accomplish this
goal.

Currently, `Regridders` comes with two implementations:
1. `TempestRegridder` uses
   [TempestRemap](https://github.com/ClimateGlobalChange/tempestremap) (through
   `ClimaCoreTempestRemap`) to perform conservative interpolation onto lat-long
   grids. `TempestRegridder` only works for single-threaded CPU runs and works
   directly with files.
2. `InterpolationsRegridder` uses
   [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) to
   perform non-conservative linear interpolation onto lat-long(-z) and x-y-z grids.
   `InterpolationsRegridder` works directly with data.

!!! note

    While the `Regridders` can be used independently, most users will find their
    needs are immediately met by the [`SpaceVaryingInputs` and
    `TimeVaryingInputs`](@ref) interfaces. These higher-level objects implement
    everything that is needed to read a file to the model grid (internally using the
    `Regridders`).

If a regridder type is not specified by the user, and multiple are available,
the `InterpolationsRegridder` will be used by default. At least one regridder
extension must be loaded to be able to use regridding.

## [`TempestRegridder`](@id tempest_regridder)

> This extension is loaded when loading `ClimaCoreTempestRemap`

`TempestRegridder` performs conservative interpolation of lat-lon grids onto
computational domains. `TempestRegridder` performs all the interpolation ahead
of time and saves the regridded fields to HDF5 files that can be read during the
simulation.

Currently, `TempestRegridder` does not support regridding on 3D domains.
The `InterpolationsRegridder` described below can be used for these cases.

### Example

Assuming `target_space` is a `ClimaCore` 2D spherical field, the input data is
the variable `u` in the file `era5_example.nc`, and we want to read the data
associated with date `target_date`.

```julia
import ClimaUtilities.Regridders
import ClimaCoreTempestRemap
# Loading ClimaCoreTempest automatically loads TempestRegridder

reg = Regridders.TempestRegridder(target_space, "regrid_output", "u", "era5_example.nc")
# When reg is created, the variable `u` is regridded and the output files
# are saved to the `regrid_output` folder

regridded_u = Regridders.regrid(reg, target_date)
```

## [`InterpolationsRegridder`](@id interp_regridder)

> This extension is loaded when loading `ClimaCore` and `Interpolations`

`InterpolationsRegridder` performs linear interpolation of input data (linear
along each direction) and returns a `ClimaCore` `Field` defined on the
`target_space`.

!!! note "Did you know?"
    With versions of ClimaUtilities after v0.1.24, you can choose the mode of
    the gridded interpolation with the keyword argument `interpolation_method`.
    For example, to do constant interpolation, you can pass
    `interpolation_method = Interpolations.Constant()` as a keyword argument
    when constructing the regridder.

Currently, `InterpolationsRegridder` only supports spherical shells and extruded
spherical shells (but it could be easily extended to other domains).

!!! note

    It is easy to change the spatial interpolation type if needed.

`InterpolationsRegridder` are created once, they are tied to a `target_space`,
but can be used with any input data. With MPI runs, every process computes the
interpolating function. This is always done on the CPU and moved to GPU for
accelerated runs.

By default, `InterpolationsRegridder` assumes you are interpolating on a globe
and the default extrapolation boundary conditions are: periodic (along
longitudes), copy (along latitude), and throwing an error (along z). These can
be changed by passing the `extrapolation_bc` to the constructor of the regridder.

## FAQ: How do I enable linear extrapolation in z?

Create the regridder like this:
```julia
import Interpolations as Intp

extrapolation_bc = (Intp.Periodic(), Intp.Flat(), Intp.Linear())
regridder = InterpolationsRegridder(target_space; extrapolation_bc)
```

### Example

Assuming `target_space` is a `ClimaCore` 2D spherical field.
```julia
import ClimaUtilities.Regridders
import ClimaCore, Interpolations
# Loading ClimaCore and Interpolations automatically loads InterpolationsRegridder

reg = Regridders.InterpolationsRegridder(target_space)

# Now we can regrid any data
lon = collect(-180:1:180)
lat = collect(-90:1:90)
# It has to be lon, lat (because this is the assumed order in the CF conventions)
dimensions = (lon, lat)

data = rand((length(lon), length(lat)))

interpolated_data = Regridders.InterpolationsRegridder(reg, data, dimensions)
interpolated_2data = Regridders.InterpolationsRegridder(reg, 2 .* data, dimensions)
```

## API

```@docs
ClimaUtilities.Regridders.TempestRegridder
ClimaUtilities.Regridders.InterpolationsRegridder
ClimaUtilities.Regridders.regrid
```
