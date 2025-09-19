## Installation of Julia and ClimaLand

ClimaLand is provided as a Julia package, so it requires having Julia installed. Information about Julia packages is available on the [Julia website](https://julialang.org/packages/).

First, download and install Julia by following the instructions at [https://julialang.org/downloads/](https://julialang.org/downloads/).
Then, you can launch a [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) by running `julia --project` and install the
ClimaLand.jl package by running the following within the REPL:

```julia
using Pkg
Pkg.add(ClimaLand)
```

## Running a simple soil model

Now that we have Julia installed, let's set up and run a simple ClimaLand simulation!
You can follow along by copying the following code segments into your REPL.

First we import the necessary Julia packages:
```julia
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains
using ClimaLand.Soil
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
using Dates
```

Choose a floating point precision, and get the parameter set, which holds constants used across CliMA models.
Note: we use SI units unless otherwise specified.
See our [Physical Units](https://clima.github.io/ClimaLand.jl/stable/physical_units/) documentation for more information.
```julia
FT = Float32
toml_dict = LP.create_toml_dict(FT);
```

We will run this simulation on a column domain with 1 meter depth, at a lat/lon location
near Pasadena, California.

```julia
zmax = FT(0)
zmin = FT(-1.0)
longlat = FT.((-118.1, 34.1))
domain = Domains.Column(; zlim = (zmin, zmax), nelements = 10, longlat);
surface_space = domain.space.surface;
```

We choose the initial and final simulation times as `DateTime`s, and a timestep in seconds.
```julia
start_date = DateTime(2008);
stop_date = start_date + Second(60 * 60 * 72);
dt = 1000.0;
```

The soil model takes in 2 forcing objects, atmosphere and radiation,
which we read in from ERA5 data.
```julia
era5_ncdata_path =
    ClimaLand.Artifacts.era5_land_forcing_data2008_path(; lowres = true);
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    era5_ncdata_path,
    surface_space,
    start_date,
    toml_dict,
    FT,
);
```

Now, we can create the [`EnergyHydrology`](@ref ClimaLand.Soil.EnergyHydrology) model.
This constructor uses default parameters and parameterizations, but these can also be
overwritten, which we'll demonstrate in later tutorials.
```julia
model = Soil.EnergyHydrology{FT}(
    domain,
    (; atmos, radiation),
    toml_dict,
);
```

Define a function to set initial conditions for the prognostic variables.
```julia
function set_ic!(Y, p, t0, model)
    Y.soil.ϑ_l .= FT(0.24);
    Y.soil.θ_i .= FT(0.0);
    T = FT(290.15);
    ρc_s = Soil.volumetric_heat_capacity.(
        Y.soil.ϑ_l,
        Y.soil.θ_i,
        model.parameters.ρc_ds,
        model.parameters.earth_param_set,
    );
    Y.soil.ρe_int .=
        Soil.volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T,
            model.parameters.earth_param_set,
        );
end
```

Now construct the `LandSimulation` object, which contains the model
and additional timestepping information.
```julia
simulation = LandSimulation(start_date, stop_date, dt, model; set_ic!, user_callbacks = ());
```

Now we can run the simulation!
```julia
solve!(simulation);
```

That's it! For more complex examples and pretty plots,
please see the tutorial section of our docs.

Atmospheric forcing data citation:
Hersbach, Hans, et al. "The ERA5 global reanalysis."
Quarterly journal of the royal meteorological society 146.730 (2020): 1999-2049.
