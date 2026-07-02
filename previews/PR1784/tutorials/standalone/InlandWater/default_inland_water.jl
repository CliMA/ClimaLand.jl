# # Default Inland Water Tutorial
# This shows how to setup a standalone inland water model in a
# single column, and run a simulation. With such a simple model, this
# is mostly useful for assessing the physicality of our equations.
import ClimaParams
using ClimaLand
import ClimaLand.Parameters as LP
using ClimaLand.InlandWater
using Dates
using ClimaDiagnostics
using ClimaLand.Domains: Point
using ClimaLand: prescribed_forcing_era5
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, Statistics
import ClimaLand.LandSimVis as LandSimVis;
# Set the floating point type and get the parameter dictionary
FT = Float32
toml_dict = LP.create_toml_dict(FT);
# Pick a location where there is inland waters - like the Caspian sea,
# and make the column domain:
longlat = FT.((50.6689, 41.9350))
z_sfc = FT(0)
domain = Point(; z_sfc, longlat)
surface_space = domain.space.surface;
# Obtain the forcing data for one year at this location:
start_date = DateTime(2008);
stop_date = start_date + Year(3)
forcing = prescribed_forcing_era5(
    start_date,
    stop_date,
    domain.space.surface,
    toml_dict,
    FT;
    use_lowres_forcing = true,
);
# Create the default slab lake/inland sea model:
lake = InlandWater.SlabLakeModel(FT, domain, forcing, toml_dict);
# We next make and solve the simulation object. Most of the defaults we won't
# explain here, but note that we pass in a (default) function which sets
# the inital condition for the lake prognostic internal energy. This
# sets the lake temperature equal to the air temperature.
set_ic! = ClimaLand.Simulations.set_lake_initial_conditions!
Δt = 450.0 # seconds
diag_writer = ClimaDiagnostics.Writers.DictWriter(); # conveniently stores the diagnostics in memory
diagnostics = ClimaLand.Diagnostics.default_diagnostics(
    lake,
    start_date;
    output_vars = ["tlake", "alake", "shflake", "lhflake", "rnlake", "qlake"],
    output_writer = diag_writer,
    reduction_period = :hourly,
);
simulation = ClimaLand.Simulations.LandSimulation(
    start_date,
    stop_date,
    Δt,
    lake;
    set_ic!,
    diagnostics,
    user_callbacks = (),
)
ClimaLand.Simulations.solve!(simulation);

# Let's visualize the results:
LandSimVis.make_timeseries(
    simulation;
    short_names = ["tlake", "alake", "shflake", "lhflake", "rnlake", "qlake"],
    plot_stem_name = "timeseries",
    spinup_date = stop_date - Year(1),
);
# ![](tlake_timeseries.png)
# ![](alake_timeseries.png)
# ![](shflake_timeseries.png)
# ![](lhflake_timeseries.png)
# ![](rnlake_timeseries.png)
# ![](qlake_timeseries.png)

# Atmospheric forcing data:
# [Hersbach2020](@citet)
