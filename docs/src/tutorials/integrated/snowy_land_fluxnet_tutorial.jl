# # Fluxnet simulations with the full land model: snow, soil, canopy

# In the
# [`previous tutorial`](@ref https://clima.github.io/ClimaLand.jl/dev/generated/soil_canopy_tutorial/),
# we demonstrated how to run the an integrated model with a soil and
# canopy component at the US-MOz fluxnet site.
# Here we add in a snow component, and run the site at the Niwot Ridge site instead.
# Again, the focus of this tutorial is to learn the steps towards setting up and
# running an integrated simulation, and less on the parameterization choices.
# As such, the default parameters are implicitly set. To experiment with
# modularity in the parameters and parameterizations, please see tutorial X.

# # Preliminary Setup
using Dates
import ClimaParams as CP
using ClimaDiagnostics
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Simulations
import ClimaLand.Parameters as LP
using DelimitedFiles
FluxnetSimulationsExt =
    Base.get_extension(ClimaLand, :FluxnetSimulationsExt).FluxnetSimulationsExt;
using CairoMakie, ClimaAnalysis, GeoMakie, Poppler_jll, Printf, StatsBase
LandSimulationVisualizationExt =
    Base.get_extension(
        ClimaLand,
        :LandSimulationVisualizationExt,
    ).LandSimulationVisualizationExt;

# Define the floating point precision desired (64 or 32 bit), and get the
# parameter set holding constants used across CliMA Models.

const FT = Float32;
earth_param_set = LP.LandParameters(FT);

# We will use prescribed atmospheric and radiative drivers from the
# US-NR1 tower, which we read in here.  We also
# read in the MODIS LAI and let that vary in time in a prescribed manner.
site_ID = "US-NR1"
time_offset = 7 # Timezone (offset from UTC in hrs)
lat = FT(40.0329) # degree
long = FT(-105.5464) # degree
atmos_h = FT(21.5); # Height of the sensor at the site (m)
(start_date, stop_date) =
    FluxnetSimulationsExt.get_data_dates(site_ID, time_offset) # in UTC
Δt = 450.0; # seconds

# Forcing data for the site - this uses our interface for working with Fluxnet data
forcing = FluxnetSimulationsExt.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    earth_param_set,
    FT,
);
# LAI for the site - this uses our interface for working with MODIS data.
modis_lai_ncdata_path = ClimaLand.Artifacts.modis_lai_multiyear_paths(;
    start_date,
    end_date = stop_date,
)
LAI = ClimaLand.prescribed_lai_modis(
    modis_lai_ncdata_path,
    domain.space.surface,
    start_date,
);
# Setup the domain for the model:
zmin = FT(-5) # in m
zmax = FT(0) # in m
domain = Column(; zlim = (zmin, zmax), nelements = 10, longlat = (long, lat))

# # Setup the integrated model

# We want to simulate the canopy-soil-snow system together, so the model type
# [`LandModel`](https://clima.github.io/ClimaLand.jl/dev/APIs/ClimaLand/#LSM-Model-Types-and-methods)
# is chosen. Here we use the highest level model constructor, which uses default parameters,
# and parameterizations, for the soil, snow, and canopy models.
# A different tutorial will show you how to change these parameters and parameterizations.
land_model = LandModel{FT}(forcing, LAI, earth_param_set, domain, Δt);
set_ic! = FluxnetSimulationsExt.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    time_offset,
    land_model,
);
output_vars = ["swu", "lwu", "shf", "lhf", "swe", "swc", "si"]
diagnostics = ClimaLand.default_diagnostics(
    land_model,
    start_date;
    output_writer = ClimaDiagnostics.Writers.DictWriter(),
    output_vars,
    average_period = :hourly,
);

simulation = Simulations.LandSimulation(
    FT,
    start_date,
    stop_date,
    Δt, # seconds
    land_model;
    set_ic!,
    user_callbacks = (),
    diagnostics,
);
@time solve!(simulation);

# # Plotting results
LandSimulationVisualizationExt.make_diurnal_timeseries(
    simulation;
    short_names = ["shf", "lhf", "swu", "lwu"],
    spinup_date = start_date + Day(20),
);
# ![](diurnal_timeseries.png)
LandSimulationVisualizationExt.make_timeseries(
    simulation;
    short_names = ["swc", "si", "swe"],
    spinup_date = start_date + Day(20),
);
