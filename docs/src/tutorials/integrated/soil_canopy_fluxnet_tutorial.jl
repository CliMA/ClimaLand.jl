# # Fluxnet simulations with an integrated soil and canopy model

# In the
# [`previous tutorial`](@ref https://clima.github.io/ClimaLand.jl/dev/generated/soil_plant_hydrology_tutorial/),
# we demonstrated how to run the canopy model in
# standalone mode using prescribed values for the inputs of soil moisture
# and ground temperature
# into the canopy hydraulics model. However, ClimaLand can also
# integrate the canopy model with a soil model and timestep the two
# components together to simulate an interacting canopy-soil system. This tutorial
# demonstrates how to setup and run an integrated simulation, again using
# initial conditions, atmospheric and radiative flux conditions, and leaf
# area index observed at the US-MOz flux tower, a flux tower located within an
# oak-hickory forest in Ozark, Missouri, USA.
# The focus of this tutorial is to learn the steps towards setting up and
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
# US-MOz tower, which we read in here.  We also
# read in the MODIS LAI and let that vary in time in a prescribed manner.
site_ID = "US-MOz"
time_offset = 7 # Timezone (offset from UTC in hrs)
lat = FT(38.7441) # degree
long = FT(-92.2000) # degree
atmos_h = FT(32); # Height of the sensor at the site (m)

start_date = DateTime("2010-05-01", "yyyy-mm-dd") # in UTC
stop_date = DateTime("2010-09-01", "yyyy-mm-dd") # in UTC
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
zmin = FT(-2) # in m
zmax = FT(0) # in m
domain = Column(; zlim = (zmin, zmax), nelements = 10, longlat = (long, lat))

# # Setup the integrated model

# We want to simulate the canopy-soil system together, so the model type
# [`SoilCanopyModel`](https://clima.github.io/ClimaLand.jl/dev/APIs/ClimaLand/#LSM-Model-Types-and-methods)
# is chosen. Here we use the highest level model constructor, which uses default parameters,
# and parameterizations, for the soil and canopy models.
# A different tutorial will show you how to change these parameters and parameterizations.
land_model = SoilCanopyModel{FT}(forcing, LAI, earth_param_set, domain);
set_ic! = FluxnetSimulationsExt.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    time_offset,
    land_model,
);
output_vars = ["gpp", "swu", "lwu", "shf", "lhf"]
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
solve!(simulation);

# # Plotting results
LandSimulationVisualizationExt.make_diurnal_timeseries(
    simulation;
    short_names = ["gpp", "shf", "lhf", "swu", "lwu"],
    spinup_date = start_date + Day(20),
);
# ![](diurnal_timeseries.png)
