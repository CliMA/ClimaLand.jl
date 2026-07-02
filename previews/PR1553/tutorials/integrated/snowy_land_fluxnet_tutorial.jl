# # [Fluxnet simulations with the full land model: snow, soil, canopy](@id soilcanopy_fluxnet)

# In the
# [SoilCanopyModel tutorial](@ref "Fluxnet simulations with an integrated soil and canopy model"),
# we demonstrated how to run the an integrated model with a soil and
# canopy component at the US-MOz fluxnet site.
# Here we add in a snow component, and run at the Niwot Ridge site instead.
# The forcing data was obtained from
# [AmeriFlux FLUXNET](https://doi.org/10.17190/AMF/1871141)
# [Blanken2022](@citet)

# The focus of this tutorial is to learn the steps towards setting up and
# running an integrated simulation, and less on the parameterization
# choices. As such, the default parameters are implicitly set.
# To experiment with modularity in the parameters and parameterizations, please see the [soil parameterizations tutorial](@ref "Changing Soil Parameterizations"),
# the [canopy parameterizations tutorial](@ref "Changing Canopy Parameterizations"),
# or the [snowy land parameterizations tutorial](@ref "Changing LandModel Parameterizations").

# # Preliminary Setup
using Dates
import ClimaParams as CP
using ClimaDiagnostics
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Simulations
import ClimaLand.Parameters as LP
using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis;

# Define the floating point precision desired (64 or 32 bit), and get the
# parameter set holding constants used across CliMA Models.

const FT = Float32;
toml_dict = LP.create_toml_dict(FT);

# We will use prescribed atmospheric and radiative forcing from the
# US-NR1 tower.  We also
# read in the MODIS LAI and let that vary in time in a prescribed manner.
site_ID = "US-NR1";
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID);
# Get the latitude and longitude in degrees, as well as the
# time offset in hours of local time from UTC
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val));
# Get the height of the sensors in m
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val));
# Set a start and stop date of the simulation in UTC, as well as
# a timestep in seconds
(start_date, stop_date) =
    FluxnetSimulations.get_data_dates(site_ID, time_offset)
Δt = 450.0;

# Setup the domain for the model. This corresponds to
# a column of 2m in depth, with 10 equally spaced layers.
# The lat and long are provided so that we can look up default parameters
# for this location using the default ClimaLand parameter maps.
zmin = FT(-2) # in m
zmax = FT(0) # in m
domain = Column(; zlim = (zmin, zmax), nelements = 10, longlat = (long, lat));

# Forcing data for the site - this uses our interface for working with Fluxnet data
forcing = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    toml_dict,
    FT,
);
# LAI for the site - this uses our interface for working with MODIS data.
LAI = ClimaLand.Canopy.prescribed_lai_modis(
    domain.space.surface,
    start_date,
    stop_date,
);

# # Setup the integrated model

# We want to simulate the canopy-soil-snow system together, so the model type
# [`LandModel`](@ref "Integrated Land Model Types and methods")
# is chosen. Here we use the highest level model constructor, which uses default parameters,
# and parameterizations, for the soil, snow, and canopy models.

land_model = LandModel{FT}(forcing, LAI, toml_dict, domain, Δt);
set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
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
    reduction_period = :hourly,
);

# Choose how often we want to update the forcing.
# Choosing a frequency > the data frequency results in linear
# interpolation in time to the intermediate times.
data_dt = Second(FluxnetSimulations.get_data_dt(site_ID));

# Now we can construct the simulation object and solve it.
simulation = Simulations.LandSimulation(
    start_date,
    stop_date,
    Δt, # seconds
    land_model;
    set_ic!,
    updateat = Second(data_dt),
    user_callbacks = (),
    diagnostics,
);
solve!(simulation);

# We can optionally save the simulation parameters to a file for later reference.
# Here we specify the filepath where we want to save the parameters, and then
# ClimaParams handles the saving.
# Note that any parameters overwritten via keyword arguments when constructing
# models will not be reflected in this file (in this example there are none).
parameter_log_file = "snowy_land_fluxnet_parameters.toml"
CP.log_parameter_information(toml_dict, parameter_log_file)

# # Plotting results
LandSimVis.make_diurnal_timeseries(
    simulation;
    short_names = ["shf", "lhf", "swu", "lwu"],
    spinup_date = start_date + Day(20),
    plot_stem_name = "US_NR1_diurnal_timeseries",
);
# ![](lwu_US_NR1_diurnal_timeseries.png)
# ![](swu_US_NR1_diurnal_timeseries.png)
# ![](shf_US_NR1_diurnal_timeseries.png)
# ![](lhf_US_NR1_diurnal_timeseries.png)
LandSimVis.make_timeseries(
    simulation;
    short_names = ["swc", "si", "swe"],
    spinup_date = start_date + Day(20),
    plot_stem_name = "US_NR1_timeseries",
);
# ![](swc_US_NR1_timeseries.png)
# ![](si_US_NR1_timeseries.png)
# ![](swe_US_NR1_timeseries.png)
