# # Changing LandModel Parameterizations
# In the [Canopy, soil, and snow](docs/src/tutorials/integrated/snowy_land_fluxnet_tutorial.jl)
# tutorial, we ran the integrated `LandModel`
# at a Fluxnet site using all of the default parameterizations and parameters.

# In two other tutorials ([Changing soil parameterizations](docs/src/tutorials/standalone/Soil/changing_soil_parameterizations.jl)
# and [Changing canopy parameterizations](docs/src/tutorials/standalone/Canopy/changing_canopy_parameterizations.jl))
# we explored how to change the parameterizations of the standalone soil and canopy models, respectively.

# This tutorial will combine these two streams of work to set up a `LandModel`
# with non-default parameterizations within the soil and canopy components.

# This time, we'll use non-default parameterizations for one canopy component:
# - radiative transfer: change from the default `TwoStreamModel` to `BeerLambertModel`
# and one snow component:
# - snow albedo: change from the default `ConstantAlbedo` to `ZenithAngleAlbedoModel`

# # Fluxnet simulations with the full land model: snow, soil, canopy

# As in the previous LandModel tutorial, we'll run at the Niwot Ridge site.
# The forcing data was obtained from
# AmeriFlux FLUXNET: https://doi.org/10.17190/AMF/1871141

# Citation: Peter D. Blanken, Russel K. Monson, Sean P. Burns,
# David R. Bowling, Andrew A. Turnipseed (2022),
# AmeriFlux FLUXNET-1F US-NR1 Niwot Ridge Forest (LTER NWT1),
# Ver. 3-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1871141

# # Preliminary Setup
using Dates
import ClimaParams as CP
using ClimaDiagnostics
using ClimaLand
using ClimaLand.Domains: Column, obtain_surface_domain
using ClimaLand.Simulations
using ClimaLand.Snow
using ClimaLand.Canopy
import ClimaLand.Parameters as LP
using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis;

# Define the floating point precision desired (64 or 32 bit), and get the
# parameter set holding constants used across CliMA Models.

const FT = Float32;
default_params_filepath =
    joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml");
toml_dict = LP.create_toml_dict(FT, default_params_filepath);
earth_param_set = LP.LandParameters(toml_dict);

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
    FluxnetSimulations.get_data_dates(site_ID, time_offset);
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
    earth_param_set,
    FT,
);
# LAI for the site - this uses our interface for working with MODIS data.
LAI = ClimaLand.Canopy.prescribed_lai_modis(
    domain.space.surface,
    start_date,
    stop_date,
);

# # Setup the integrated model

# First, we need to set up the component models we won't be using the defaults for.
# Since we are constructing these outside of the `LandModel` constructor,
# we need to provide some additional inputs that were omitted in the previous
# `LandModel` tutorial:
# - `surface_domain`: the surface of this simulation domain, which the snow and canopy models will use
# - `prognostic_land_components`: the prognostic land components, which must be consistent across all components
# - `ground`: the canopy ground conditions, which are prognostic since we're running with the soil
surface_domain = obtain_surface_domain(domain);
prognostic_land_components = (:canopy, :snow, :soil, :soilco2);
ground = ClimaLand.PrognosticGroundConditions{FT}();

# First, we will set up the snow model.
# We will use the `ZenithAngleAlbedoModel` for the snow albedo parameterization.
# This parameterization uses the zenith angle of the sun to determine the albedo,
# as opposed to the default `ConstantAlbedoModel` which uses a temporally and
# spatially constant snow albedo.
α_0 = FT(0.6) # parameter controlling the minimum snow albedo
Δα = FT(0.06) # parameter controlling the snow albedo when θs = 90∘
k = FT(2) # rate at which albedo drops to its minimum value with zenith angle
α_snow = Snow.ZenithAngleAlbedoModel(α_0, Δα, k);

# Now we can create the `SnowModel` model with the specified snow albedo parameterization.
snow = Snow.SnowModel(
    FT,
    surface_domain,
    forcing,
    toml_dict,
    Δt;
    prognostic_land_components,
    α_snow,
);

# Next, let's set up the canopy model using the `BeerLambertModel` radiative transfer parameterization
# with custom parameters. We already did this in the [Changing canopy parameterizations](docs/src/tutorials/standalone/Canopy/changing_canopy_parameterizations.jl)
# tutorial, so it should be familiar. In that example we showed three different ways to construct the `BeerLambertModel`.
# Here, we will use the third method, which takes the parameters object directly.
G_Function = Canopy.ConstantGFunction(FT(0.5)); # leaf angle distribution value 0.5
α_PAR_leaf = 0.2; # albedo in the PAR band
α_NIR_leaf = 0.3; # albedo in the NIR band
Ω = 1; # clumping index
radiative_transfer_parameters = Canopy.BeerLambertParameters(
    toml_dict;
    G_Function,
    α_PAR_leaf,
    α_NIR_leaf,
    Ω,
);
radiative_transfer = Canopy.BeerLambertModel(radiative_transfer_parameters);

# Now we can create the `CanopyModel` model with the specified radiative transfer
# parameterizations passed as [keyword arguments](https://docs.julialang.org/en/v1/manual/functions/#Keyword-Arguments).
(; atmos, radiation) = forcing;
canopy = Canopy.CanopyModel{FT}(
    surface_domain,
    (; atmos, radiation, ground),
    LAI,
    toml_dict;
    prognostic_land_components,
    radiative_transfer,
);

# Now we can construct the integrated `LandModel`. Since we want to use the
# defaults for the soil model, we don't need to provide anything for it.
# For the canopy and snow models, we'll provide the models we just set up.

land_model = LandModel{FT}(forcing, LAI, toml_dict, domain, Δt; snow, canopy);
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
    average_period = :hourly,
);

# Choose how often we want to update the forcing.
data_dt = Second(FluxnetSimulations.get_data_dt(site_ID));
updateat = Array(start_date:data_dt:stop_date);

# Now we can construct the simulation object and solve it.
simulation = Simulations.LandSimulation(
    start_date,
    stop_date,
    Δt, # seconds
    land_model;
    set_ic!,
    updateat,
    user_callbacks = (),
    diagnostics,
);
solve!(simulation);

# # Plotting results
LandSimVis.make_diurnal_timeseries(
    simulation;
    short_names = ["shf", "lhf", "swu", "lwu"],
    spinup_date = start_date + Day(20),
    plot_stem_name = "US_NR1_diurnal_timeseries_parameterizations",
);
# ![](lwu_US_NR1_diurnal_timeseries_parameterizations.png)
# ![](swu_US_NR1_diurnal_timeseries_parameterizations.png)
# ![](shf_US_NR1_diurnal_timeseries_parameterizations.png)
# ![](lhf_US_NR1_diurnal_timeseries_parameterizations.png)
LandSimVis.make_timeseries(
    simulation;
    short_names = ["swc", "si", "swe"],
    spinup_date = start_date + Day(20),
    plot_stem_name = "US_NR1_timeseries_parameterizations",
);
# ![](swc_US_NR1_timeseries_parameterizations.png)
# ![](si_US_NR1_timeseries_parameterizations.png)
# ![](swe_US_NR1_timeseries_parameterizations.png)

# Now you can compare these plots to those generated in the default LandModel Fluxnet tutorial.
# How are the results different? How are they the same?
