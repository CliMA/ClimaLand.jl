# # Fluxnet forcing data: LAI, radiation, and atmospheric variables

# In this tutorial, we will demonstrate how we read in forcing data at a Fluxnet site
# using our ClimaLand infrastructure. To see an example running a simulation at a
# Fluxnet site, please see the corresponding tutorials for the [SoilCanopyModel](docs/src/tutorials/integrated/soil_canopy_fluxnet_tutorial.jl) 
# or [LandModel](docs/src/tutorials/integrated/snowy_land_fluxnet_tutorial.jl)
# To access the forcing data (LAI from MODIS, SW\_d, LW\_d, T\_air, q\_air,
# P\_air, and precipitation from fluxtower data), you first need the
# the fluxtower site ID.
# Currently, ClimaLand provides an interface for working with four
# fluxtower sites; adding support for a much larger set of sites is
# in progress. The sites we support are Vaira Ranch (US_Var),
# Missouri Ozark (US-MOz), Niwot Ridge (US-NR1), and Harvard Forest
# (US-Ha1).

using Dates
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains: Column
import ClimaLand.Parameters as LP
using DelimitedFiles
using CairoMakie
import ClimaLand.FluxnetSimulations as FluxnetSimulations

# Define the floating point precision desired (64 or 32 bit), and get the
# parameter set holding constants used across CliMA Models.
const FT = Float32;
earth_param_set = LP.LandParameters(FT);

# Pick a site ID; convert the dash to an underscore:
site_ID = "US-NR1";
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)
# The functions we use below use multiple dispatch, treating the
# site_ID_val as a Julia type.

# First, we need the latitude and longitude of the site. These are used
# to get the zenith angle as a function of time, and to look up
# default parameters using the global ClimaLand parameter maps. We also
# need the offset of the local time of the site in hours from UTC. This is
# because ClimaLand simulations are carried out in UTC, while the fluxtower
# data is reported in local time.
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))

# ClimaLand also needs to know the height at which the atmospheric data
# was recorded. 
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

# It is also useful to know the bounds of the data,
# in UTC, to use as the start and stop date of the simulation.
(start_date, stop_date) =
    FluxnetSimulations.get_data_dates(site_ID, time_offset)

# Now we can construct the forcing objects. Under the hood, this
# function finds the local path to the fluxtower data (and downloads it
# if it is not present) and reads the data in. It then creates
# two objects, one called `atmos`, of type `PrescibedAtmosphere`, and
# one called `radiation`, of type `PrescribedRadiativeFluxes`. These
# encode the data in interpolating functions which allow us to
# estimate the forcing at any time during the simulation using linear
# interpolation across gaps in the data.
(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    earth_param_set,
    FT,
);
# The atmosphere object holds the air temperature, pressure, specific
# humidity, wind speed, and liquid and solid precipitation fluxes.
# Since many fluxtower sites do not measure the snow and rain fluxes
# separately, we estimate them internally. This is optional, and by providing
# the kwarg `split_precip = false`, you can change the behavior to
# return all measured precip as a liquid water flux (no snow).
# The radiation object holds the downwelling short and long wave
# radiative fluxes. The diffuse fraction is estimated internally
# using an empirical function, and the zenith angle is computed
# analytically.

# The simulation time is measured in seconds since the `start_date`. We
# can get the atmospheric temperature at the `start_date` as follows:
sim_time = 0.0
T = [NaN]
evaluate!(T, atmos.T, sim_time);
@show T
# Note that `evaluate!` updates `T` in place. This is important: in our
# simulations, we allocate memory for the forcing, and then update the values
# at each step. If we did not do this, the dynamic memory allocation
# would cause the simulation to be incredibly slow.
# We can plot the interpolated air temperature for the first day:
sim_times = 0.0:1800.0:86400.0 # one day in seconds
air_temps = [];
for sim_time in sim_times
    evaluate!(T, atmos.T, sim_time)
    push!(air_temps, T[1])
end
fig = CairoMakie.Figure()
ax = CairoMakie.Axis(
    fig[1, 1],
    ylabel = "Temperature (K)",
    xlabel = "Date",
    title = "Near-surface air temperature at Niwot Ridge",
)
lines!(ax, Second.(sim_times) .+ start_date, air_temps)
CairoMakie.save("air_temp.png", fig);
# ![](air_temp.png)

# We do something very similar with LAI, but now the data is coming
# from MODIS. In this case, we use nearest-neighbor interpolatation
# to move from the gridded
# MODIS data to the LAI at the latitude and longitude of the site.
# To do spatial interpolation, we need to create a ClimaLand domain.
domain =
    Column(; zlim = (FT(-3.0), FT(0.0)), nelements = 10, longlat = (long, lat))
surface_space = domain.space.surface;
# Get the paths to each year of MODIS data within the start and stop dates.
LAI = ClimaLand.prescribed_lai_modis(surface_space, start_date, stop_date);

# Just like with the air temperature, the LAI is an object that we can use
# to linearly interpolate observed LAI to any simulation time.

# It can also be useful to know the maximum LAI at a site. To do so, we
# can call, for the first year of data:
maxLAI =
    FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long)

# # Fluxnet comparison data

# To assess performance of the simulation, we need comparison data from the
# site. This is complicated since different sites provided different
# data (e.g. soil moisture and temperature and depth). ClimaLand
# provides a function to get the comparison data for a select number of
# variables:

# | Variable  | Unit       | ClimaLand short name | Fluxnet Column name
# |:---------| :--------:| :------------:| -----------:|
# |Gross primary prod. | mol/m^2/s |"gpp"|"GPP\_DT\_VUT\_REF"|
# |Latent heat flux | W/m^2/s |"lhf"|"LE\_CORR"|
# |Sensible heat flux | W/m^2/s |"shf"|"H\_CORR"|
# |Upwelling SW | W/m^2/s |"swu"|"SW\_OUT"|
# |Upwelling LW | W/m^2/s |"lwu"|"LW\_OUT"|
# |Soil water content |  |"swc"|"SWC\_F\_MDS\_1"|
# |Soil temperature |  K|"tsoil"|"TS\_F\_MDS\_1"|
# |Total precipitation | m/s|"precip"|"P\_F"|

comparison_data = FluxnetSimulations.get_comparison_data(site_ID, time_offset);
@show propertynames(comparison_data)
# If a column is missing, the column is not returned. Missing values
# are replaced with the mean of the not missing values.

# The data we use to force the simulations and to compare the results against
# were obtained from Ameriflux:

# US-Moz: https://doi.org/10.17190/AMF/1854370
# Citation: Jeffrey Wood, Lianhong Gu (2025), AmeriFlux FLUXNET-1F US-MOz Missouri Ozark
# Site, Ver. 5-7, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1854370

# US-NR1: https://doi.org/10.17190/AMF/1871141
# Citation: Peter D. Blanken, Russel K. Monson, Sean P. Burns,
# David R. Bowling, Andrew A. Turnipseed (2022),
# AmeriFlux FLUXNET-1F US-NR1 Niwot Ridge Forest (LTER NWT1),
# Ver. 3-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1871141

# US-Var: https://doi.org/10.17190/AMF/1993904
# Citation: Siyan Ma, Liukang Xu, Joseph Verfaillie, Dennis Baldocchi (2023),
# AmeriFlux FLUXNET-1F US-Var Vaira Ranch- Ione, Ver. 3-5, AmeriFlux AMP, (Dataset).
# https://doi.org/10.17190/AMF/1993904

# US-Ha1: https://doi.org/10.17190/AMF/1871137
# Citation: J. William Munger (2022), AmeriFlux FLUXNET-1F US-Ha1
# Harvard Forest EMS Tower (HFR1), Ver. 3-5, AmeriFlux AMP, (Dataset).
# https://doi.org/10.17190/AMF/1871137
