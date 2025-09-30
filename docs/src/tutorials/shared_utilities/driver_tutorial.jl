# # Using atmospheric and radiative drivers
# The goal of this is to outline how to set up simulations driven by
# prescribed forcing data (``drivers"). These are grouped into
# radiative forcing and atmospheric forcing. We will first cover
# the types of forcing we support, followed by how to specify the
# driver structs given the forcing data and how to update the values
# used during a simulation.

# # Types of forcing data

# We currently support site-level simulations and have two site-level
# driver types, `PrescribedAtmosphere` and `PrescribedRadiativeFluxes`.

# The atmosphere driver stores the atmospheric state data as a
# function of time, including the liquid precipitation rate (m/s), the
# snow precipitation rate converted into an equivalent rate of liquid
# water (m/s), the atmopheric pressure (Pa), specific humidity, horizontal
# wind speed (m/s), temperature (K), CO2 concentration (mol/mol), and the
# height at which these measurements were taken (currently assumed to be the
# same value for all variables).

# The radiative fluxes driver stores the data required to specify the
# radiative forcing. We currently support only a single downwelling
# shortwave and longwave flux (W/m^2). The radiative driver is also where
# a function which computes the zenith angle for the site is stored.

# Both drivers store the start date for the data/simulation.
# This is the DateTime object which corresponds to the time at which t=0
# in the simulation. Additionally, for site-level runs, both drivers store the
# forcing data as a spline function fit to the data which takes the time
# `t` as an argument, where `t` is the simulation time measured in seconds since
# the start date. The start date should be in UTC.

# Note: for coupled runs, corresponding types `CoupledAtmosphere`
# and `CoupledRadiativeFluxes` exist. However, these are not defined
# in ClimaLand, but rather inside of the Clima Coupler repository.


# # Creating site-level drivers for radiation

# First, assume that we have data stored for the longwave and shortwave radiation
# at a particular site, and that we have read it in to an array, along with the
# times at which the observations were made and the latitude and longitude of the site.
using Dates
using Insolation # for computing zenith angle given lat, lon, time.
using ClimaLand
import ClimaLand.Parameters as LP
import ClimaParams

# Assume the local_datetime array is read in from the data file.
local_datetime = DateTime(2013):Dates.Hour(1):DateTime(2013, 1, 7); # one week, hourly data
# Timezone (offset of local time from UTC in hrs) for the Missouri Ozark site
time_offset = -6;
# Site latitude and longitude
lat = 38.7441; # degree
long = -92.2000; # degree
# Compute the start date in UTC, and convert local datetime
# vector into a vector of seconds since the start date
start_date = local_datetime[1] + Dates.Hour(time_offset);
data_dt = 3600.0;
seconds = 0:data_dt:((length(local_datetime) - 1) * data_dt);

# Assume the downwelling long and shortwave radiation, as well as the diffuse fraction
# are read in from the file
# and are measured at the times in local_datetime. Here, we'll just make them up
# periodic on daily timescales:
T = @. 298.15 + 5.0 * sin(2π * (seconds - 3600 * 6) / (3600 * 24));
LW_d = 5.67 * 10^(-8) .* T .^ 4;
SW_d = @. max(1400 * sin(2π * (seconds - 3600 * 6) / (3600 * 24)), 0.0);
diffuse_fraction = 0.5 .+ zeros(length(seconds));

# Next, fit interpolators to the data. These interpolators are what are stored in
# the driver function. Then we can evaluate the radiative forcing
# at any simulation time (and not just at times coinciding with measurements).
# By default, linear interpolation is used.
LW_d = TimeVaryingInput(seconds, LW_d)
SW_d = TimeVaryingInput(seconds, SW_d);
diffuse_fraction = TimeVaryingInput(seconds, diffuse_fraction);

# Finally, for many models we also need to specify the function
# for computing the zenith angle as a function of simulation time.
# To do so, we use the `Insolation` package as follows:
toml_dict = LP.create_toml_dict(Float64);
earth_param_set = LP.LandParameters(toml_dict);
insol_params = earth_param_set.insol_params # parameters of Earth's orbit required to compute the insolation
zenith_angle =
    (t, s) -> default_zenith_angle(
        t,
        s;
        insol_params,
        longitude = long,
        latitude = lat,
    );

# Lastly, we store the interpolators for downwelling fluxes and the zenith angle function
# in the `PrescribedRadiativeFluxes` struct.
radiation = ClimaLand.PrescribedRadiativeFluxes(
    Float64,
    SW_d,
    LW_d,
    start_date;
    θs = zenith_angle,
    frac_diff = diffuse_fraction,
);

# # Updating the driver variables during the simulation
# The values for LW_d, SW_d, the diffuse fraction,
# and the cosine of the zenith angle θ_s are stored
# in the simulation/model cache `p` under the name `drivers`.
# Since we did not pass the diffuse fraction directly,
# it will be computed using an empirical function
# of the atmospheric state, the zenith angle, and the earth
# parameter set.

# When you `initialize` the variables and cache of a model,
# the cache `p` will be returned with memory allocated but all
# values set to zero:
p = (; drivers = (LW_d = [0.0], SW_d = [0.0], cosθs = [0.0], frac_diff = [0.0]));

# In order to update them, we can make use of default update functions:
update_radiation! = ClimaLand.make_update_drivers(radiation)
t0 = seconds[1] # midnight local time
update_radiation!(p, t0);
@show(p.drivers);

# During a simulation, the drivers are updated in place in `p.drivers`
# via a "callback", which is a function which is called a specified times or when
# certain criteria are met during a simulation.
# In general, then, we don't update drivers every timestep, but less frequently.
# For example, the simulation timestep may be 10 minutes, but we may only update
# the drivers every three hours:
updatefunc = update_radiation!;
cb = ClimaLand.DriverUpdateCallback(updatefunc, 3600.0 * 3, t0);

# This callback must then be provided to the simulation [`solve`](https://docs.sciml.ai/DiffEqCallbacks/stable/) function.

# # Using ERA5 data
# If you wish to force your ClimaLand simulation with ERA5 reanalysis data, there is a helper function
# makes this easier than specifying each atmospheric variable as individual `TimeVaryingInput`
# objects, and then making the `PrescribedAtmosphere` struct.
# ClimaLand provides two ERA5 NetCDF datasets:
# - A high-resolution dataset with data on a 1 degree x 1 degree grid, available from 1979 to 2024.
#  This dataset is about 8.5 GB per year of data.
# - A low-resolution dataset with data on a 8 degree x 8 degree grid, available for the year 2008.
#  This dataset is about 350 MB in size.
# The low-resolution dataset is best for local simulations and for use in examples
# in the documentation, but for production runs on compute clusters we recommend using the
# high-resolution data.

# To use the ERA5 data, you need to provide the `start_date` and `stop_date` of your simulation,
# which should be within the range of dates covered by the dataset you are using (unless you are using
# the low-resolution dataset for 2008, which will be reused for each year of simulation).
# You'll also provide a flag `use_lowres_forcing` which is true if you want to use the low-resolution
# dataset, and false otherwise.
# Finally, you need the `surface_space` of your simulation (corresponding to the grid being used),
# the parameter `toml_dict`, and the floating point type of the simulation `FT`.
# Then you can access the atmospheric and radiative drivers like this:

# ```julia
# atmos, radiation = ClimaLand.prescribed_forcing_era5(start_date,
#                                                      stop_date,
#                                                      surface_space,
#                                                      earth_param_set,
#                                                      FT;
#                                                      use_lowres_forcing = true)
# ```

# Each forcing dataset contains the following variables:
# - "tp" Total precipitation as a mass/m^2/hour (accumulated over an hour)
# - "sf" Snow precipitation as a mass/m^2/hour (accumulated over an hour)
# - "u10n", "v10n", Neutral wind speed components in the horizontal, at 10m, in m/s
# - "d2m", Dewpoint temperature at 2m in K
# - "t2m", Air temperature at 2m in K
# - "sp", Surface pressure in Pa
# - "ssrd", Downwelling shortwave radiation in J/m^2/hour (accumulated over an hour)
# - "strd", Downwelling longwave radiation J/m^2/hour (accumulated over an hour)

# Please note that the default era5 forcing uses linear interpolation in space and time.
# If your simulation encompasses a time that is beyond the extrema of the data,
# the corresponding day/time from the last year of data will be repeated.
# This corresponds to a value of `time_interpolation_method = LinearInterpolation(PeriodicCalendar())
# `time_interpolation_method = LinearInterpolation(PeriodicCalendar(Dates.Year(1), DateTime(Dates.year(stop_date))))`.
# To repeat the first year of data instead, you can specify to repeat the start date's year, i.e.:
# `time_interpolation_method = LinearInterpolation(PeriodicCalendar(Dates.Year(1), DateTime(Dates.year(start_date))))`.
# Or, to repeat the entire time series of data you can call `PeriodicCalendar` without any arguments:
# `time_interpolation_method = LinearInterpolation(PeriodicCalendar())`.
# This behavior can be changed by passing in a time_interpolation_method. Another
# option that be be useful is `time_interpolation_method = LinearInterpolation(Flat())`,
# which uses the last (first) value of the data repeatedly when outside the bounds
# of the data. For more details, please see the documentation string for
# ` ClimaLand.prescribed_forcing_era5` or ClimaUtilities.
