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

# Both drivers store the reference time for the data/simulation.
# This is the DateTime object which corresponds to the time at which t=0
# in the simulation. Additionally, for site-level runs, both drivers store the
# forcing data as a spline function fit to the data which takes the time
# `t` as an argument, where `t` is the simulation time measured in seconds since
# the reference time. The reference time should be in UTC.

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
include(joinpath(pkgdir(ClimaLand), "parameters", "create_parameters.jl"));

# Assume the local_datetime array is read in from the data file.
local_datetime = DateTime(2013):Dates.Hour(1):DateTime(2013, 1, 7); # one week, hourly data
# Timezone (offset of local time from UTC in hrs)
time_offset = 7;
# Site latitude and longitude
lat = 38.7441; # degree
long = -92.2000; # degree
# Compute the reference time in UTC, and convert local datetime
# vector into a vector of seconds since the reference time
ref_time = local_datetime[1] + Dates.Hour(time_offset);
data_dt = 3600.0;
seconds = 0:data_dt:((length(local_datetime) - 1) * data_dt);

# Assume the downwelling long and shortwave radiation are read in from the file
# and are measured at the times in local_datetime. Here, we'll just make them up
# periodic on daily timescales:
T = @. 298.15 + 5.0 * sin(2π * (seconds - 3600 * 6) / (3600 * 24));
LW_d = 5.67 * 10^(-8) .* T .^ 4;
SW_d = @. max(1400 * sin(2π * (seconds - 3600 * 6) / (3600 * 24)), 0.0);

# Next, fit interpolators to the data. These interpolators are what are stored in
# the driver function. Then we can evaluate the radiative forcing
# at any simulation time (and not just at times coinciding with measurements).
# By default, linear interpolation is used.
LW_d = TimeVaryingInput(seconds, LW_d)
SW_d = TimeVaryingInput(seconds, SW_d);

# Finally, for many models we also need to specify the function
# for computing the zenith angle as a function of simulation time.
# To do so, we use the `Insolation` package as follows:
earth_param_set = create_lsm_parameters(Float64);
insol_params = earth_param_set.insol_params # parameters of Earth's orbit required to compute the insolation
function zenith_angle(
    t,
    ref_time;
    latitude = lat,
    longitude = long,
    insol_params = insol_params,
)
    current_datetime = ref_time + Dates.Second(round(t)) # Time in UTC

    d, δ, η_UTC = (Insolation.helper_instantaneous_zenith_angle(
        current_datetime,
        ref_time,
        insol_params,
    ))


    return Insolation.instantaneous_zenith_angle(
        d,
        δ,
        η_UTC,
        longitude,
        latitude,
    )[1]
end;

# Lastly, we store the interpolators for downwelling fluxes and the zenith angle function
# in the `PrescribedRadiativeFluxes` struct.
radiation = ClimaLand.PrescribedRadiativeFluxes(
    Float64,
    SW_d,
    LW_d,
    ref_time;
    θs = zenith_angle,
);

# # Updating the driver variables during the simulation
# The values for LW_d, SW_d, and zenith angle θ_s are stored
# in the simulation/model cache `p` under the name `drivers`.
# When you `initialize` the variables and cache of a model,
# the cache `p` will be returned with memory allocated but all
# values set to zero:
p = (; drivers = (LW_d = [0.0], SW_d = [0.0], θs = [0.0]));

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
updateat = collect(seconds[1]:(3600 * 3):seconds[end]);
updatefunc = update_radiation!;
cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc);

# This callback must then be provided to the simulation [`solve`](https://docs.sciml.ai/DiffEqCallbacks/stable/) function.
