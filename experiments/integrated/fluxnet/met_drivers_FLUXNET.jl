using DelimitedFiles
using Thermodynamics
using Dates
using Format
using HTTP
using JSON
using Insolation

using ClimaLand

import ClimaComms

context = ClimaComms.context()

data_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)
driver_data = readdlm(data_path, ',')

LOCAL_DATETIME = DateTime.(format.(driver_data[2:end, 1]), "yyyymmddHHMM")
UTC_DATETIME = LOCAL_DATETIME .+ Dates.Hour(time_offset)
DATA_DT = Second(LOCAL_DATETIME[2] - LOCAL_DATETIME[1]).value # seconds

# Label of every data column to be collected from the driver data file
labels = (
    :TA,
    :VPD,
    :PA,
    :P,
    :WS,
    :LW_IN,
    :SW_IN,
    :G,
    :GPP,
    :LE,
    :H,
    :LW_OUT,
    :SW_OUT,
    :SWC,
    :TS,
    :CO2,
)

# For every data column to be collected, name of the column in the file,
# transformation function to desired unit, and string of final unit
collect_args = [
    ("TA_F", (x) -> x .+ 273.15, "K")
    ("VPD_F", (x) -> x .* 100, "Pa")
    ("PA_F", (x) -> x .* 1000, "Pa")
    ("P_F", (x) -> x ./ (1000 * DATA_DT), "m/s")
    ("WS_F", (x) -> x, "m/s")
    ("LW_IN_F", (x) -> x, "W/m^2")
    ("SW_IN_F", (x) -> x, "W/m^2")
    ("G_F_MDS", (x) -> x, "W/m^2")
    ("GPP_DT_VUT_REF", (x) -> x .* 1e-6, "mol/m^2/s")
    ("LE_CORR", (x) -> x, "W/m^2")
    ("H_CORR", (x) -> x, "W/m^2")
    ("LW_OUT", (x) -> x, "W/m^2")
    ("SW_OUT", (x) -> x, "W/m^2")
    ("SWC_F_MDS_1", (x) -> x ./ 100, "m^3/m^3")
    ("TS_F_MDS_1", (x) -> x .+ 273.15, "K")
    ("CO2_F_MDS", (x) -> x .* 1e-6, "mol")
]

# Named tuple mapping every label to a DataColumn with the correct transformed
# unit and data status. Automatically fills gaps in the data
drivers = (;
    zip(
        labels,
        [
            transform_column(
                filter_column(driver_data, name, ""),
                transform,
                unit,
            ) for (name, transform, unit) in collect_args
        ],
    )...
)
# Check that all required driver data is present, otherwise throw error
required = [
    "TA" => drivers.TA,
    "VPD" => drivers.VPD,
    "PA" => drivers.PA,
    "P" => drivers.P,
    "WS" => drivers.WS,
    "LW_IN" => drivers.LW_IN,
    "SW_IN" => drivers.SW_IN,
]
missing_drivers = []
for (name, column) in required
    if column.status == absent
        push!(missing_drivers, name)
    end
end
if length(missing_drivers) != 0
    error("Driver data missing for columns: $([missing_drivers[i] * " " for
        i in 1:length(missing_drivers)]...)")
end

thermo_params = LP.thermodynamic_parameters(earth_param_set)
esat =
    Thermodynamics.saturation_vapor_pressure.(
        Ref(thermo_params),
        FT.(drivers.TA.values),
        Ref(Thermodynamics.Liquid()),
    )
e = @. esat - drivers.VPD.values
q = @. 0.622 * e ./ (drivers.PA.values - 0.378 * e)
RH = @. e / esat

"""
    snow_precip_fraction(air_temp, hum)

Estimate the fraction of precipitation that is in snow form,
given the air temperature at the surface in K and the relative humidity
(between 0 and 1).

See Jennings, K.S., Winchell, T.S., Livneh, B. et al.
Spatial variation of the rain–snow temperature threshold across the
Northern Hemisphere. Nat Commun 9, 1148 (2018).
https://doi.org/10.1038/s41467-018-03629-7
"""
function snow_precip_fraction(air_temp, hum)
    air_temp_C = air_temp - 273.15
    α = -10.04
    β = 1.41
    γ = 0.09
    snow_frac = (1.0 / (1.0 + exp(α + β * air_temp_C + γ * hum)))
    return snow_frac
end
snow_frac = snow_precip_fraction.(drivers.TA.values[:], RH[:])
# Create interpolators for each atmospheric driver needed for PrescribedAtmosphere and for
# PrescribedRadiation
seconds = FT.(0:DATA_DT:((length(UTC_DATETIME) - 1) * DATA_DT));

P_liq = -FT.(drivers.P.values[:] .* (1 .- snow_frac))
precip = TimeVaryingInput(seconds, P_liq; context) # m/s
atmos_q = TimeVaryingInput(seconds, FT.(q[:]); context)
atmos_T = TimeVaryingInput(seconds, FT.(drivers.TA.values[:]); context)
atmos_p = TimeVaryingInput(seconds, FT.(drivers.PA.values[:]); context)
atmos_co2 = TimeVaryingInput(seconds, FT.(drivers.CO2.values[:]); context)
atmos_u = TimeVaryingInput(seconds, FT.(drivers.WS.values[:]); context)
LW_IN = TimeVaryingInput(seconds, FT.(drivers.LW_IN.values[:]); context)
SW_IN = TimeVaryingInput(seconds, FT.(drivers.SW_IN.values[:]); context)
P_snow = -FT.(drivers.P.values[:] .* snow_frac)
snow_precip = TimeVaryingInput(seconds, P_snow; context) # m/s

# Construct the drivers. For the start date we will use the UTC time at the
# start of the simulation
atmos = ClimaLand.PrescribedAtmosphere(
    precip,
    snow_precip,
    atmos_T,
    atmos_u,
    atmos_q,
    atmos_p,
    UTC_DATETIME[1],
    atmos_h,
    earth_param_set;
    c_co2 = atmos_co2,
)

zenith_angle =
    (t, s) -> default_zenith_angle(
        t,
        s;
        insol_params = earth_param_set.insol_params,
        longitude = long,
        latitude = lat,
    )
radiation = ClimaLand.PrescribedRadiativeFluxes(
    FT,
    SW_IN,
    LW_IN,
    UTC_DATETIME[1];
    θs = zenith_angle,
    earth_param_set = earth_param_set,
)

# Desired start and end dates of data in MODIS format
modis_start_date = DateTime(
    "$(Dates.year(UTC_DATETIME[1]))-$(Dates.month(UTC_DATETIME[1]))-$(Dates.day(UTC_DATETIME[1]))T00:00:00.0",
)
modis_end_date = DateTime(
    "$(Dates.year(UTC_DATETIME[end]))-$(Dates.month(UTC_DATETIME[end]))-$(Dates.day(UTC_DATETIME[end]))T00:00:00.0",
)

MODIS_LAI_path = ClimaLand.Artifacts.get_modis_lai_fluxnet_data(site_ID)
MODIS_LAI_raw, header = readdlm(MODIS_LAI_path, ',', header = true)
dates = DateTime.(MODIS_LAI_raw[:, 1])

# Trim MODIS_LAI to start_date and end_date
indices = findall(d -> modis_start_date <= d <= modis_end_date, dates)
filtered_rows = MODIS_LAI_raw[indices, :]
filtered_dates = dates[indices]

# Create final matrix with properly formatted data for computation
MODIS_LAI = hcat(filtered_dates, filtered_rows[:, 2])

LAI_dt = Second(MODIS_LAI[2, 1] - MODIS_LAI[1, 1]).value
LAI_seconds = FT.(0:LAI_dt:((length(MODIS_LAI[:, 1]) - 1) * LAI_dt))

# LAI function for radiative transfer
LAIfunction = TimeVaryingInput(LAI_seconds, FT.(MODIS_LAI[:, 2]); context)

# Necessary inputs from LAI Data
# Note that f_root_to_shoot, capacity, and h_leaf are site-specific parameters
# defined in the parameters file for each site
maxLAI = FT(maximum(MODIS_LAI[:, 2]))
# RAI = maxLAI * f_root_to_shoot
# capacity = plant_ν * maxLAI * h_leaf * FT(1000)
