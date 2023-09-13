#Citation: Dave Billesbach, Ryan Sullivan (2021), AmeriFlux BASE US-A10 ARM-NSA-Barrow, Ver. 4-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1498753
#Atmospheric Radiation Measurement (ARM) user facility. 1998. ARM Best Estimate Data Products (ARMBECLDRAD). 2011-01-01 to2013-12-31
#North Slope Alaska (NSA), Central Facility, Barrow AK (C1)
# Compiled by C. Xiao and X. Shaocheng. ARM Data Center. Data set accessed http://dx.doi.org/10.5439/1333228.
using DelimitedFiles
using Dierckx
using Thermodynamics
using Dates
using NCDatasets
function replace_missing_with_mean!(field, flag)
    good_indices = (flag .== 0) .|| (flag .== 1)
    fill_value = mean(field[good_indices])
    field[.~good_indices] .= fill_value
    return field
end


function replace_missing_with_mean_by_value!(field)
    good_indices = .~(field .== -9999)
    fill_value = mean(field[good_indices])
    field[.~good_indices] .= fill_value
    return field
end

function replace_missing_with_zero_by_value!(field)
    good_indices = .~(field .== -9999)
    field[.~good_indices] .= 0.0
    return field
end

# Use ARM data when possible
arm_data = NCDataset(
    "/Users/katherinedeck/Desktop/Barrow/nsaarmbeatmC1.c1.20130101.003000.custom.nc",
)
WS = @. sqrt(arm_data["u_wind_sfc"][:]^2 + arm_data["v_wind_sfc"][:]^2)
# wind speed is at 10m, other measurements at 2m
# linearly interpolate to ground at zero speed
WS_2m = (2 / 10) .* WS
WS_2m[typeof.(WS_2m) .== Missing] .= mean(WS_2m[.~(typeof.(WS_2m) .== Missing)])

arm_LOCAL_DATETIME = arm_data["time"][:]
arm_UTC_DATETIME = arm_LOCAL_DATETIME .+ Dates.Hour(8)
TA = arm_data["temperature_sfc"][:]
TA[typeof.(TA) .== Missing] .= mean(TA[.~(typeof.(TA) .== Missing)])
P = arm_data["precip_rate_sfc"][:] ./ (1000 * 3600) # convert mm/hr to m/s
P[typeof.(P) .== Missing] .= 0.0f0
PA = arm_data["pressure_sfc"][:] .* 100
PA[typeof.(PA) .== Missing] .= mean(PA[.~(typeof.(PA) .== Missing)])
RH = arm_data["relative_humidity_sfc"][:]
RH[typeof.(RH) .== Missing] .= mean(RH[.~(typeof.(RH) .== Missing)])
thermo_params = LSMP.thermodynamic_parameters(earth_param_set)
esat =
    Thermodynamics.saturation_vapor_pressure.(
        Ref(thermo_params),
        TA,
        Ref(Thermodynamics.Liquid()),
    )
e = @. RH * esat / 100
q = @. 0.622 * e ./ (PA - 0.378 * e)
arm_DATA_DT = 3600


#Make a bunch of splines
arm_seconds = FT.(1800:arm_DATA_DT:(length(arm_UTC_DATETIME) * arm_DATA_DT));
p_spline = Spline1D(arm_seconds, -P[:]) # m/s
atmos_q = Spline1D(arm_seconds, q[:])
atmos_T = Spline1D(arm_seconds, TA[:])
atmos_p = Spline1D(arm_seconds, PA[:])
atmos_u = Spline1D(arm_seconds, WS_2m[:])
atmos_h = FT(2)
precipitation_function(t::FT) where {FT} = p_spline(t) < 0.0 ? p_spline(t) : 0.0 # m/s
snow_precip(t) = eltype(t)(0) # this is likely not correct


# Get other data from Fluxnet
dataset_path = "/Users/katherinedeck/Desktop/Barrow/AMF_US-A10_BASE-BADM_4-5"
data = joinpath(dataset_path, "AMF_US-A10_BASE_HH_4-5.csv");
driver_data = readdlm(data, ',')
column_names = driver_data[3, :]

indices = 1:1:length(driver_data[4:end, 1])
LOCAL_DATETIME = DateTime.(string.(driver_data[4:end, 1]), "yyyymmddHHMM")
subset = indices[Dates.Year.(LOCAL_DATETIME) .== Dates.Year(2013)]
LOCAL_DATETIME = LOCAL_DATETIME[subset]
UTC_DATETIME = LOCAL_DATETIME .+ Dates.Hour(8)
DATA_DT = Second(LOCAL_DATETIME[2] - LOCAL_DATETIME[1]).value # seconds

CO2 = driver_data[2:end, column_names .== "CO2"][subset]
replace_missing_with_mean_by_value!(CO2)
CO2 .= CO2 .* 1e-6 # μmol to mol

LW_IN = driver_data[2:end, column_names .== "LW_IN"][subset]
replace_missing_with_mean_by_value!(LW_IN)
SW_IN = driver_data[2:end, column_names .== "SW_IN"][subset]
replace_missing_with_zero_by_value!(SW_IN)


SWC_1 = driver_data[2:end, column_names .== "SWC_PI_1_1_A"][subset]
SWC_2 = driver_data[2:end, column_names .== "SWC_PI_1_2_A"][subset]
replace_missing_with_zero_by_value!(SWC_1)
replace_missing_with_zero_by_value!(SWC_2)
SWC_1 .= SWC_1 ./ 100
SWC_2 .= SWC_2 ./ 100

TS_1 = driver_data[2:end, column_names .== "TS_PI_1_1_A"][subset]
TS_2 = driver_data[2:end, column_names .== "TS_PI_2_2_A"][subset]
TS_3 = driver_data[2:end, column_names .== "TS_PI_2_3_A"][subset]
replace_missing_with_zero_by_value!(TS_1)
replace_missing_with_zero_by_value!(TS_2)
replace_missing_with_zero_by_value!(TS_3)

TS_1 .= TS_1 .+ 273.15;# convert C to K
TS_2 .= TS_2 .+ 273.15;# convert C to K
TS_3 .= TS_3 .+ 273.15;# convert C to K

LE = driver_data[2:end, column_names .== "LE"][subset]
H = driver_data[2:end, column_names .== "H"][subset]
G = driver_data[2:end, column_names .== "G_PI_1_1_A"][subset]
replace_missing_with_zero_by_value!(LE)
replace_missing_with_zero_by_value!(H)
replace_missing_with_zero_by_value!(G)


LW_OUT = driver_data[2:end, column_names .== "LW_OUT"][subset]
SW_OUT = driver_data[2:end, column_names .== "SW_OUT"][subset]
replace_missing_with_zero_by_value!(SW_OUT)
replace_missing_with_zero_by_value!(LW_OUT)


#Make a bunch of splines
seconds = FT.(0:DATA_DT:((length(UTC_DATETIME) - 1) * DATA_DT));
atmos_co2 = Spline1D(seconds, CO2[:])
LW_IN_spline = Spline1D(seconds, LW_IN[:])
SW_IN_spline = Spline1D(seconds, SW_IN[:])

# Construct the drivers. For the reference time we will use the UTC time at the
# start of the simulation
atmos = ClimaLSM.PrescribedAtmosphere(
    precipitation_function,
    snow_precip,
    atmos_T,
    atmos_u,
    atmos_q,
    atmos_p,
    UTC_DATETIME[1],
    atmos_h;
    c_co2 = atmos_co2,
)

lat = FT(72.32193) # degree
long = FT(-156.608) # degree

function zenith_angle(
    t,
    orbital_data,
    ref_time;
    latitude = lat,
    longitude = long,
    insol_params = earth_param_set.insol_params,
)
    FT = Float64
    # This should be time in UTC
    dt = ref_time + Dates.Second(round(t))
    FT(
        instantaneous_zenith_angle(
            dt,
            orbital_data,
            longitude,
            latitude,
            insol_params,
        )[1],
    )
end


radiation = ClimaLSM.PrescribedRadiativeFluxes(
    FT,
    SW_IN_spline,
    LW_IN_spline,
    UTC_DATETIME[1];
    θs = zenith_angle,
    orbital_data = Insolation.OrbitalData(),
)

LAIfunction = (t) -> 0.0
