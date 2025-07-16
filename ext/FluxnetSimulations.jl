module FluxnetSimulations
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput,
    LinearInterpolation,
    PeriodicCalendar
import ClimaUtilities.TimeManager: ITime, date
using Thermodynamics
using Dates
using DelimitedFiles
using DocStringExtensions
using Insolation
import ClimaLand.Parameters as LP
using ClimaLand
export prescribed_forcing_fluxnet
"""
     prescribed_forcing_fluxnet(site_ID,
                                lat,
                                long,
                                hour_offset_from_UTC,
                                start_date,
                                earth_param_set,
                                FT;
                                gustiness=1,
                                c_co2 = TimeVaryingInput((t) -> 4.2e-4),
                                time_interpolation_method = LinearInterpolation(PeriodicCalendar())) # not used??

A helper function which constructs the `PrescribedAtmosphere` and `PrescribedRadiativeFluxes`
from a file path pointing to the Fluxnet data in a csv file, the start date, latitude, longitude,
the hour offset of the site from UTC (local_time + offset = time in UTC),
and the earth_param_set.
"""
function prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    hour_offset_from_UTC,
    atmos_h,
    start_date, # in UTC
    earth_param_set,
    FT;
    split_precip = true,
    gustiness = 1,
    c_co2 = TimeVaryingInput((t) -> 4.2e-4),
    time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
)
    fluxnet_csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID);
    driver_data = readdlm(fluxnet_csv_path, ',')
    columns = driver_data[1,:]
    labels = ("TA_F", "VPD_F", "PA_F", "P_F", "WS_F", "LW_IN_F", "SW_IN_F")
    indices = [findfirst(columns .== label) for label in labels]

    if any(indices .== nothing)
        nothing_id = findall(indices .== nothing)
        @error("$(labels[nothing_id]) is missing in the data, but required.")
    end
    
    index_map = Dict(zip(labels, indices))
    LOCAL_DATETIME = DateTime.(format.(driver_data[2:end, 1]), "yyyymmddHHMM")
    UTC_DATETIME = LOCAL_DATETIME .+ Dates.Hour(hour_offset_from_UTC)
    DATA_DT = Second(LOCAL_DATETIME[2] - LOCAL_DATETIME[1]).value # seconds
    seconds = [Second(UTC - start_date).value for UTC in UTC_DATETIME]
    
    atmos_T_data = driver_data[2:end, index_map["TA_F"]];
    not_missing_mask = .~ ismissing.(atmos_T_data)
    atmos_T = TimeVaryingInput(seconds[not_missing_mask], Float64.(atmos_T_data[not_missing_mask]) .+ 273.15)

    atmos_P_data = driver_data[2:end, index_map["PA_F"]];
    not_missing_mask = .~ ismissing.(atmos_P_data)
    atmos_P = TimeVaryingInput(seconds[not_missing_mask], Float64.(atmos_P_data[not_missing_mask]) .* 1000)
    
    atmos_u_data = driver_data[2:end, index_map["WS_F"]];
    not_missing_mask = .~ ismissing.(atmos_u_data)
    atmos_u = TimeVaryingInput(seconds[not_missing_mask], Float64.(atmos_u_data[not_missing_mask]))

    LW_d_data = driver_data[2:end, index_map["LW_IN_F"]];
    not_missing_mask = .~ ismissing.(LW_d_data)
    LW_d = TimeVaryingInput(seconds[not_missing_mask], Float64.(LW_d_data[not_missing_mask]))

    SW_d_data = driver_data[2:end, index_map["SW_IN_F"]];
    not_missing_mask = .~ ismissing.(SW_d_data)
    SW_d = TimeVaryingInput(seconds[not_missing_mask], Float64.(SW_d_data[not_missing_mask]))

    atmos_VPD_data = driver_data[2:end, index_map["VPD_F"]];

    # For q, we need to compute the values first
    not_missing_mask = .~ ismissing.(atmos_P_data) .&& .~ ismissing(atmos_VPD_data) .&& .~ ismissing(atmos_T_data)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    esat =
        Thermodynamics.saturation_vapor_pressure.(
            Ref(thermo_params),
            atmos_T_data[not_missing_mask],
            Ref(Thermodynamics.Liquid()),
        )
    e = @. esat - (atmos_VPD_data[not_missing_mask] .* 1000) # unit conversion of VPD to Pa
    q = @. 0.622 * e ./ (atmos_P_data[not_missing_mask] - 0.378 * e)

    atmos_q = TimeVaryingInput(seconds[not_missing_mask], Float64.(q))

    # Split precip into rain and snow
    atmos_precip_data = driver_data[2:end, index_map["P_F"]];
    if split_precip
        not_missing_mask = .~ ismissing.(atmos_T_data) .&& .~ ismissing(atmos_precip_data) .&& .~ ismissing(atmos_VPD_data)
        esat =
            Thermodynamics.saturation_vapor_pressure.(
                Ref(thermo_params),
                atmos_T_data[not_missing_mask],
                Ref(Thermodynamics.Liquid()),
            )
        e = @. esat - (atmos_VPD_data[not_missing_mask] .* 1000) # unit conversion of VPD to Pa
        RH = @. e / esat
        snow_frac = snow_precip_fraction.(atmos_T_data[not_missing_mask], RH)
        P_liq_data = -1 .* Float64.(atmos_precip_data[not_missing_mask]) .* (1 .- snow_frac)
        P_snow_data = -1 .* Float64.(atmos_precip_data[not_missing_mask]) .* snow_frac
        P_liq = TimeVaryingInput(seconds[not_missing_mask], P_liq_data ./ 1000 ./ DATA_DT)
        P_snow = TimeVaryingInput(seconds[not_missing_mask], P_snow_data ./ 1000 ./ DATA_DT)
    else
        not_missing_mask = .~ ismissing(atmos_precip_data)
        P_liq_data = -1 .* Float64.(atmos_precip_data[not_missing_mask])
        P_liq = TimeVaryingInput(seconds[not_missing_mask], P_liq_data ./ 1000 ./ DATA_DT)
        P_snow = TimeVaryingInput(seconds[not_missing_mask],zeros(length(seconds[not_missing_mask])))
    end

# Construct the drivers. For the start date we will use the UTC time at the
# start of the simulation
atmos = ClimaLand.PrescribedAtmosphere(
    P_liq,
    P_snow,
    atmos_T,
    atmos_u,
    atmos_q,
    atmos_P,
    start_date,
    atmos_h,
    earth_param_set;
    c_co2
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
    start_date,
    θs = zenith_angle,
    earth_param_set = earth_param_set,
)

    # Desired start and end dates of data in MODIS format
    modis_start_date = DateTime(
        "$(Dates.year(start_date))-$(Dates.month(start_date))-$(Dates.day(start_date))T00:00:00.0",
    )
    modis_end_date = DateTime(
        "$(Dates.year(UTC_DATETIME[end]))-$(Dates.month(UTC_DATETIME[end]))-$(Dates.day(UTC_DATETIME[end]))T00:00:00.0",
    )

    MODIS_LAI_path = ClimaLand.Artifacts.get_modis_lai_fluxnet_data(site_ID)
    MODIS_LAI_raw, header = readdlm(MODIS_LAI_path, ',', header = true)
    dates = DateTime.(MODIS_LAI_raw[:, 1])
    
    # Trim MODIS_LAI to start_date and end_date - this shouldnt be necessary
    indices = findall(d -> modis_start_date <= d <= modis_end_date, dates)
    filtered_rows = MODIS_LAI_raw[indices, :]
    filtered_dates = dates[indices]
    
    # Create final matrix with properly formatted data for computation
    MODIS_LAI = hcat(filtered_dates, filtered_rows[:, 2])
    
    LAI_dt = Second(MODIS_LAI[2, 1] - MODIS_LAI[1, 1]).value
    LAI_seconds = FT.(0:LAI_dt:((length(MODIS_LAI[:, 1]) - 1) * LAI_dt))
    
    # LAI function for radiative transfer
    LAIfunction = TimeVaryingInput(LAI_seconds, FT.(MODIS_LAI[:, 2]))
    maxLAI = FT(maximum(MODIS_LAI[:, 2]))

    return (; atmos, radiation, LAI, maxLAI)
end

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
end
