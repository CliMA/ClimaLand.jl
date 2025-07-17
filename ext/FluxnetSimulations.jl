module FluxnetSimulations
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.TimeManager: ITime, date
using Thermodynamics
using Dates
using DelimitedFiles
using Format
using DocStringExtensions
using Insolation
import ClimaLand.Parameters as LP
using ClimaLand
export prescribed_forcing_fluxnet, set_fluxnet_ic!
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
    fluxnet_csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)
    driver_data = readdlm(fluxnet_csv_path, ',')
    columns = driver_data[1, :]
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
    seconds_since_start_date =
        [Second(UTC - start_date).value for UTC in UTC_DATETIME]

    atmos_T_data = driver_data[2:end, index_map["TA_F"]]
    not_missing_mask = .~ismissing.(atmos_T_data)
    atmos_T = TimeVaryingInput(
        seconds_since_start_date[not_missing_mask],
        Float64.(atmos_T_data[not_missing_mask]) .+ 273.15,
    )

    atmos_P_data = driver_data[2:end, index_map["PA_F"]]
    not_missing_mask = .~ismissing.(atmos_P_data)
    atmos_P = TimeVaryingInput(
        seconds_since_start_date[not_missing_mask],
        Float64.(atmos_P_data[not_missing_mask]) .* 1000,
    )

    atmos_u_data = driver_data[2:end, index_map["WS_F"]]
    not_missing_mask = .~ismissing.(atmos_u_data)
    atmos_u = TimeVaryingInput(
        seconds_since_start_date[not_missing_mask],
        Float64.(atmos_u_data[not_missing_mask]),
    )

    LW_d_data = driver_data[2:end, index_map["LW_IN_F"]]
    not_missing_mask = .~ismissing.(LW_d_data)
    LW_d = TimeVaryingInput(
        seconds_since_start_date[not_missing_mask],
        Float64.(LW_d_data[not_missing_mask]),
    )

    SW_d_data = driver_data[2:end, index_map["SW_IN_F"]]
    not_missing_mask = .~ismissing.(SW_d_data)
    SW_d = TimeVaryingInput(
        seconds_since_start_date[not_missing_mask],
        Float64.(SW_d_data[not_missing_mask]),
    )

    atmos_VPD_data = driver_data[2:end, index_map["VPD_F"]]

    # For q, we need to compute the values first
    not_missing_mask =
        .~ismissing.(atmos_P_data) .&&
        .~ismissing(atmos_VPD_data) .&&
        .~ismissing(atmos_T_data)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    esat =
        Thermodynamics.saturation_vapor_pressure.(
            Ref(thermo_params),
            atmos_T_data[not_missing_mask] .+ 273.15,
            Ref(Thermodynamics.Liquid()),
        )
    e = @. esat - (atmos_VPD_data[not_missing_mask] .* 1000) # unit conversion of VPD to Pa
    q = @. 0.622 * e ./ (atmos_P_data[not_missing_mask] - 0.378 * e)

    atmos_q = TimeVaryingInput(
        seconds_since_start_date[not_missing_mask],
        Float64.(q),
    )

    # Split precip into rain and snow
    atmos_precip_data = driver_data[2:end, index_map["P_F"]]
    if split_precip
        not_missing_mask =
            .~ismissing.(atmos_T_data) .&&
            .~ismissing.(atmos_precip_data) .&&
            .~ismissing.(atmos_VPD_data)
        esat =
            Thermodynamics.saturation_vapor_pressure.(
                Ref(thermo_params),
                atmos_T_data[not_missing_mask] .+ 273.15,
                Ref(Thermodynamics.Liquid()),
            )
        e = @. esat - (atmos_VPD_data[not_missing_mask] .* 1000) # unit conversion of VPD to Pa
        RH = @. e / esat
        snow_frac = snow_precip_fraction.(atmos_T_data[not_missing_mask], RH)
        P_liq_data =
            -1 .* Float64.(atmos_precip_data[not_missing_mask]) .*
            (1 .- snow_frac)
        P_snow_data =
            -1 .* Float64.(atmos_precip_data[not_missing_mask]) .* snow_frac
        P_liq = TimeVaryingInput(
            seconds_since_start_date[not_missing_mask],
            P_liq_data ./ 1000 ./ DATA_DT,
        )
        P_snow = TimeVaryingInput(
            seconds_since_start_date[not_missing_mask],
            P_snow_data ./ 1000 ./ DATA_DT,
        )
    else
        not_missing_mask = .~ismissing.(atmos_precip_data)
        P_liq_data = -1 .* Float64.(atmos_precip_data[not_missing_mask])
        P_liq = TimeVaryingInput(
            seconds_since_start_date[not_missing_mask],
            P_liq_data ./ 1000 ./ DATA_DT,
        )
        P_snow = TimeVaryingInput(
            seconds_since_start_date[not_missing_mask],
            zeros(length(seconds_since_start_date[not_missing_mask])),
        )
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
        c_co2,
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
        SW_d,
        LW_d,
        start_date,
        θs = zenith_angle,
        earth_param_set = earth_param_set,
    )

    MODIS_LAI_path = ClimaLand.Artifacts.get_modis_lai_fluxnet_data(site_ID)
    MODIS_LAI, header = readdlm(MODIS_LAI_path, ',', header = true)
    dates = DateTime.(MODIS_LAI[:, 1])

    LAI_dt = Second(dates[2] - dates[1]).value
    LAI_seconds_since_start_date =
        [Second(date - start_date).value for date in dates]
    LAI = TimeVaryingInput(LAI_seconds_since_start_date, FT.(MODIS_LAI[:, 2]))
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

function set_fluxnet_ic!(
    Y,
    site_ID,
    start_date,
    hour_offset_from_UTC,
    model::LandModel,
)
    fluxnet_csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)
    driver_data = readdlm(fluxnet_csv_path, ',')
    columns = driver_data[1, :]
    labels = ("SWC_F_MDS_1", "TS_F_MDS_1", "TA_F")
    indices = [findfirst(columns .== label) for label in labels]
    index_map = Dict(zip(labels, indices))

    LOCAL_DATETIME = DateTime.(format.(driver_data[2:end, 1]), "yyyymmddHHMM")
    UTC_DATETIME = LOCAL_DATETIME .+ Dates.Hour(hour_offset_from_UTC)
    idx_start = findmin(abs.(UTC_DATETIME .- start_date))

    if index_map["SWC_F_MDS_1"] isa Nothing
        @. Y.soil.ϑ_l = model.soil.parameters.ν / 2
    else
        soil_data = driver_data[2:end, index_map["SWC_F_MDS_1"]]
        not_missing_mask = .~ismissing.(soil_data)
        idx_start = argmin(abs.(UTC_DATETIME[not_missing_mask] .- start_date))
        Y.soil.ϑ_l .= soil_data[not_missing_mask][idx_start] ./ 100
    end
    Y.soil.θ_i .= 0

    if index_map["TS_F_MDS_1"] isa Nothing
        atmos_T_data = driver_data[2:end, index_map["TA_F"]]
        not_missing_mask = .~ismissing.(atmos_T_data)
        idx_start = argmin(abs.(UTC_DATETIME[not_missing_mask] .- start_date))
        T_0 = atmos_T_data[not_missing_mask][idx_start] + 273.15
    else
        soil_T_data = driver_data[2:end, index_map["TS_F_MDS_1"]]
        not_missing_mask = .~ismissing.(soil_T_data)
        idx_start = argmin(abs.(UTC_DATETIME[not_missing_mask] .- start_date))
        T_0 = soil_T_data[not_missing_mask][idx_start] + 273.15
    end
    ρc_s =
        ClimaLand.Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            model.soil.parameters.ρc_ds,
            model.soil.parameters.earth_param_set,
        )
    Y.soil.ρe_int =
        ClimaLand.Soil.volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T_0,
            model.soil.parameters.earth_param_set,
        )
    Y.soilco2.C .= 0.000412 # set to atmospheric co2, mol co2 per mol air
    ψ_stem_0 = -1e5 / 9800 # pressure in the leaf divided by rho_liquid*gravitational acceleration [m]
    ψ_leaf_0 = -2e5 / 9800
    n_stem = model.canopy.hydraulics.n_stem
    n_leaf = model.canopy.hydraulics.n_leaf
    ψ_comps = n_stem > 0 ? [ψ_stem_0, ψ_leaf_0] : ψ_leaf_0
    S_l_ini =
        ClimaLand.Canopy.PlantHydraulics.inverse_water_retention_curve.(
            model.canopy.hydraulics.parameters.retention_model,
            ψ_comps,
            model.canopy.hydraulics.parameters.ν,
            model.canopy.hydraulics.parameters.S_s,
        )
    for i in 1:(n_stem + n_leaf)
        Y.canopy.hydraulics.ϑ_l.:($i) .=
            ClimaLand.Canopy.PlantHydraulics.augmented_liquid_fraction.(
                model.canopy.hydraulics.parameters.ν,
                S_l_ini[i],
            )
    end

    Y.canopy.energy.T .= T_0 # use the same initial condition as the soil

    Y.snow.S .= 0.0
    Y.snow.S_l .= 0.0
    Y.snow.U .= 0.0
end
end
