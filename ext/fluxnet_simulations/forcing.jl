using NCDatasets

"""
     prescribed_forcing_fluxnet(site_ID,
                                lat,
                                long,
                                hour_offset_from_UTC,
                                start_date,
                                earth_param_set,
                                FT;
                                gustiness=1,
                                split_precip= true,
                                c_co2 = TimeVaryingInput((t) -> 4.2e-4),
)
A helper function which constructs the `PrescribedAtmosphere` and `PrescribedRadiativeFluxes`
from a file path pointing to the Fluxnet data in a csv file, the start date, latitude, longitude,
the hour offset of the site from UTC (local_time + offset = time in UTC),
and the earth_param_set.

This requires (1) reading in the data, (2) removing missing values,
 (3) converting units, (4) computing the specific humidity and percent of
precipitation that is in snow (if split_precip ==true), and (5)
making the TimeVaryingInput objects.

Note that the TimeVaryingInput objects can be used to interpolate in time,
which is why we drop missing data.

This assumes that the first row of the CSV file is the list of of column names,
and that these names are:
"TA_F" (air temperature in C)
"VPD_F" (air vapor pressure deficit in hPa)
"PA_F" (air pressure in kPa)
"P_F" (accumulated precipitation in mm)
"WS_F" (wind speed in m/s)
"LW_IN_F" (downwelling LW radiation in W/m^2)
"SW_IN_F" (downwelling SW radiation in W/m^2)
"""
function FluxnetSimulations.prescribed_forcing_fluxnet(
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
)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)

    fluxnet_csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)
    (data, columns) = readdlm(fluxnet_csv_path, ','; header = true)

    # Determine which column index corresponds to which varname
    varnames = ("TA_F", "VPD_F", "PA_F", "P_F", "WS_F", "LW_IN_F", "SW_IN_F")
    column_name_map = Dict(
        varname => findfirst(columns[:] .== varname) for varname in varnames
    )
    # If any of these are missing, error, because we need all of them
    # to run a simulation
    nothing_id = findall(collect(values(column_name_map)) .== nothing)
    if !isempty(nothing_id)
        @error("$(labels[nothing_id]) is missing in the data, but required.")
    end

    # Convert the local timestamp to UTC
    # Since it was read in as Float64 type, convert to a string before
    # converting to a DateTime
    local_datetime = DateTime.(string.(Int.(data[:, 1])), "yyyymmddHHMM")
    UTC_datetime = local_datetime .+ Dates.Hour(hour_offset_from_UTC)

    # The TimeVaryingInput interface for columns expects the time in seconds
    # from the start date of the simulation
    seconds_since_start_date =
        [Second(UTC - start_date).value for UTC in UTC_datetime]

    # Create the TVI objects
    atmos_T = time_varying_input_from_data(
        data,
        "TA_F",
        column_name_map,
        seconds_since_start_date;
        preprocess_func = (x) -> x + 273.15,
    )
    atmos_P = time_varying_input_from_data(
        data,
        "PA_F",
        column_name_map,
        seconds_since_start_date;
        preprocess_func = (x) -> x * 1000,
    )
    atmos_u = time_varying_input_from_data(
        data,
        "WS_F",
        column_name_map,
        seconds_since_start_date,
    )
    LW_d = time_varying_input_from_data(
        data,
        "LW_IN_F",
        column_name_map,
        seconds_since_start_date,
    )
    SW_d = time_varying_input_from_data(
        data,
        "SW_IN_F",
        column_name_map,
        seconds_since_start_date,
    )

    # Specific humidity is computed from P, VPD, and T using `compute_q`
    # This is computed from multiple columns of data; the argument order of q_closure, compute_q, must match the
    # column order in the second argument of the time_varying_input function.
    q_closure(P, VPD, T) = compute_q(P, VPD, T; thermo_params)
    atmos_q = time_varying_input_from_data(
        data,
        ["PA_F", "VPD_F", "TA_F"],
        column_name_map,
        seconds_since_start_date;
        preprocess_func = q_closure,
    )

    # Next is precipitation, which is reported as an accumulation over
    # the time between observations in the data
    # Convert to a flux using data_dt. Also, change the sign, because in ClimaLand
    # precipitation is a downward flux (negative), and convert from accumulated mm to m/s
    data_dt = Second(local_datetime[2] - local_datetime[1]).value # seconds
    if split_precip
        compute_rain(T, VPD, precip; thermo_params = thermo_params) =
            -1 * precip / 1000 / data_dt *
            (1 - snow_precip_fraction(T, VPD; thermo_params))
        compute_snow(T, VPD, precip; thermo_params = thermo_params) =
            -1 * precip / 1000 / data_dt *
            snow_precip_fraction(T, VPD; thermo_params)
        atmos_P_liq = time_varying_input_from_data(
            data,
            ["TA_F", "VPD_F", "P_F"],
            column_name_map,
            seconds_since_start_date;
            preprocess_func = compute_rain,
        )
        atmos_P_snow = time_varying_input_from_data(
            data,
            ["TA_F", "VPD_F", "P_F"],
            column_name_map,
            seconds_since_start_date;
            preprocess_func = compute_snow,
        )
    else
        atmos_P_liq = time_varying_input_from_data(
            data,
            "P_F",
            column_name_map,
            seconds_since_start_date;
            preprocess_func = (x) -> -x / 1000 / data_dt,
        )
        atmos_P_snow = time_varying_input_from_data(
            data,
            "P_F",
            column_name_map,
            seconds_since_start_date;
            preprocess_func = (x) -> zero(x),
        ) # no snow
    end

    # Construct the drivers. For the start date we use the UTC time at the
    # start of the simulation
    atmos = ClimaLand.PrescribedAtmosphere(
        atmos_P_liq,
        atmos_P_snow,
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
    return (; atmos, radiation)
end

"""
    get_maxLAI_at_site(date, lat, long;
                       ncd_path = ClimaLand.Artifacts.modis_lai_single_year_path(;year = Year(date)))

A helper function to get the maximum LAI at a site from the MODIS LAI data for the 
year corresponding to `date`. This is used in some simulations to 
set the root area index.

By default, this uses MODIS data for the year desired; the default file 
at `ncd_path` is expected to contain latitude and longitude in variables
"lat" and "lon", and LAI in the variable "lai". The function finds the maximum LAI
at the closest latitude and longitude to the given `lat` and `long` values over all dates
in the dataset.

The user can override the default by passing a different ncd_path, in which case `date` is not used.
"""
function FluxnetSimulations.get_maxLAI_at_site(
    date,
    lat,
    long;
    ncd_path = ClimaLand.Artifacts.modis_lai_single_year_path(;
        year = Year(date),
    ),
)
    NCDataset(ncd_path) do ds
        # Find the indices of the closest latitude and longitude in the dataset
        # to the given lat and long
        j = findmin(x -> abs(x - lat), ds["lat"][:])[2]
        i = findmin(x -> abs(x - long), ds["lon"][:])[2]

        # Get the maximum LAI value for the given site and time period
        max_lai = maximum(ds["lai"][i, j, :])
        return max_lai
    end
end

"""
    get_data_dt(site_ID)

A helper function to get the difference in time between observations;
this is used in making some plots.
"""
function FluxnetSimulations.get_data_dt(site_ID)
    fluxnet_csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)
    (data, columns) = readdlm(fluxnet_csv_path, ','; header = true)
    local_datetime = DateTime.(string.(Int.(data[:, 1])), "yyyymmddHHMM")
    # Convert to seconds
    Δts = [
        Second(Δdate).value for
        Δdate in (local_datetime[2:end] .- local_datetime[1:(end - 1)])
    ]
    # Make sure all the Δts are the same
    @assert all(Δts .== Δts[1])
    return Float64(Δts[1])
end

"""
    get_data_dates(
        site_ID,
        hour_offset_from_UTC;
        duration::Union{Nothing, Period} = nothing,
        start_offset::Period = Second(0),
    )

A helper function to get the first and last dates, in UTC, for which we have
Fluxnet data at `site_ID`, given the offset in hours of local time
from UTC. If `duration` is provided, it is used to determine the end date,
otherwise the end date is the last date in the data. The `start_offset` is
added to the start date, and must be non-negative.
"""
function FluxnetSimulations.get_data_dates(
    site_ID,
    hour_offset_from_UTC;
    duration::Union{Nothing, Period} = nothing,
    start_offset::Period = Second(0),
)
    fluxnet_csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)
    (data, columns) = readdlm(fluxnet_csv_path, ','; header = true)
    local_datetime = DateTime.(string.(Int.(data[:, 1])), "yyyymmddHHMM")
    UTC_datetime = local_datetime .+ Dates.Hour(hour_offset_from_UTC)
    earliest_date, latest_date = extrema(UTC_datetime)
    Dates.value(start_offset) < 0 && error("start_offset must be non-negative")
    if !isnothing(duration) && Dates.value(duration) < 0
        error("If duration is not provided, it must be non-negative.")
    end
    duration_available_ms = latest_date - earliest_date
    start_date = earliest_date + start_offset
    stop_date = isnothing(duration) ? latest_date : start_date + duration
    if !isnothing(duration) && (stop_date > latest_date)
        error(
            "The sum of the requested duration of $duration and start_offset of $start_offset \
            is greater than the available $duration_available_ms of data.",
        )
    end
    return (start_date, stop_date)
end

"""
    snow_precip_fraction(T_in_C, VPD_in_hPa; thermo_params)

Estimate the fraction of precipitation that is in snow form,
given the air temperature at the surface in C and the vapor pressure
deficit in hPa.

See Jennings, K.S., Winchell, T.S., Livneh, B. et al.
Spatial variation of the rain–snow temperature threshold across the
Northern Hemisphere. Nat Commun 9, 1148 (2018).
https://doi.org/10.1038/s41467-018-03629-7
"""
function snow_precip_fraction(T_in_C, VPD_in_hPa; thermo_params)
    T_in_K = T_in_C + 273.15
    VPD_in_Pa = VPD_in_hPa * 100
    esat = Thermodynamics.saturation_vapor_pressure(
        thermo_params,
        T_in_K,
        Thermodynamics.Liquid(),
    )
    e = esat - VPD_in_Pa
    RH = e / esat
    α = -10.04
    β = 1.41
    γ = 0.09
    snow_frac = (1.0 / (1.0 + exp(α + β * T_in_C + γ * RH)))
    return snow_frac
end

"""
    compute_q(P, VPD, T; thermo_params)

Returns the specific humidity (kg/kg) given the
air pressure `P`, in kPa, vapor pressure deficit
`VPD`, in hPa, and temperature `T` in degrees C;
these units are what Fluxnet sites report the data
in.
"""
function compute_q(P, VPD, T; thermo_params)
    # Convert units
    T = T + 273.15  # C to K
    VPD = VPD * 100.0 # hPa to Pa
    P = P * 1000.0 # kPa to Pa
    # Compute q
    esat = Thermodynamics.saturation_vapor_pressure(
        thermo_params,
        T,
        Thermodynamics.Liquid(),
    )
    e = esat - VPD
    q = 0.622 * e / (P - 0.378 * e)
    return q
end

"""
     get_comparison_data(
        site_ID,
        hour_offset_from_UTC;
        val = -9999
)

Gets and returns the a NamedTuple with the comparison
data for the Fluxnet `site_ID`, given its hour offset from
UTC, and the value used to indicate missing data (`val`)
"""
function FluxnetSimulations.get_comparison_data(
    site_ID::String,
    hour_offset_from_UTC::Int;
    val = -9999,
)
    fluxnet_csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)
    (data, columns) = readdlm(fluxnet_csv_path, ','; header = true)

    # Determine which column index corresponds to which varname
    varnames = (
        "GPP_DT_VUT_REF",
        "LE_CORR",
        "H_CORR",
        "SW_OUT",
        "LW_OUT",
        "SWC_F_MDS_1",
        "TS_F_MDS_1",
        "P_F",
    )
    column_name_map = Dict(
        varname => findfirst(columns[:] .== varname) for varname in varnames
    )

    # Convert the local timestamp to UTC
    local_datetime = DateTime.(string.(Int.(data[:, 1])), "yyyymmddHHMM")
    UTC_datetime = local_datetime .+ Dates.Hour(hour_offset_from_UTC)
    data_dt = Second(local_datetime[2] - local_datetime[1]).value # seconds

    gpp = FluxnetSimulations.get_comparison_data(
        data,
        "GPP_DT_VUT_REF",
        column_name_map,
        "gpp";
        preprocess_func = (x) -> x * 1e-6, # converts from μmol/m^2/s to mol/m^2/s
        val,
    )
    lhf = FluxnetSimulations.get_comparison_data(
        data,
        "LE_CORR",
        column_name_map,
        "lhf";
        val,
    )
    shf = FluxnetSimulations.get_comparison_data(
        data,
        "H_CORR",
        column_name_map,
        "shf";
        val,
    )
    swu = FluxnetSimulations.get_comparison_data(
        data,
        "SW_OUT",
        column_name_map,
        "swu";
        val,
    )
    lwu = FluxnetSimulations.get_comparison_data(
        data,
        "LW_OUT",
        column_name_map,
        "lwu";
        val,
    )
    swc = FluxnetSimulations.get_comparison_data(
        data,
        "SWC_F_MDS_1",
        column_name_map,
        "swc";
        preprocess_func = (x) -> x / 100, # converts a percent to an absolute number
        val,
    )
    tsoil = FluxnetSimulations.get_comparison_data(
        data,
        "TS_F_MDS_1",
        column_name_map,
        "tsoil";
        preprocess_func = (x) -> x + 273.15, # converts degrees C to K
        val,
    )
    precip = FluxnetSimulations.get_comparison_data(
        data,
        "P_F",
        column_name_map,
        "precip";
        preprocess_func = (x) -> -x / 1000 / data_dt, # converts accumulated mm to m/s
        val,
    )
    return merge(
        (; UTC_datetime = UTC_datetime),
        gpp,
        lhf,
        shf,
        swu,
        lwu,
        swc,
        tsoil,
        precip,
    )
end
