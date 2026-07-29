using NCDatasets
import ClimaParams as CP


"""
     prescribed_forcing_fluxnet(site_ID,
                                lat,
                                long,
                                hour_offset_from_UTC,
                                start_date,
                                toml_dict::CP.ParamDict,
                                FT;
                                gustiness=1,
                                split_precip= true,
)
A helper function which constructs the `PrescribedAtmosphere` and `PrescribedRadiativeFluxes`
from a file path pointing to the Fluxnet data in a csv file, the start date, latitude, longitude,
the hour offset of the site from UTC (local_time + offset = time in UTC),
and the `toml_dict`.

This requires (1) reading in the data, (2) removing missing values,
 (3) converting units, (4) computing the specific humidity and percent of
precipitation that is in snow (if split_precip ==true), and (5)
making the TimeVaryingInput objects.

Note that the TimeVaryingInput objects can be used to interpolate in time,
which is why we drop missing data. For the timestamp of the observation, we use
the halfway point between the timestamp start and timestamp end:
see https://fluxnet.org/data/fluxnet2015-dataset/fullset-data-product/ for details.

This assumes that the first row of the CSV file is the list of of column names,
and that these names are:
"TA_F" (air temperature in C)
"VPD_F" (air vapor pressure deficit in hPa)
"PA_F" (air pressure in kPa)
"P_F" (accumulated precipitation in mm)
"WS_F" (wind speed in m/s)
"LW_IN_F" (downwelling LW radiation in W/m^2)
"SW_IN_F" (downwelling SW radiation in W/m^2)
"CO2_F_MDS" (CO2 concentration in μmol/mol)
"""
function FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    hour_offset_from_UTC,
    atmos_h,
    start_date, # in UTC
    toml_dict::CP.ParamDict,
    FT;
    split_precip = true,
    gustiness = 1,
)
    earth_param_set = LP.LandParameters(toml_dict)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)

    # Read data and process timestamps
    (data, columns) = read_fluxnet_data(site_ID)

    # Determine which column index corresponds to which varname
    varnames = (
        "TIMESTAMP_START",
        "TIMESTAMP_END",
        "TA_F",
        "VPD_F",
        "PA_F",
        "P_F",
        "WS_F",
        "LW_IN_F",
        "SW_IN_F",
        "CO2_F_MDS",
    )
    column_name_map =
        get_column_name_map(varnames, columns; error_on_missing = false)

    # Convert the local timestamps to UTC, using the midpoint of each averaging
    # period as the corresponding datetime
    UTC_datetimes_start = get_UTC_datetimes(
        hour_offset_from_UTC,
        data,
        column_name_map;
        timestamp_name = "TIMESTAMP_START",
    )
    UTC_datetimes_end = get_UTC_datetimes(
        hour_offset_from_UTC,
        data,
        column_name_map;
        timestamp_name = "TIMESTAMP_END",
    )
    UTC_datetimes =
        UTC_datetimes_start .+ (UTC_datetimes_end .- UTC_datetimes_start) ./ 2

    # The TimeVaryingInput interface for columns expects the time in seconds
    # from the start date of the simulation
    seconds_since_start_date =
        [Second(UTC - start_date).value for UTC in UTC_datetimes]

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
    c_co2 = time_varying_input_from_data(
        data,
        "CO2_F_MDS",
        column_name_map,
        seconds_since_start_date;
        preprocess_func = (x) -> x * 1e-6, # convert from μmol/mol to mol/mol
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
    data_dt = Second(UTC_datetimes[2] - UTC_datetimes[1]).value # seconds
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
        toml_dict;
        c_co2,
    )

    cos_zenith_angle =
        (t, s) -> default_cos_zenith_angle(
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
        cosθs = cos_zenith_angle,
        toml_dict = toml_dict,
    )
    return (; atmos, radiation)
end

"""
    get_maxLAI_at_site(date, lat, long;
                       ncd_path = ClimaLand.Artifacts.modis_lai_single_year_path(;year = Dates.year(date)))

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
        year = Dates.year(date),
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
    (data, columns) = read_fluxnet_data(site_ID)

    varnames = ("TIMESTAMP_START",)
    column_name_map =
        get_column_name_map(varnames, columns; error_on_missing = true)

    # Since we only need to know the difference in time between observations,
    # we can use a dummy hour offset from UTC, giving us the local times.
    dummy_hour_offset = 0
    local_datetimes = get_UTC_datetimes(
        dummy_hour_offset,
        data,
        column_name_map;
        timestamp_name = varnames[1],
    )

    # Get the difference in time between each timestamp pair, in seconds
    Δts = [
        Second(Δdate).value for
        Δdate in (local_datetimes[2:end] .- local_datetimes[1:(end - 1)])
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

Please note that we use date halfway between the TIMESTAMP_START
and TIMESTAMP_END date of the Fluxnet averaging window as the timestamp of
the data for forcing.
"""
function FluxnetSimulations.get_data_dates(
    site_ID,
    hour_offset_from_UTC;
    duration::Union{Nothing, Period} = nothing,
    start_offset::Period = Second(0),
)
    (data, columns) = read_fluxnet_data(site_ID)

    # If any of these are missing, error, because we need all of them
    varnames = ("TIMESTAMP_START", "TIMESTAMP_END")
    column_name_map =
        get_column_name_map(varnames, columns; error_on_missing = true)

    # Convert the local timestamps to UTC
    UTC_datetimes_start = get_UTC_datetimes(
        hour_offset_from_UTC,
        data,
        column_name_map;
        timestamp_name = "TIMESTAMP_START",
    )
    UTC_datetimes_end = get_UTC_datetimes(
        hour_offset_from_UTC,
        data,
        column_name_map;
        timestamp_name = "TIMESTAMP_END",
    )
    # Get the midpoint of the averaging period
    UTC_datetimes =
        UTC_datetimes_start .+ (UTC_datetimes_end .- UTC_datetimes_start) ./ 2
    earliest_date, latest_date = extrema(UTC_datetimes)

    # Compute the start and stop dates, applying the start_offset and duration if provided
    @assert Dates.value(start_offset) >= 0 "start_offset must be non-negative, got $start_offset"
    start_date = earliest_date + start_offset
    stop_date = isnothing(duration) ? latest_date : start_date + duration

    # If duration is provided, check that it is non-negative and that the stop date
    # is not greater than the latest date in the data
    if !isnothing(duration)
        @assert Dates.value(duration) >= 0 "If duration is provided, it must be non-negative."
        @assert stop_date <= latest_date "The sum of the requested duration $duration and start_offset $start_offset \
            is greater than the available duration $(latest_date - earliest_date) of data."
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
    _T_freeze = Thermodynamics.Parameters.T_freeze(thermo_params)
    # Convert units
    T = T + _T_freeze  # C to K
    VPD = VPD * 100.0 # hPa to Pa
    P = P * 1000.0 # kPa to Pa
    # Compute q
    esat = Thermodynamics.saturation_vapor_pressure(
        thermo_params,
        T,
        Thermodynamics.Liquid(),
    )
    e = esat - VPD
    ρ = Thermodynamics.air_density(thermo_params, T, P) # this assumes the gas constant for dry air, which induces a small error
    q = Thermodynamics.q_vap_from_p_vap(thermo_params, T, ρ, e)
    return q
end

"""
     get_comparison_data(
        site_ID,
        hour_offset_from_UTC;
        val = -9999,
)

Gets and returns the a NamedTuple with the comparison data for the Fluxnet `site_ID`,
given its hour offset from UTC, and the value used to indicate missing data (`val`).

Note that the time associated with the comparison data is the start of the averaging
period used in Fluxnet, which aligns with convention used for saving ClimaLand diagnostics.
"""
function FluxnetSimulations.get_comparison_data(
    site_ID::String,
    hour_offset_from_UTC::Int;
    val = -9999,
)
    # Read data and process timestamps
    (data, columns) = read_fluxnet_data(site_ID)

    # Determine which column index corresponds to which varname
    varnames = (
        "TIMESTAMP_START",
        "TIMESTAMP_END",
        "GPP_DT_VUT_REF",
        "LE_CORR",
        "H_CORR",
        "SW_OUT",
        "LW_OUT",
        "SWC_F_MDS_1",
        "SWC_F_MDS_2",
        "SWC_F_MDS_3",
        "SWC_F_MDS_4",
        "SWC_F_MDS_5",
        "TS_F_MDS_1",
        "TS_F_MDS_2",
        "TS_F_MDS_3",
        "TS_F_MDS_4",
        "TS_F_MDS_5",
        "P_F",
    )

    # For comparison data, we don't need to error on missing variables
    column_name_map =
        get_column_name_map(varnames, columns; error_on_missing = false)
    UTC_datetimes = get_UTC_datetimes(
        hour_offset_from_UTC,
        data,
        column_name_map,
        timestamp_name = "TIMESTAMP_START",
    )
    data_dt = Second(UTC_datetimes[2] - UTC_datetimes[1]).value # seconds

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

    swc2 = FluxnetSimulations.get_comparison_data(
        data,
        "SWC_F_MDS_2",
        column_name_map,
        "swc2";
        preprocess_func = (x) -> x / 100, # converts a percent to an absolute number
        val,
    )
    swc3 = FluxnetSimulations.get_comparison_data(
        data,
        "SWC_F_MDS_3",
        column_name_map,
        "swc3";
        preprocess_func = (x) -> x / 100, # converts a percent to an absolute number
        val,
    )
    swc4 = FluxnetSimulations.get_comparison_data(
        data,
        "SWC_F_MDS_4",
        column_name_map,
        "swc4";
        preprocess_func = (x) -> x / 100, # converts a percent to an absolute number
        val,
    )
    swc5 = FluxnetSimulations.get_comparison_data(
        data,
        "SWC_F_MDS_5",
        column_name_map,
        "swc5";
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

    tsoil2 = FluxnetSimulations.get_comparison_data(
        data,
        "TS_F_MDS_2",
        column_name_map,
        "tsoil2";
        preprocess_func = (x) -> x + 273.15, # converts degrees C to K
        val,
    )
    tsoil3 = FluxnetSimulations.get_comparison_data(
        data,
        "TS_F_MDS_3",
        column_name_map,
        "tsoil3";
        preprocess_func = (x) -> x + 273.15, # converts degrees C to K
        val,
    )
    tsoil4 = FluxnetSimulations.get_comparison_data(
        data,
        "TS_F_MDS_4",
        column_name_map,
        "tsoil4";
        preprocess_func = (x) -> x + 273.15, # converts degrees C to K
        val,
    )
    tsoil5 = FluxnetSimulations.get_comparison_data(
        data,
        "TS_F_MDS_5",
        column_name_map,
        "tsoil5";
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
        (; UTC_datetime = UTC_datetimes),
        gpp,
        lhf,
        shf,
        swu,
        lwu,
        swc,
        swc2,
        swc3,
        swc4,
        swc5,
        tsoil,
        tsoil2,
        tsoil3,
        tsoil4,
        tsoil5,
        precip,
    )
end

# Returns a `tvi(varnames; compose, preprocess)` closure which wraps each variable in a
# `ClimaUtilities.DataHandling.MultiColumnDataHandler`: one source file per column, read
# onto the matching column of `surface_space`. Shared by `prescribed_forcing_netcdf` and
# `prescribed_LAI_netcdf`.
#
# The handler matches sources to columns by each file's scalar `longitude`/`latitude`, so
# `met_nc_paths` need not be in column order. Because of that, the local-standard-time to
# UTC shift is carried by each column's `DataSource` `time_transform` rather than applied
# by column index: the offset travels with its file through the matching.
function _netcdf_tvi_builder(
    met_nc_paths,
    surface_space,
    start_date;
    hour_offset_from_UTC = 0,
    time_interpolation_method = LinearInterpolation(),
)
    ncolumns = length(met_nc_paths)
    offsets =
        hour_offset_from_UTC isa AbstractVector ? hour_offset_from_UTC :
        fill(hour_offset_from_UTC, ncolumns)
    length(offsets) == ncolumns || error(
        "hour_offset_from_UTC must be a scalar or a length-$ncolumns vector",
    )
    # A comprehension binds a fresh `offset` per column, so each closure captures its own.
    time_transforms =
        [date -> date - hour_offset_to_period(offset) for offset in offsets]

    function tvi(varnames; compose = identity, preprocess = identity)
        varnames isa AbstractString && (varnames = [varnames])
        # The outer vector selects the variable, and its order sets the argument order of
        # `compose`; each inner vector holds that variable's columns in `met_nc_paths` order.
        dataset_sources = [
            [
                DataSource(
                    met_nc_paths[column],
                    varname;
                    time_transform = time_transforms[column],
                ) for column in 1:ncolumns
            ] for varname in varnames
        ]
        data_handler = DataHandling.MultiColumnDataHandler(
            dataset_sources,
            surface_space;
            start_date,
            compose_function = compose,
            file_reader_kwargs = (; preprocess_func = preprocess),
        )
        return TimeVaryingInput(
            data_handler;
            method = time_interpolation_method,
        )
    end
    return tvi
end

"""
    prescribed_forcing_netcdf(met_nc_paths,
                              surface_space,
                              atmos_h,
                              start_date,
                              toml_dict::CP.ParamDict,
                              FT;
                              hour_offset_from_UTC = 0,
                              split_precip = true,
                              gustiness = 1,
                              time_interpolation_method = LinearInterpolation())

Constructs the `PrescribedAtmosphere` and `PrescribedRadiativeFluxes` from site-level
FLUXNET-style NetCDF forcing files, read through
`ClimaUtilities.DataHandling.MultiColumnDataHandler`.

Unlike `prescribed_forcing_fluxnet`, which reads a CSV into memory and builds
time series that carry no spatial information, this reads the files lazily onto
`surface_space` and so returns forcing that varies in space as well as time. It is
therefore restricted to domains whose surface is a `ClimaCore.Spaces.PointCloudSpace`,
i.e. a `ClimaLand.Domains.ColumnEnsemble` — a single-column `ColumnEnsemble` is the
one-site case, and a plain `ClimaLand.Domains.Column` is not supported
because its surface is a `PointSpace`.

`met_nc_paths` is one path (one column) or a vector of paths, one per column. The
handler assigns each file to a column by the file's scalar `longitude`/`latitude`
variables, matched within a tolerance of 1e-3 degrees, so the paths need not be in
column order — but each column of `surface_space` must have a file whose coordinates
agree with it to that tolerance.

Expected variables in each NetCDF file:
- `Tair`   : air temperature (K)
- `Wind`   : wind speed (m/s)
- `Qair`   : specific humidity (kg/kg)
- `Psurf`  : surface pressure (Pa)
- `SWdown` : downwelling shortwave radiation (W/m²)
- `LWdown` : downwelling longwave radiation (W/m²)
- `CO2air` : atmospheric CO₂ concentration (ppm); converted to mol/mol
- `Precip` : precipitation rate (kg/m²/s); converted to a volume flux and sign-flipped
- `VPD`    : vapor pressure deficit (hPa); used for the snow fraction

The NetCDF `time` is assumed to be in FLUXNET local standard time and is converted to
UTC as `UTC = local - hour_offset_from_UTC`, matching `get_UTC_datetimes` and the
`time_offset` reported by `get_location` (so CET is `+1` and CST is `-6`).
`hour_offset_from_UTC` is a scalar applied to every column, or a vector giving each
column's offset in `met_nc_paths` order. `start_date` must be in UTC.

Unlike the CSV path, entries flagged as missing cannot be dropped, since the handler
interpolates between the dates present in the file. Gap-filled forcing is expected;
pass a `preprocess_func` through `file_reader_kwargs` upstream if a file needs masking.
"""
function FluxnetSimulations.prescribed_forcing_netcdf(
    met_nc_paths,
    surface_space,
    atmos_h,
    start_date, # in UTC
    toml_dict::CP.ParamDict,
    FT;
    hour_offset_from_UTC = 0,
    split_precip = true,
    gustiness = 1,
    time_interpolation_method = LinearInterpolation(),
)
    met_nc_paths isa AbstractString && (met_nc_paths = [met_nc_paths])
    earth_param_set = LP.LandParameters(toml_dict)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)

    tvi = _netcdf_tvi_builder(
        met_nc_paths,
        surface_space,
        start_date;
        hour_offset_from_UTC,
        time_interpolation_method,
    )

    atmos_T = tvi("Tair")
    atmos_u = tvi("Wind")
    atmos_q = tvi("Qair")
    atmos_P = tvi("Psurf")
    SW_d = tvi("SWdown")
    LW_d = tvi("LWdown")
    c_co2 = tvi("CO2air"; preprocess = x -> x * 1e-6) # ppm -> mol/mol

    if split_precip
        # `Precip` is a mass flux in kg/m^2/s; dividing by ρ_water (1000 kg/m^3) gives the
        # volume flux in m/s that `PrescribedAtmosphere` expects, negated because it is
        # directed into the land. `snow_precip_fraction` takes T in °C and VPD in hPa, and
        # is broadcast explicitly rather than with `@.` because it takes a keyword argument.
        atmos_P_liq = tvi(
            ["Precip", "Tair", "VPD"];
            compose = (P, T, VPD) ->
                -(P ./ 1000) .*
                (1 .- snow_precip_fraction.(T .- 273.15, VPD; thermo_params)),
        )
        atmos_P_snow = tvi(
            ["Precip", "Tair", "VPD"];
            compose = (P, T, VPD) ->
                -(P ./ 1000) .*
                snow_precip_fraction.(T .- 273.15, VPD; thermo_params),
        )
    else
        atmos_P_liq = tvi("Precip"; preprocess = x -> -x / 1000)
        atmos_P_snow = tvi("Precip"; preprocess = x -> zero(x))
    end

    atmos = ClimaLand.PrescribedAtmosphere(
        atmos_P_liq,
        atmos_P_snow,
        atmos_T,
        atmos_u,
        atmos_q,
        atmos_P,
        start_date,
        atmos_h,
        toml_dict;
        c_co2,
        gustiness = FT(gustiness),
    )

    # The zenith angle is per column, so latitude and longitude come from the surface
    # space rather than from scalar site coordinates.
    cos_zenith_angle =
        (t, s) -> default_cos_zenith_angle(
            t,
            s;
            insol_params = earth_param_set.insol_params,
            longitude = ClimaLand.Domains.get_long(surface_space),
            latitude = ClimaLand.Domains.get_lat(surface_space),
        )
    radiation = ClimaLand.PrescribedRadiativeFluxes(
        FT,
        SW_d,
        LW_d,
        start_date;
        cosθs = cos_zenith_angle,
        toml_dict,
    )
    return (; atmos, radiation)
end

"""
    prescribed_LAI_netcdf(met_nc_paths,
                          surface_space,
                          start_date;
                          varname = "LAI_alternative",
                          hour_offset_from_UTC = 0,
                          time_interpolation_method = LinearInterpolation())

Constructs a surface `TimeVaryingInput` of leaf area index from the `varname` variable of
the site-level NetCDF files in `met_nc_paths`, read onto `surface_space` through the same
`MultiColumnDataHandler` machinery as `prescribed_forcing_netcdf`.

Not every FLUXNET-style met file carries an LAI variable; `varname` selects which one to
read, and the file must contain it.
"""
function FluxnetSimulations.prescribed_LAI_netcdf(
    met_nc_paths,
    surface_space,
    start_date;
    varname = "LAI_alternative",
    hour_offset_from_UTC = 0,
    time_interpolation_method = LinearInterpolation(),
)
    met_nc_paths isa AbstractString && (met_nc_paths = [met_nc_paths])
    tvi = _netcdf_tvi_builder(
        met_nc_paths,
        surface_space,
        start_date;
        hour_offset_from_UTC,
        time_interpolation_method,
    )
    return tvi(varname)
end

"""
    get_location_netcdf(FT, met_nc_path)

Returns the `(; long, lat)` of a site-level FLUXNET-style NetCDF forcing file, read from
its `longitude` and `latitude` variables.

`prescribed_forcing_netcdf` matches each file to a column of the surface space by
these coordinates within a tolerance of 1e-3 degrees, so use this to build the
`ClimaLand.Domains.ColumnEnsemble` when a site's coordinates are not otherwise
known to the same precision as its file.
"""
function FluxnetSimulations.get_location_netcdf(FT, met_nc_path)
    NCDataset(met_nc_path, "r") do ds
        return (;
            long = FT(only(ds["longitude"])),
            lat = FT(only(ds["latitude"])),
        )
    end
end
