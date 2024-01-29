export setup_drivers

function setup_drivers(site_ID)
    context = ClimaComms.context()

    af = ArtifactFile(
        url = data_link,
        filename = "AMF_$(site_ID)_FLUXNET_FULLSET.csv",
    )
    dataset = ArtifactWrapper(
        "$ClimaLandSimulations_dir/src/fluxnet/$site_ID",
        "ameriflux_data",
        ArtifactFile[af],
    )
    dataset_path = get_data_folder(dataset)
    data = joinpath(dataset_path, "AMF_$(site_ID)_FLUXNET_FULLSET.csv")
    driver_data = readdlm(data, ',')

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
        :RECO,
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
        ("RECO_DT_VUT_REF", (x) -> x .* 1e-6, "mol/m^-2/s")
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
        error(
            "Driver data missing for columns: $([missing_drivers[i] * " " for 
          i in 1:length(missing_drivers)]...)",
        )
    end

    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    esat =
        Thermodynamics.saturation_vapor_pressure.(
            Ref(thermo_params),
            drivers.TA.values,
            Ref(Thermodynamics.Liquid()),
        )
    e = @. esat - drivers.VPD.values
    q = @. 0.622 * e ./ (drivers.PA.values - 0.378 * e)

    # Create interpolators for each atmospheric driver needed for PrescribedAtmosphere and
    # for PrescribedRadiation
    seconds = FT.(0:DATA_DT:((length(UTC_DATETIME) - 1) * DATA_DT))

    precip = TimeVaryingInput(seconds, -FT.(drivers.P.values[:]); context) # m/s
    atmos_q = TimeVaryingInput(seconds, FT.(q[:]); context)
    atmos_T = TimeVaryingInput(seconds, FT.(drivers.TA.values[:]); context)
    atmos_p = TimeVaryingInput(seconds, FT.(drivers.PA.values[:]); context)
    atmos_co2 = TimeVaryingInput(seconds, FT.(drivers.CO2.values[:]); context)
    atmos_u = TimeVaryingInput(seconds, FT.(drivers.WS.values[:]); context)
    LW_IN = TimeVaryingInput(seconds, FT.(drivers.LW_IN.values[:]); context)
    SW_IN = TimeVaryingInput(seconds, FT.(drivers.SW_IN.values[:]); context)
    snow_precip = TimeVaryingInput((t) -> FT(0))

    # Construct the drivers. For the reference time we will use the UTC time at the
    # start of the simulation
    atmos = ClimaLand.PrescribedAtmosphere(
        precip,
        snow_precip,
        atmos_T,
        atmos_u,
        atmos_q,
        atmos_p,
        UTC_DATETIME[1],
        atmos_h;
        c_co2 = atmos_co2,
    )

    function zenith_angle(
        t,
        ref_time;
        latitude = lat,
        longitude = long,
        insol_params::Insolation.Parameters.InsolationParameters{FT} = earth_param_set.insol_params,
    ) where {FT}
        # This should be time in UTC
        current_datetime = ref_time + Dates.Second(round(t))

        # Orbital Data uses Float64, so we need to convert to our sim FT
        d, δ, η_UTC =
            FT.(
                Insolation.helper_instantaneous_zenith_angle(
                    current_datetime,
                    ref_time,
                    insol_params,
                )
            )

        FT(
            Insolation.instantaneous_zenith_angle(
                d,
                δ,
                η_UTC,
                longitude,
                latitude,
            )[1],
        )
    end

    radiation = ClimaLand.PrescribedRadiativeFluxes(
        FT,
        SW_IN,
        LW_IN,
        UTC_DATETIME[1];
        θs = zenith_angle,
    )

    # Start and end dates of data in MODIS format
    modis_start_date = "A$(Dates.year(UTC_DATETIME[1]))$(lpad(Dates.dayofyear(UTC_DATETIME[1]), 3, "0"))"
    modis_end_date = "A$(Dates.year(UTC_DATETIME[end]))$(lpad(Dates.dayofyear(UTC_DATETIME[end]), 3, "0"))"

    MODIS_LAI = single_col_data_matrix(
        parse_response(
            check_response(
                send_get_subset(
                    "MCD15A2H",
                    modis_start_date,
                    modis_end_date,
                    site_ID,
                    band = "Lai_500m",
                ),
            ),
        ),
    )

    LAI_dt = Second(MODIS_LAI[2, 1] - MODIS_LAI[1, 1]).value
    LAI_seconds = FT.(0:LAI_dt:((length(MODIS_LAI[:, 1]) - 1) * LAI_dt))

    # LAI function for radiative transfer
    LAIfunction = TimeVaryingInput(LAI_seconds, FT.(MODIS_LAI[:, 2]); context)

    # Necessary inputs from LAI Data
    # Note that f_root_to_shoot, capacity, and h_leaf are site-specific parameters
    # defined in the parameters file for each site
    maxLAI = FT(maximum(MODIS_LAI[:, 2]))
    RAI = maxLAI * f_root_to_shoot
    plant_ν = capacity / (maxLAI * h_leaf) / FT(1000)

    return (
        LOCAL_DATETIME,
        atmos_co2,
        DATA_DT,
        drivers,
        atmos,
        radiation,
        LAIfunction,
        maxLAI,
        RAI,
        plant_ν,
    )
end
