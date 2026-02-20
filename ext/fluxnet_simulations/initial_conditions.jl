"""
    FluxnetSimulations.make_set_fluxnet_initial_conditions(
        site_ID,
        start_date,
        hour_offset_from_UTC,
        model,
    )

Creates and returns a function `set_ic!(Y,p,t,model)` which
updates `Y` in place with an estimated set of initial conditions
based on the fluxnet observations at `site_ID` at the `start_date` in UTC,
and the type of the `model`.
In order to convert between local time and UTC, the hour offset from
UTC is required.
"""
function FluxnetSimulations.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    hour_offset_from_UTC,
    model,
)
    set_ic!(Y, p, t, model) =
        set_fluxnet_ic!(Y, site_ID, start_date, hour_offset_from_UTC, model)
    return set_ic!
end


"""
     set_fluxnet_ic!(
        Y,
        site_ID,
        start_date,
        hour_offset_from_UTC,
        model::ClimaLand.AbstractLandModel,
    )

Sets the initial conditions of `Y` using observations from the site `site_ID`, if available,
using the observations closest to the start_date (in UTC). Since the data from Fluxnet sites
is provided in local time, we require the offset from UTC in hours `hour_offset_from_UTC`.
The `model` indicates which how to update it `Y` from these observations,
via different methods of `set_fluxnet_ic!`.
"""
function set_fluxnet_ic!(
    Y,
    site_ID,
    start_date,
    hour_offset_from_UTC,
    model::ClimaLand.AbstractLandModel,
)
    fluxnet_csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)

    # Read the data and get the column name map
    (data, columns) = readdlm(fluxnet_csv_path, ','; header = true)

    # Get the UTC datetimes of the data
    varnames = ("TIMESTAMP_START",)
    column_name_map =
        get_column_name_map(varnames, columns; error_on_missing = true)
    UTC_datetimes = get_UTC_datetimes(
        hour_offset_from_UTC,
        data,
        column_name_map;
        timestamp_name = "TIMESTAMP_START",
    )
    Δ_date = UTC_datetimes .- start_date
    for component in ClimaLand.land_components(model)
        set_fluxnet_ic!(Y, data, columns, Δ_date, getproperty(model, component))
    end
end

"""
     set_fluxnet_ic!(Y, data, columns, Δ_date, model::ClimaLand.Soil.EnergyHydrology)

Sets the values of Y.soil in place with:
- \vartheta_l: observed value of SWC at the surface at the observation date closest to the start date, unless this is larger than 90% of porosity.
- θ_i: no ice (θ_i = 0)
- \rho e_int: an internal energy computed using the above θ_l, θ_i, and the temperature of the soil
  in the first layer, at the observation date closest to the start date. If the soil
  temperature is not available, the air temperature is used.

Here, `Y` is the prognostic field vector, `data` is the raw data for the site read from
a CSV file, `columns` is the list of column names,
`Δ_date` is the vector of date differences between the observations (in UTC) and the
start date (in UTC), and `model` indicates which part of `Y` we are updating, and how to update it,
via different methods of `set_fluxnet_ic!`.
"""
function set_fluxnet_ic!(
    Y,
    data,
    columns,
    Δ_date,
    model::ClimaLand.Soil.EnergyHydrology;
    val = -9999,
)
    # Determine which column index corresponds to which varname
    varnames = ("SWC_F_MDS_1", "TS_F_MDS_1", "TA_F")
    column_name_map = Dict(
        varname => findfirst(columns[:] .== varname) for varname in varnames
    )
    FT = eltype(Y.soil.ρe_int)
    tmp_ic =
        model.parameters.θ_r +
        (model.parameters.ν - model.parameters.θ_r) * FT(0.95)
    if isnothing(column_name_map["SWC_F_MDS_1"])
        θ_l_0 = tmp_ic
    elseif unique(data[:, column_name_map["SWC_F_MDS_1"]]) == val
        θ_l_0 = tmp_ic
    else
        θ_l_0 =
            min.(
                FT(
                    get_data_at_start_date(
                        data[:, column_name_map["SWC_F_MDS_1"]],
                        Δ_date;
                        preprocess_func = x -> x / 100,
                        val,
                    ),
                ),
                tmp_ic,
            )
    end
    Y.soil.ϑ_l .= θ_l_0
    Y.soil.θ_i .= 0

    if isnothing(column_name_map["TS_F_MDS_1"])
        T_soil_0 = get_data_at_start_date(
            data[:, column_name_map["TA_F"]],
            Δ_date;
            preprocess_func = x -> x + 273.15,
            val,
        )
    elseif unique(data[:, column_name_map["TS_F_MDS_1"]]) == [val]
        T_soil_0 = get_data_at_start_date(
            data[:, column_name_map["TA_F"]],
            Δ_date;
            preprocess_func = x -> x + 273.15,
            val,
        )
    else
        T_soil_0 = get_data_at_start_date(
            data[:, column_name_map["TS_F_MDS_1"]],
            Δ_date;
            preprocess_func = x -> x + 273.15,
            val,
        )
    end

    ρc_s =
        ClimaLand.Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            model.parameters.ρc_ds,
            model.parameters.earth_param_set,
        )
    Y.soil.ρe_int =
        ClimaLand.Soil.volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            FT(T_soil_0),
            model.parameters.earth_param_set,
        )
end

"""
    set_fluxnet_ic!(Y, data, columns, Δ_date, model::ClimaLand.Canopy.CanopyModel)

Sets Y.canopy.energy.T to the air temperature at the observation date closest to the start
date of the model; sets the potential in the stem and leaf to -0.1 and -0.2 MPa, respectively,
and the computes the resulting water content Y.canopy.hydraulics.ϑ_l using the retention curve
of the plant.

If ony a leaf compartment is used, only the leaf ψ is used.
"""
function set_fluxnet_ic!(
    Y,
    data,
    columns,
    Δ_date,
    model::ClimaLand.Canopy.CanopyModel;
    val = -9999,
)
    # Determine which column index corresponds to air temperature
    idx = findfirst(columns[:] .== "TA_F")
    T_air_0 = get_data_at_start_date(
        data[:, idx],
        Δ_date;
        preprocess_func = x -> x + 273.15,
        val,
    )

    Y.canopy.energy.T .= T_air_0
    FT = eltype(Y.canopy.energy.T)
    ψ_stem_0 = FT(-1e5 / 9800) # pressure in the leaf divided by rho_liquid*gravitational acceleration [m]
    ψ_leaf_0 = FT(-2e5 / 9800)
    hydraulics = model.hydraulics
    n_stem = hydraulics.n_stem
    n_leaf = hydraulics.n_leaf
    ψ_comps = n_stem > 0 ? [ψ_stem_0, ψ_leaf_0] : ψ_leaf_0
    S_l_ini =
        ClimaLand.Canopy.PlantHydraulics.inverse_water_retention_curve.(
            hydraulics.parameters.retention_model,
            ψ_comps,
            hydraulics.parameters.ν,
            hydraulics.parameters.S_s,
        )
    for i in 1:(n_stem + n_leaf)
        Y.canopy.hydraulics.ϑ_l.:($i) .=
            ClimaLand.Canopy.PlantHydraulics.augmented_liquid_fraction.(
                hydraulics.parameters.ν,
                S_l_ini[i],
            )
    end
end

"""
    set_fluxnet_ic!(Y, data, columns, Δ_date, model::ClimaLand.Snow.SnowModel)

Sets Y.snow.S, Y.snow.S_l, and Y.snow.U in place to be zero at the start of the simulation
(no snow).

Note that the Snow NeuralDensity model has additional prognostic variables which also must be set
to zero; another method may work well for that case.
"""
function set_fluxnet_ic!(
    Y,
    data,
    columns,
    Δ_date,
    model::ClimaLand.Snow.SnowModel,
)
    Y.snow.S .= 0.0
    Y.snow.S_l .= 0.0
    Y.snow.U .= 0.0
end

"""
     set_fluxnet_ic!(Y, data, columns, Δ_date, model::ClimaLand.Soil.Biogeochemistry.SoilCO2Model)

Sets Y.soilco2.CO2, Y.soilco2.O2_f, and Y.soilco2.SOC in place with initial values:
- CO2: atmospheric CO2 concentration in mol co2 per mol air
- O2_f: volumetric fraction of O2 in soil air (~0.21)
- SOC: soil organic carbon concentration (kg C/m³)
"""
function set_fluxnet_ic!(
    Y,
    data,
    columns,
    Δ_date,
    model::ClimaLand.Soil.Biogeochemistry.SoilCO2Model,
)
    Y.soilco2.CO2 .= 0.000412
    Y.soilco2.O2_f .= 0.21  # Atmospheric O2 volumetric fraction
    Y.soilco2.SOC .= 5.0    # Default SOC concentration (kg C/m³)
end
