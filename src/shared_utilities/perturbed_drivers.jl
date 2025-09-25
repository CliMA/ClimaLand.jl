export prescribed_perturbed_temperature_era5, prescribed_perturbed_rh_era5

"""
    perturbed_temp_specific_humidity_from_dewpoint(dewpoint_temperature, temperature, air_pressure, earth_param_set, ΔT)

Estimates the perturbed specific humidity given the true dewpoint temperature, true temperature of the air
in Kelvin, and true air pressure in Pa, along with a perturbation to the temperature ΔT, and the ClimaLand
`earth_param_set`. We first compute relative humidity, and then using RH and the perturbed air temperature,
equal to  true temperature+ΔT, compute a perturbed `q`.

Relative humidity is computed using the Magnus formula.

For more information on the Magnus formula, see e.g.
Lawrence, Mark G. (1 February 2005).
"The Relationship between Relative Humidity and the Dewpoint Temperature in Moist Air:
A Simple Conversion and Applications".
Bulletin of the American Meteorological Society. 86 (2): 225–234.
"""
function perturbed_temp_specific_humidity_from_dewpoint(
    T_dew_air::data_FT,
    T_air::data_FT,
    P_air::data_FT,
    earth_param_set,
    ΔT,
) where {data_FT <: Real}
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    _T_freeze = LP.T_freeze(earth_param_set)
    sim_FT = typeof(_T_freeze)
    # Obtain the relative humidity. This function requires temperatures in Celsius
    rh::sim_FT = rh_from_dewpoint(
        sim_FT(T_dew_air) - _T_freeze,
        sim_FT(T_air) - _T_freeze,
    )
    q = Thermodynamics.q_vap_from_RH_liquid(
        thermo_params,
        sim_FT(P_air),
        sim_FT(T_air + ΔT),
        rh,
    )
    return q
end


"""
    perturbed_rh_specific_humidity_from_dewpoint(dewpoint_temperature, temperature, air_pressure, earth_param_set, Δrh)

Estimates the perturbed specific humidity given the true dewpoint temperature, true temperature of the air
in Kelvin, and true air pressure in Pa, along with a perturbation to the relative humidity Δrh, and the ClimaLand
`earth_param_set`. We first compute relative humidity, and then perturb it to rh + Δrh, compute a perturbed `q`.

Relative humidity is computed using the Magnus formula.

For more information on the Magnus formula, see e.g.
Lawrence, Mark G. (1 February 2005).
"The Relationship between Relative Humidity and the Dewpoint Temperature in Moist Air:
A Simple Conversion and Applications".
Bulletin of the American Meteorological Society. 86 (2): 225–234.
"""
function perturbed_rh_specific_humidity_from_dewpoint(
    T_dew_air::data_FT,
    T_air::data_FT,
    P_air::data_FT,
    earth_param_set,
    Δrh,
) where {data_FT <: Real}
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    _T_freeze = LP.T_freeze(earth_param_set)
    sim_FT = typeof(_T_freeze)
    # Obtain the relative humidity. This function requires temperatures in Celsius
    rh_true::sim_FT = rh_from_dewpoint(
        sim_FT(T_dew_air) - _T_freeze,
        sim_FT(T_air) - _T_freeze,
    )
    rh = min(max(rh_true + Δrh, sqrt(eps(sim_FT))), sim_FT(1))
    q = Thermodynamics.q_vap_from_RH_liquid(
        thermo_params,
        sim_FT(P_air),
        sim_FT(T_air),
        rh,
    )
    return q
end

"""
    prescribed_perturbed_temperature_era5(era5_ncdata_path,
                             surface_space,
                             start_date,
                             toml_dict::CP.ParamDict,
                             ΔT,
                             FT;
                             gustiness=1,
                             max_wind_speed = nothing,
                             c_co2 = TimeVaryingInput((t) -> 4.2e-4),
                             time_interpolation_method = LinearInterpolation(PeriodicCalendar(Dates.Year(1), DateTime(Dates.year(stop_date)))),
                             regridder_type = :InterpolationsRegridder,
                             interpolation_method = Interpolations.Constant(),)

A helper function which constructs the `PrescribedAtmosphere` and `PrescribedRadiativeFluxes`
from a file path pointing to the ERA5 data in a netcdf file, the surface_space, the start date,
and the `toml_dict`, applying a change in the instantaneous temperature at each point
of ΔT, while keeping the relative humidity fixed. The LW_d shifts as LW_d -> LW_d + aΔT, with
a= 2W/m^2/K a surface climate sensitivity parameter.

Please see the documentation for `prescribed_forcing_era5` for more information.
"""
function prescribed_perturbed_temperature_era5(
    era5_ncdata_path,
    surface_space,
    start_date,
    toml_dict::CP.ParamDict,
    ΔT,
    FT;
    gustiness = 1,
    max_wind_speed = nothing,
    c_co2 = TimeVaryingInput((t) -> 4.2e-4),
    time_interpolation_method = LinearInterpolation(
        PeriodicCalendar(Dates.Year(1), DateTime(Dates.year(stop_date))),
    ),
    regridder_type = :InterpolationsRegridder,
    interpolation_method = Interpolations.Constant(),
)
    earth_param_set = LP.LandParameters(toml_dict)
    # Pass a list of files in all cases
    era5_ncdata_path isa String && (era5_ncdata_path = [era5_ncdata_path])
    _ρ_liq = LP.ρ_cloud_liq(earth_param_set)
    # Precip is provided as a mass flux; convert to volume flux of liquid water with ρ = 1000 kg/m^3
    precip = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path],
        ["mtpr", "msr"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (; preprocess_func = (data) -> -data / _ρ_liq,),
        method = time_interpolation_method,
        compose_function = (mtpr, msr) -> min.(mtpr .- msr, Float32(0)),
    )
    # Precip is provided as a mass flux; convert to volume flux of liquid water with ρ = 1000 kg/m^3
    snow_precip = TimeVaryingInput(
        era5_ncdata_path,
        "msr",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (; preprocess_func = (data) -> -data / _ρ_liq,),
        method = time_interpolation_method,
    )
    if max_wind_speed isa Nothing
        wind_function = (u, v) -> sqrt.(u .^ 2 .+ v .^ 2)
    else
        wind_function = (u, v) -> min.(sqrt.(u .^ 2 .+ v .^ 2), max_wind_speed)
    end

    u_atmos = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path],
        ["u10", "v10"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        compose_function = wind_function,
        method = time_interpolation_method,
    )
    specific_humidity(Td, T, P; params = earth_param_set, ΔT = ΔT) =
        ClimaLand.perturbed_temp_specific_humidity_from_dewpoint.(
            Td,
            T,
            P,
            params,
            ΔT,
        )
    q_atmos = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path, era5_ncdata_path],
        ["d2m", "t2m", "sp"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        compose_function = specific_humidity,
        method = time_interpolation_method,
    )
    P_atmos = TimeVaryingInput(
        era5_ncdata_path,
        "sp",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
    )
    perturb_temp(T; ΔT = ΔT) = T .+ ΔT
    T_atmos = TimeVaryingInput(
        era5_ncdata_path,
        "t2m",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (; preprocess_func = perturb_temp),
        method = time_interpolation_method,
    )
    h_atmos = FT(10)

    atmos = PrescribedAtmosphere(
        precip,
        snow_precip,
        T_atmos,
        u_atmos,
        q_atmos,
        P_atmos,
        start_date,
        h_atmos,
        toml_dict;
        gustiness = FT(gustiness),
        c_co2 = c_co2,
    )

    SW_d = TimeVaryingInput(
        era5_ncdata_path,
        "msdwswrf",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
    )
    function compute_diffuse_fraction(total, direct)
        diff = max(total - direct, Float32(0))
        return min(diff / (total + eps(Float32)), Float32(1))
    end
    function compute_diffuse_fraction_broadcasted(total, direct)
        return @. compute_diffuse_fraction(total, direct)
    end

    frac_diff = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path],
        ["msdwswrf", "msdrswrf"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
        compose_function = compute_diffuse_fraction_broadcasted,
    )
    perturb_lw_d(LW_d; ΔT = ΔT, a = 2) = LW_d .+ (a * ΔT)
    LW_d = TimeVaryingInput(
        era5_ncdata_path,
        "msdwlwrf",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
        file_reader_kwargs = (; preprocess_func = perturb_lw_d),
    )
    zenith_angle =
        (t, s) -> default_zenith_angle(
            t,
            s;
            latitude = ClimaCore.Fields.coordinate_field(surface_space).lat,
            longitude = ClimaCore.Fields.coordinate_field(surface_space).long,
            insol_params = earth_param_set.insol_params,
        )

    radiation = PrescribedRadiativeFluxes(
        FT,
        SW_d,
        LW_d,
        start_date;
        θs = zenith_angle,
        toml_dict = toml_dict,
        frac_diff = frac_diff,
    )
    return (; atmos, radiation)
end


"""
     prescribed_perturbed_rh_era5(era5_ncdata_path,
                             surface_space,
                             start_date,
                             toml_dict::CP.ParamDict,
                             Δrh,
                             FT;
                             gustiness=1,
                             max_wind_speed = nothing,
                             c_co2 = TimeVaryingInput((t) -> 4.2e-4),
                             time_interpolation_method = LinearInterpolation(PeriodicCalendar(Dates.Year(1), DateTime(Dates.year(stop_date)))),
                             regridder_type = :InterpolationsRegridder,
                             interpolation_method = Interpolations.Constant(),)

A helper function which constructs the `PrescribedAtmosphere` and `PrescribedRadiativeFluxes`
from a file path pointing to the ERA5 data in a netcdf file, the surface_space, the start date,
and the `toml_dict`, applying a change in the instantaneous change to relative
humidity at each point
of Δrh. The perturbed rh is clipped to be within the range (0,1].

Please see the documentation for `prescribed_forcing_era5` for more information.
"""
function prescribed_perturbed_rh_era5(
    era5_ncdata_path,
    surface_space,
    start_date,
    toml_dict::CP.ParamDict,
    Δrh,
    FT;
    gustiness = 1,
    max_wind_speed = nothing,
    c_co2 = TimeVaryingInput((t) -> 4.2e-4),
    time_interpolation_method = LinearInterpolation(
        PeriodicCalendar(Dates.Year(1), DateTime(Dates.year(stop_date))),
    ),
    regridder_type = :InterpolationsRegridder,
    interpolation_method = Interpolations.Constant(),
)
    earth_param_set = LP.LandParameters(toml_dict)
    # Pass a list of files in all cases
    era5_ncdata_path isa String && (era5_ncdata_path = [era5_ncdata_path])
    _ρ_liq = LP.ρ_cloud_liq(earth_param_set)
    # Precip is provided as a mass flux; convert to volume flux of liquid water with ρ = 1000 kg/m^3
    precip = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path],
        ["mtpr", "msr"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (; preprocess_func = (data) -> -data / _ρ_liq,),
        method = time_interpolation_method,
        compose_function = (mtpr, msr) -> min.(mtpr .- msr, Float32(0)),
    )
    # Precip is provided as a mass flux; convert to volume flux of liquid water with ρ = 1000 kg/m^3
    snow_precip = TimeVaryingInput(
        era5_ncdata_path,
        "msr",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (; preprocess_func = (data) -> -data / _ρ_liq,),
        method = time_interpolation_method,
    )
    if max_wind_speed isa Nothing
        wind_function = (u, v) -> sqrt.(u .^ 2 .+ v .^ 2)
    else
        wind_function = (u, v) -> min.(sqrt.(u .^ 2 .+ v .^ 2), max_wind_speed)
    end

    u_atmos = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path],
        ["u10", "v10"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        compose_function = wind_function,
        method = time_interpolation_method,
    )
    specific_humidity(Td, T, P; params = earth_param_set, Δrh = Δrh) =
        ClimaLand.perturbed_rh_specific_humidity_from_dewpoint.(
            Td,
            T,
            P,
            params,
            Δrh,
        )
    q_atmos = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path, era5_ncdata_path],
        ["d2m", "t2m", "sp"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        compose_function = specific_humidity,
        method = time_interpolation_method,
    )
    P_atmos = TimeVaryingInput(
        era5_ncdata_path,
        "sp",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
    )
    T_atmos = TimeVaryingInput(
        era5_ncdata_path,
        "t2m",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
    )
    h_atmos = FT(10)

    atmos = PrescribedAtmosphere(
        precip,
        snow_precip,
        T_atmos,
        u_atmos,
        q_atmos,
        P_atmos,
        start_date,
        h_atmos,
        toml_dict;
        gustiness = FT(gustiness),
        c_co2 = c_co2,
    )

    SW_d = TimeVaryingInput(
        era5_ncdata_path,
        "msdwswrf",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
    )
    function compute_diffuse_fraction(total, direct)
        diff = max(total - direct, Float32(0))
        return min(diff / (total + eps(Float32)), Float32(1))
    end
    function compute_diffuse_fraction_broadcasted(total, direct)
        return @. compute_diffuse_fraction(total, direct)
    end

    frac_diff = TimeVaryingInput(
        [era5_ncdata_path, era5_ncdata_path],
        ["msdwswrf", "msdrswrf"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
        compose_function = compute_diffuse_fraction_broadcasted,
    )
    LW_d = TimeVaryingInput(
        era5_ncdata_path,
        "msdwlwrf",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
    )
    zenith_angle =
        (t, s) -> default_zenith_angle(
            t,
            s;
            latitude = ClimaCore.Fields.coordinate_field(surface_space).lat,
            longitude = ClimaCore.Fields.coordinate_field(surface_space).long,
            insol_params = earth_param_set.insol_params,
        )

    radiation = PrescribedRadiativeFluxes(
        FT,
        SW_d,
        LW_d,
        start_date;
        θs = zenith_angle,
        toml_dict = toml_dict,
        frac_diff = frac_diff,
    )
    return (; atmos, radiation)
end
