function soil_integrator(;
    t0 = 0.0,
    tf = 86400.0,
    dt = 450.0,
    stepper = CTS.ARS111(),
    dlat_degrees = 1.0,
    FT = Float32,
    n_vertical_elements = 15,
    diagnostics = false,
    update_drivers = false,
    info = true,
)

    context = ClimaComms.context()
    ClimaComms.init(context)
    device = ClimaComms.device()
    info && @info "Running on $device"

    n_horizontal_elements, effective_resolution, num_columns =
        resolution(; dlat_degrees)
    info &&
        @info "Running with $(n_horizontal_elements) horizontal elements ($(round(effective_resolution, sigdigits = 2)) degrees, $num_columns columns)"

    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    radius = FT(6378.1e3)
    depth = FT(3.5)
    domain = ClimaLand.Domains.SphericalShell(;
        radius = radius,
        depth = depth,
        nelements = (n_horizontal_elements, n_vertical_elements),
        npolynomial = 0,
        dz_tuple = FT.((1.0, 0.05)),
    )
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    start_date = DateTime(2008)
    stop_date = start_date + Second(tf - t0)
    time_interpolation_method = LinearInterpolation(Throw())
    # Forcing data
    era5_ncdata_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_path(; context)
    forcing = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        domain.space.surface,
        start_date,
        earth_param_set,
        FT;
        max_wind_speed = 25.0,
        time_interpolation_method,
    )
    model = ClimaLand.Soil.EnergyHydrology{FT}(domain, forcing, toml_dict)
    function set_ic!(Y, p, t, soil)
        soil_params = soil.parameters
        θ_r = soil_params.θ_r
        ν = soil_params.ν
        init_soil(ν, θ_r) = θ_r + (ν - θ_r) / 2
        Y.soil.ϑ_l .= init_soil.(ν, θ_r)
        Y.soil.θ_i .= FT(0.0)
        T = FT(276.85)
        ρc_s =
            Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l,
                Y.soil.θ_i,
                soil_params.ρc_ds,
                soil_params.earth_param_set,
            )
        Y.soil.ρe_int .=
            Soil.volumetric_internal_energy.(
                Y.soil.θ_i,
                ρc_s,
                T,
                soil_params.earth_param_set,
            )
    end

    sim = ClimaLand.Simulations.LandSimulation(
        start_date,
        stop_date,
        dt,
        model;
        diagnostics = (),
        updateat = update_drivers ? Second(dt) : Second(2 * (tf - t0)), # we still want to update drivers on init
        user_callbacks = (),
        set_ic!,
    )
    return sim
end
