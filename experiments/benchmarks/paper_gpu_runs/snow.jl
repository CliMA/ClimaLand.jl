function snow_integrator(;
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
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    subsurface_space = domain.space.subsurface

    start_date = DateTime(2008)
    stop_date = start_date + Second(tf - t0)
    time_interpolation_method = LinearInterpolation(Throw())
    # Forcing data
    era5_ncdata_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_path(; context)
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT;
        max_wind_speed = 25.0,
        time_interpolation_method,
    )
    forcing = (; atmos, radiation)
    snow = Snow.SnowModel(FT, surface_domain, forcing, toml_dict, dt;)


    function set_ic!(Y, p, t, model)
        SWE = [FT(1)]

        Y.snow.S .= FT(SWE[1]) # first data point
        Y.snow.U .= ClimaLand.Snow.energy_from_q_l_and_swe(
            FT(SWE[1]),
            FT(0),
            model.parameters,
        ) # with q_l = 0
    end

    sim = ClimaLand.Simulations.LandSimulation(
        start_date,
        stop_date,
        dt,
        snow;
        diagnostics = (),
        updateat = update_drivers ? Second(dt) : Second(2 * (tf - t0)), # we still want to update drivers on init
        user_callbacks = (),
        set_ic!,
    )
    return sim
end
