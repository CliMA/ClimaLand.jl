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

    earth_param_set = LP.LandParameters(FT)
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
    time_interpolation_method = LinearInterpolation(Throw())
    # Forcing data
    era5_artifact_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(; context)
    era5_ncdata_path = joinpath(era5_artifact_path, "era5_2008_1.0x1.0.nc")
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT;
        time_interpolation_method = time_interpolation_method,
        regridder_type = :InterpolationsRegridder,
    )

    snow_parameters = SnowParameters{FT}(dt; earth_param_set = earth_param_set)
    snow = Snow.SnowModel(
        parameters = snow_parameters,
        domain = ClimaLand.obtain_surface_domain(domain),
        boundary_conditions = ClimaLand.Snow.AtmosDrivenSnowBC(
            atmos,
            radiation,
        ),
    )


    Y, p, coords = ClimaLand.initialize(snow)
    SWE = [FT(1)]

    Y.snow.S .= FT(SWE[1]) # first data point
    Y.snow.U .= ClimaLand.Snow.energy_from_q_l_and_swe(
        FT(SWE[1]),
        FT(0),
        snow_parameters,
    ) # with q_l = 0

    set_initial_cache! = make_set_initial_cache(snow)
    set_initial_cache!(p, Y, t0)
    exp_tendency! = make_exp_tendency(snow)

    callbacks = tuple()

    if info
        walltime_info = WallTimeInfo()
        every10steps(u, t, integrator) = mod(integrator.step, 10) == 0
        report = let wt = walltime_info
            (integrator) -> report_walltime(wt, integrator)
        end
        report_cb = SciMLBase.DiscreteCallback(every10steps, report)
        callbacks = (callbacks..., report_cb)
    end

    if diagnostics
        outdir = mktempdir(pwd())
        info && @info "Output directory: $outdir"
        nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
            subsurface_space,
            outdir;
            start_date,
        )

        diags = ClimaLand.default_diagnostics(
            land,
            start_date;
            output_writer = nc_writer,
            output_vars = :short,
            average_period = :monthly,
        )

        diagnostic_handler =
            ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = dt)

        diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)
        callbacks = (callbacks..., diag_cb)
    end

    if update_drivers
        info && @info "Updating drivers every 3 hours"
        updateat = Array(t0:(3600 * 3):tf)
        drivers = ClimaLand.get_drivers(snow)
        updatefunc = ClimaLand.make_update_drivers(drivers)
        driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
        callbacks = (callbacks..., driver_cb)
    end

    ode_algo = CTS.IMEXAlgorithm(
        stepper,
        CTS.NewtonsMethod(
            max_iters = 3,
            update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        ),
    )
    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(T_exp! = exp_tendency!, dss! = ClimaLand.dss!),
        Y,
        (t0, tf),
        p,
    )
    callback = SciMLBase.CallbackSet(callbacks...)

    integrator = SciMLBase.init(prob, ode_algo; dt, callback)
    return integrator
end
