function bucket_integrator(;
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
    # Set up parameters
    σS_c = FT(0.2)
    W_f = FT(0.2)
    z_0m = FT(1e-3)
    z_0b = FT(1e-3)
    κ_soil = FT(1.5)
    ρc_soil = FT(2e6)
    τc = FT(dt)
    α_snow = FT(0.8)
    albedo = PrescribedBaregroundAlbedo{FT}(α_snow, surface_space)
    bucket_parameters = BucketModelParameters(FT; albedo, z_0m, z_0b, τc)
    bucket = BucketModel(
        parameters = bucket_parameters,
        domain = domain,
        atmosphere = atmos,
        radiation = radiation,
    )

    temp_anomaly_amip(coord) = 40 * cosd(coord.lat)^4
    Y, p, cds = initialize(bucket)
    # Set temperature IC including anomaly, based on atmospheric setup
    T_sfc_0 = FT(271.0)
    @. Y.bucket.T = T_sfc_0 + temp_anomaly_amip(cds.subsurface)
    Y.bucket.W .= FT(0.15)
    Y.bucket.Ws .= FT(0.0)
    Y.bucket.σS .= FT(0.0)

    set_initial_cache! = make_set_initial_cache(bucket)
    set_initial_cache!(p, Y, t0)
    exp_tendency! = make_exp_tendency(bucket)

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
        drivers = ClimaLand.get_drivers(bucket)
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
