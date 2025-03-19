function richards_integrator(;
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
    spatially_varying_soil_params =
        ClimaLand.default_spatially_varying_soil_parameters(
            subsurface_space,
            surface_space,
            FT,
        )
    (; ν, hydrology_cm, K_sat, S_s, θ_r, f_max) = spatially_varying_soil_params
    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
        f_over = f_over,
        f_max = f_max,
        R_sb = R_sb,
    )
    soil_params = ClimaLand.Soil.RichardsParameters(;
        hydrology_cm = hydrology_cm,
        ν = ν,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
    )

    era5_artifact_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(; context)

    # Below, the preprocess_func argument is used to
    # 1. Convert precipitation to be negative (as it is downwards)
    # 2. Convert mass flux to equivalent liquid water flux
    # Precipitation:
    precip = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2008_1.0x1.0.nc"),
        "mtpr",
        surface_space;
        start_date,
        regridder_type = :InterpolationsRegridder,
        file_reader_kwargs = (; preprocess_func = (data) -> -data / 1000),
    )
    atmos = ClimaLand.PrescribedPrecipitation{FT, typeof(precip)}(precip)
    bottom_bc = ClimaLand.Soil.WaterFluxBC((p, t) -> 0.0)
    bc = (;
        top = ClimaLand.Soil.RichardsAtmosDrivenFluxBC(atmos, runoff_model),
        bottom = bottom_bc,
    )
    model = ClimaLand.Soil.RichardsModel{FT}(;
        parameters = soil_params,
        domain = domain,
        boundary_conditions = bc,
        sources = (),
        lateral_flow = false,
    )

    Y, p, cds = initialize(model)
    z = ClimaCore.Fields.coordinate_field(cds.subsurface).z
    lat = ClimaCore.Fields.coordinate_field(domain.space.subsurface).lat
    function hydrostatic_profile(
        lat::FT,
        z::FT,
        ν::FT,
        θ_r::FT,
        α::FT,
        n::FT,
        S_s::FT,
        fmax,
    ) where {FT}
        m = 1 - 1 / n
        zmin = FT(-50.0)
        zmax = FT(0.0)

        z_∇ = FT(zmin / 5.0 + (zmax - zmin) / 2.5 * (fmax - 0.35) / 0.7)
        if z > z_∇
            S = FT((FT(1) + (α * (z - z_∇))^n)^(-m))
            ϑ_l = S * (ν - θ_r) + θ_r
        else
            ϑ_l = -S_s * (z - z_∇) + ν
        end
        return FT(ϑ_l)
    end

    # Set initial state values
    vg_α = hydrology_cm.α
    vg_n = hydrology_cm.n
    Y.soil.ϑ_l .= hydrostatic_profile.(lat, z, ν, θ_r, vg_α, vg_n, S_s, f_max)
    # Create model update functions
    set_initial_cache! = make_set_initial_cache(model)
    exp_tendency! = make_exp_tendency(model)
    imp_tendency! = make_imp_tendency(model)
    jacobian! = make_jacobian(model)

    set_initial_cache!(p, Y, t0)

    # set up jacobian info
    jac_kwargs = (;
        jac_prototype = ClimaLand.FieldMatrixWithSolver(Y),
        Wfact = jacobian!,
    )

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
        drivers = ClimaLand.get_drivers(model)
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
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
            dss! = ClimaLand.dss!,
        ),
        Y,
        (t0, tf),
        p,
    )
    callback = SciMLBase.CallbackSet(callbacks...)

    integrator = SciMLBase.init(prob, ode_algo; dt, callback)
    return integrator
end
