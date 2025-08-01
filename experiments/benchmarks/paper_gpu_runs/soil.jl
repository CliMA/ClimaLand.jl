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
        FT,
    )

    spatially_varying_soil_params =
        ClimaLand.default_spatially_varying_soil_parameters(
            subsurface_space,
            surface_space,
            FT,
        )
    (;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        PAR_albedo_dry,
        NIR_albedo_dry,
        PAR_albedo_wet,
        NIR_albedo_wet,
        f_max,
    ) = spatially_varying_soil_params
    soil_params = Soil.EnergyHydrologyParameters(
        FT;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        PAR_albedo_dry = PAR_albedo_dry,
        NIR_albedo_dry = NIR_albedo_dry,
        PAR_albedo_wet = PAR_albedo_wet,
        NIR_albedo_wet = NIR_albedo_wet,
    )

    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
        f_over = f_over,
        f_max = f_max,
        R_sb = R_sb,
    )

    surface_water_flux = WaterFluxBC((p, t) -> .-K_sat ./ 10)
    bottom_water_flux = WaterFluxBC((p, t) -> 0.0)

    surface_heat_flux = HeatFluxBC((p, t) -> 0.0)
    bottom_heat_flux = HeatFluxBC((p, t) -> 0.0)

    boundary_fluxes = (;
        top = WaterHeatBC(;
            water = surface_water_flux,
            heat = surface_heat_flux,
        ),
        bottom = WaterHeatBC(;
            water = bottom_water_flux,
            heat = bottom_heat_flux,
        ),
    )


    # Spatially varying canopy parameters from CLM
    soilco2_ps = Soil.Biogeochemistry.SoilCO2ModelParameters(FT)

    soil = Soil.EnergyHydrology{FT}(;
        domain = domain,
        parameters = soil_params,
        boundary_conditions = boundary_fluxes,
        sources = (),
    )

    Y, p, cds = initialize(soil)

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

    set_initial_cache! = make_set_initial_cache(soil)
    exp_tendency! = make_exp_tendency(soil)
    imp_tendency! = ClimaLand.make_imp_tendency(soil)
    jacobian! = ClimaLand.make_jacobian(soil)
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
            soil,
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
        drivers = ClimaLand.get_drivers(soil)
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
