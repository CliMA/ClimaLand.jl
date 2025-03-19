import StaticArrays: SVector
using ClimaLand.Canopy.PlantHydraulics

function lighter_canopy_integrator(;
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
    domain = ClimaLand.Domains.SphericalSurface(;
        radius = radius,
        nelements = n_horizontal_elements,
        npolynomial = 0,
    )
    surface_space = domain.space.surface

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

    rooting_depth = FT(0.5)
    Vcmax25 = FT(5e-5)
    is_c3 = FT(1.0)

    # Energy Balance model
    ac_canopy = FT(2.5e3)
    # Plant Hydraulics and general plant parameters
    SAI = FT(0.0) # m2/m2
    f_root_to_shoot = FT(3.5)
    RAI = FT(1.0)
    K_sat_plant = FT(5e-9) # m/s # seems much too small?
    ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
    Weibull_param = FT(4) # unitless, Holtzman's original c param value
    a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
    conductivity_model =
        Canopy.PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
    retention_model = Canopy.PlantHydraulics.LinearRetentionCurve{FT}(a)
    plant_ν = FT(1.44e-4)
    plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
    n_stem = 0
    n_leaf = 1
    h_stem = FT(0.0)
    h_leaf = FT(1.0)
    zmax = FT(0.0)
    h_canopy = h_stem + h_leaf
    compartment_midpoints =
        n_stem > 0 ? [h_stem / 2, h_stem + h_leaf / 2] : [h_leaf / 2]
    compartment_surfaces =
        n_stem > 0 ? [zmax, h_stem, h_canopy] : [zmax, h_leaf]

    z0_m = FT(0.13) * h_canopy
    z0_b = FT(0.1) * z0_m

    # Now we set up the canopy model, which we set up by component:
    # Component Types
    AR_params = Canopy.AutotrophicRespirationParameters(FT)
    AR_model = AutotrophicRespirationModel{FT}(AR_params)

    rt_params = TwoStreamParameters(
        FT;
        G_Function = ConstantGFunction(FT(0.5)),
        α_PAR_leaf = FT(0.1),
        α_NIR_leaf = FT(0.45),
        τ_PAR_leaf = FT(0.05),
        τ_NIR_leaf = FT(0.25),
        Ω = FT(0.69),
        λ_γ_PAR = FT(5e-7),
    )
    rt_model = TwoStreamModel{FT}(rt_params)
    cond_params = MedlynConductanceParameters(FT; g1 = FT(141.0))
    stomatal_model = MedlynConductanceModel{FT}(cond_params)

    photo_params = Canopy.FarquharParameters(FT, is_c3; Vcmax25 = Vcmax25)
    photosynthesis_model = FarquharModel{FT}(photo_params)

    # Set up plant hydraulics
    modis_lai_artifact_path =
        ClimaLand.Artifacts.modis_lai_forcing_data_path(; context)
    modis_lai_ncdata_path =
        joinpath(modis_lai_artifact_path, "Yuan_et_al_2008_1x1.nc")
    LAIfunction = ClimaLand.prescribed_lai_modis(
        modis_lai_ncdata_path,
        surface_space,
        start_date;
        time_interpolation_method = time_interpolation_method,
    )
    ai_parameterization =
        Canopy.PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

    plant_hydraulics_ps = Canopy.PlantHydraulics.PlantHydraulicsParameters(;
        ai_parameterization = ai_parameterization,
        ν = plant_ν,
        S_s = plant_S_s,
        rooting_depth = rooting_depth,
        conductivity_model = conductivity_model,
        retention_model = retention_model,
    )

    plant_hydraulics = Canopy.PlantHydraulicsModel{FT}(
        parameters = plant_hydraulics_ps,
        n_stem = n_stem,
        n_leaf = n_leaf,
        compartment_midpoints = compartment_midpoints,
        compartment_surfaces = compartment_surfaces,
    )
    energy = Canopy.BigLeafEnergyModel{FT}(
        Canopy.BigLeafEnergyParameters{FT}(ac_canopy),
    )

    shared_params = Canopy.SharedCanopyParameters{FT, typeof(earth_param_set)}(
        z0_m,
        z0_b,
        earth_param_set,
    )

    ψ_soil0 = FT(0.0)

    soil_driver = PrescribedGroundConditions(
        FT;
        root_depths = SVector{10, FT}(-(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0),
        ψ = t -> ψ_soil0,
        α_PAR = FT(0.2),
        α_NIR = FT(0.4),
        T = t -> 298.0,
        ϵ = FT(0.99),
    )

    # Canopy component args
    canopy = ClimaLand.Canopy.CanopyModel{FT}(;
        parameters = shared_params,
        domain = domain,
        autotrophic_respiration = AR_model,
        radiative_transfer = rt_model,
        photosynthesis = photosynthesis_model,
        conductance = stomatal_model,
        hydraulics = plant_hydraulics,
        energy,
        boundary_conditions = Canopy.AtmosDrivenCanopyBC(
            atmos,
            radiation,
            soil_driver,
        ),
    )


    Y, p, cds = initialize(canopy)

    # # Provide initial conditions for the canopy hydraulics model

    # ψ_stem_0 = FT(-1e5 / 9800)
    # ψ_leaf_0 = FT(-2e5 / 9800)

    # S_l_ini =
    #     inverse_water_retention_curve.(
    #         retention_model,
    #         [ψ_stem_0, ψ_leaf_0],
    #         plant_ν,
    #         plant_S_s,
    #     )

    # for i in 1:2
    #     Y.canopy.hydraulics.ϑ_l.:($i) .= augmented_liquid_fraction.(plant_ν, S_l_ini[i])
    # end;

    set_initial_cache! = make_set_initial_cache(canopy)
    exp_tendency! = make_exp_tendency(canopy)
    imp_tendency! = ClimaLand.make_imp_tendency(canopy)
    jacobian! = ClimaLand.make_jacobian(canopy)
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
            canopy,
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
        drivers = ClimaLand.get_drivers(canopy)
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
