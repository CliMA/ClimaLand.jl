export run_fluxnet

"""
    fluxnet_simulation(site_ID; kwargs)

Run ClimaLand at a fluxnet site.
"""
function run_fluxnet(
    site_ID;
    FT = Float64,
    context = ClimaComms.SingletonCommsContext(),
    setup = make_setup(site_ID),
    domain = make_domain(setup, FT),
    config = make_config(site_ID),
    params = make_parameters(site_ID),
    drivers = make_drivers(site_ID, setup, config, params, context),
    timestepper = make_timestepper(setup),
)
    #earth_param_set = create_lsm_parameters(FT)

    # Now we set up the model. For the soil model, we pick
    # a model type and model args:
    soil_domain = domain.land_domain
    PAR_albedo =
        FT(0.5) * params.soil.PAR_albedo_wet +
        FT(0.5) * params.soil.PAR_albedo_dry
    NIR_albedo =
        FT(0.5) * params.soil.NIR_albedo_wet +
        FT(0.5) * params.soil.NIR_albedo_dry
    soil_ps = Soil.EnergyHydrologyParameters(
        FT;
        ν = params.soil.ν,
        ν_ss_om = params.soil.ν_ss_om,
        ν_ss_quartz = params.soil.ν_ss_quartz,
        ν_ss_gravel = params.soil.ν_ss_gravel,
        hydrology_cm = params.soil.hydrology_cm,
        K_sat = params.soil.K_sat,
        S_s = params.soil.S_s,
        θ_r = params.soil.θ_r,
        z_0m = params.soil.z_0m,
        z_0b = params.soil.z_0b,
        emissivity = params.soil.emissivity,
        PAR_albedo = PAR_albedo,
        NIR_albedo = NIR_albedo,
    )

    soil_args = (domain = soil_domain, parameters = soil_ps)
    soil_model_type = Soil.EnergyHydrology{FT}

    # Soil microbes model
    soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}

    soilco2_ps = SoilCO2ModelParameters(
        FT;
        D_ref = params.hetero_resp.D_ref,
        D_liq = params.hetero_resp.D_liq,
        # DAMM
        α_sx = params.hetero_resp.α_sx,
        Ea_sx = params.hetero_resp.Ea_sx,
        kM_sx = params.hetero_resp.kM_sx,
        kM_o2 = params.hetero_resp.kM_o2,
        O2_a = params.hetero_resp.O2_a,
        D_oa = params.hetero_resp.D_oa,
        p_sx = params.hetero_resp.p_sx,
        earth_param_set = earth_param_set,
    )

    # soil microbes args
    Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))
    soilco2_args = (; domain = soil_domain, parameters = soilco2_ps)

    # Now we set up the canopy model, which we set up by component:
    # Component Types
    canopy_component_types = (;
        autotrophic_respiration = Canopy.AutotrophicRespirationModel{FT},
        radiative_transfer = Canopy.TwoStreamModel{FT},
        photosynthesis = Canopy.FarquharModel{FT},
        conductance = Canopy.MedlynConductanceModel{FT},
        hydraulics = Canopy.PlantHydraulicsModel{FT},
        energy = Canopy.BigLeafEnergyModel{FT},
    )
    # Individual Component arguments
    # Set up autotrophic respiration
    autotrophic_respiration_args = (;
        parameters = AutotrophicRespirationParameters{FT}(;
            ne = params.auto_resp.ne,
            ηsl = params.auto_resp.ηsl,
            σl = params.auto_resp.σl,
            μr = params.auto_resp.μr,
            μs = params.auto_resp.μs,
            Rel = params.auto_resp.Rel,
        )
    )
    # Set up radiative transfer
    radiative_transfer_args = (;
        parameters = TwoStreamParameters(
            FT;
            Ω = params.radiative_transfer.Ω,
            G_Function = params.radiative_transfer.G_Function,
            α_PAR_leaf = params.radiative_transfer.α_PAR_leaf,
            λ_γ_PAR = params.radiative_transfer.λ_γ_PAR,
            τ_PAR_leaf = params.radiative_transfer.τ_PAR_leaf,
            α_NIR_leaf = params.radiative_transfer.α_NIR_leaf,
            τ_NIR_leaf = params.radiative_transfer.τ_NIR_leaf,
            n_layers = params.radiative_transfer.n_layers,
            ϵ_canopy = params.radiative_transfer.ϵ_canopy,
        )
    )
    # Set up conductance
    conductance_args = (;
        parameters = MedlynConductanceParameters(
            FT;
            g1 = params.conductance.g1,
            Drel = params.conductance.Drel,
            g0 = params.conductance.g0,
        )
    )
    # Set up photosynthesis
    photosynthesis_args = (;
        parameters = FarquharParameters(
            FT,
            FT(1);
            oi = params.photosynthesis.oi,
            ϕ = params.photosynthesis.ϕ,
            θj = params.photosynthesis.θj,
            fC3 = params.photosynthesis.fC3,
            sc = params.photosynthesis.sc,
            pc = params.photosynthesis.pc,
            Vcmax25 = params.photosynthesis.Vcmax25,
            Γstar25 = params.photosynthesis.Γstar25,
            Kc25 = params.photosynthesis.Kc25,
            Ko25 = params.photosynthesis.Ko25,
            To = params.photosynthesis.To,
            ΔHkc = params.photosynthesis.ΔHkc,
            ΔHko = params.photosynthesis.ΔHko,
            ΔHVcmax = params.photosynthesis.ΔHVcmax,
            ΔHΓstar = params.photosynthesis.ΔHΓstar,
            ΔHJmax = params.photosynthesis.ΔHJmax,
            ΔHRd = params.photosynthesis.ΔHRd,
            fC4 = params.photosynthesis.fC4,
            Q10 = params.photosynthesis.Q10,
            E = params.photosynthesis.E,
            s1 = params.photosynthesis.s1,
            s2 = params.photosynthesis.s2,
            s3 = params.photosynthesis.s3,
            s4 = params.photosynthesis.s4,
            s5 = params.photosynthesis.s5,
            s6 = params.photosynthesis.s6,
        )
    )
    # Set up plant hydraulics
    ai_parameterization = PrescribedSiteAreaIndex{FT}(
        drivers.LAIfunction,
        params.plant_hydraulics.SAI,
        drivers.RAI,
    )


    plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
        ai_parameterization = ai_parameterization,
        ν = drivers.plant_ν,
        S_s = params.plant_hydraulics.S_s,
        rooting_depth = params.plant_hydraulics.rooting_depth,
        conductivity_model = params.plant_hydraulics.conductivity_model,
        retention_model = params.plant_hydraulics.retention_model,
    )
    plant_hydraulics_args = (
        parameters = plant_hydraulics_ps,
        n_stem = setup.n_stem,
        n_leaf = setup.n_leaf,
        compartment_midpoints = domain.compartment_midpoints,
        compartment_surfaces = domain.compartment_surfaces,
    )

    energy_args = (
        parameters = Canopy.BigLeafEnergyParameters{FT}(
            params.canopy_energy_balance.ac_canopy,
        ),
    )

    # Canopy component args
    canopy_component_args = (;
        autotrophic_respiration = autotrophic_respiration_args,
        radiative_transfer = radiative_transfer_args,
        photosynthesis = photosynthesis_args,
        conductance = conductance_args,
        hydraulics = plant_hydraulics_args,
        energy = energy_args,
    )

    z0_m = FT(0.13) * domain.h_canopy
    z0_b = FT(0.1) * z0_m

    # Other info needed
    shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
        z0_m,
        z0_b,
        earth_param_set,
    )

    canopy_model_args =
        (; parameters = shared_params, domain = domain.canopy_domain)

    # Integrated plant hydraulics and soil model
    land_input = (
        atmos = drivers.atmos,
        radiation = drivers.radiation,
        soil_organic_carbon = Csom,
    )
    land = SoilCanopyModel{FT}(;
        soilco2_type = soilco2_type,
        soilco2_args = soilco2_args,
        land_args = land_input,
        soil_model_type = soil_model_type,
        soil_args = soil_args,
        canopy_component_types = canopy_component_types,
        canopy_component_args = canopy_component_args,
        canopy_model_args = canopy_model_args,
    )
    Y, p, cds = initialize(land)
    exp_tendency! = make_exp_tendency(land)
    imp_tendency! = make_imp_tendency(land)
    jacobian! = make_jacobian(land)
    jac_kwargs = (;
        jac_prototype = ClimaLand.FieldMatrixWithSolver(Y),
        Wfact = jacobian!,
    )

    #Initial conditions
    Y.soil.ϑ_l =
        drivers.drivers.SWC.status != absent ?
        drivers.drivers.SWC.values[1 + Int(round(setup.t0 / drivers.DATA_DT))] :
        params.soil.soil_ν / 2 # Get soil water content at t0
    # Both data and simulation are reference to 2005-01-01-00 (LOCAL)
    # or 2005-01-01-06 (UTC)
    Y.soil.θ_i = FT(0.0)
    T_0 =
        drivers.drivers.TS.status != absent ?
        drivers.drivers.TS.values[1 + Int(round(setup.t0 / drivers.DATA_DT))] :
        drivers.drivers.TA.values[1 + Int(round(setup.t0 / drivers.DATA_DT))] +
        40# Get soil temperature at t0
    ρc_s =
        volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            land.soil.parameters.ρc_ds,
            land.soil.parameters.earth_param_set,
        )
    Y.soil.ρe_int =
        volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T_0,
            land.soil.parameters.earth_param_set,
        )
    Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
    ψ_stem_0 = FT(-1e5 / 9800) # pressure in the leaf divided by rho_liquid*gravitational acceleration [m]
    ψ_leaf_0 = FT(-2e5 / 9800)
    ψ_comps = setup.n_stem > 0 ? [ψ_stem_0, ψ_leaf_0] : ψ_leaf_0

    S_l_ini =
        inverse_water_retention_curve.(
            params.plant_hydraulics.retention_model,
            ψ_comps,
            drivers.plant_ν,
            params.plant_hydraulics.S_s,
        )

    for i in 1:(setup.n_stem + setup.n_leaf)
        Y.canopy.hydraulics.ϑ_l.:($i) .=
            augmented_liquid_fraction.(drivers.plant_ν, S_l_ini[i])
    end

    Y.canopy.energy.T =
        drivers.drivers.TA.values[1 + Int(round(setup.t0 / drivers.DATA_DT))] # Get atmos temperature at t0

    set_initial_cache! = make_set_initial_cache(land)
    set_initial_cache!(p, Y, setup.t0)

    # Simulation
    sv = (;
        t = Array{Float64}(undef, length(timestepper.saveat)),
        saveval = Array{NamedTuple}(undef, length(timestepper.saveat)),
    )
    saving_cb = ClimaLand.NonInterpSavingCallback(sv, timestepper.saveat)
    ## How often we want to update the drivers. Note that this uses the defined `t0` and `tf`
    ## defined in the simulatons file
    updateat = Array((setup.t0):(drivers.DATA_DT):(timestepper.tf))
    land_drivers = ClimaLand.get_drivers(land)
    updatefunc = ClimaLand.make_update_drivers(land_drivers)
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    cb = SciMLBase.CallbackSet(driver_cb, saving_cb)

    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
            dss! = ClimaLand.dss!,
        ),
        Y,
        (setup.t0, timestepper.tf),
        p,
    )
    sol = SciMLBase.solve(
        prob,
        timestepper.ode_algo;
        dt = setup.dt,
        callback = cb,
        adaptive = false,
        saveat = timestepper.saveat,
    )

    if isfile(
        joinpath(
            climaland_dir,
            "lib/ClimaLandSimulations/src/utilities/$site_ID/Artifacts.toml",
        ),
    )
        rm(
            joinpath(
                climaland_dir,
                "lib/ClimaLandSimulations/src/utilities/$site_ID/Artifacts.toml",
            ),
        )
    else
        nothing
    end

    return sv, sol, Y, p
end
