function run(ID; paramaters = parameters, domain = domain) # these args do not exist, for now suggestion 

    # Now we set up the model. For the soil model, we pick
    # a model type and model args:
    soil_domain = land_domain
    soil_ps = Soil.EnergyHydrologyParameters{FT}(;
        κ_dry = κ_dry,
        κ_sat_frozen = κ_sat_frozen,
        κ_sat_unfrozen = κ_sat_unfrozen,
        ρc_ds = ρc_ds,
        ν = soil_ν,
        ν_ss_om = ν_ss_om,
        ν_ss_quartz = ν_ss_quartz,
        ν_ss_gravel = ν_ss_gravel,
        hydrology_cm = vanGenuchten(; α = soil_vg_α, n = soil_vg_n),
        K_sat = soil_K_sat,
        S_s = soil_S_s,
        θ_r = θ_r,
        earth_param_set = earth_param_set,
        z_0m = z_0m_soil,
        z_0b = z_0b_soil,
        emissivity = soil_ϵ,
        PAR_albedo = soil_α_PAR,
        NIR_albedo = soil_α_NIR,
    );

    soil_args = (domain = soil_domain, parameters = soil_ps)
    soil_model_type = Soil.EnergyHydrology{FT}

    # Soil microbes model
    soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}

    soilco2_ps = SoilCO2ModelParameters{FT}(;
        ν = soil_ν, # same as soil
        θ_a100 = θ_a100,
        D_ref = D_ref,
        b = b,
        D_liq = D_liq,
        # DAMM
        α_sx = α_sx,
        Ea_sx = Ea_sx,
        kM_sx = kM_sx,
        kM_o2 = kM_o2,
        O2_a = O2_a,
        D_oa = D_oa,
        p_sx = p_sx,
        earth_param_set = earth_param_set,
    )

    # soil microbes args
    Csom = (z, t) -> eltype(z)(5.0)

    soilco2_top_bc = Soil.Biogeochemistry.SoilCO2StateBC((p, t) -> atmos_co2(t))
    soilco2_bot_bc = Soil.Biogeochemistry.SoilCO2FluxBC((p, t) -> 0.0) # no flux
    soilco2_sources = (MicrobeProduction{FT}(),)

    soilco2_boundary_conditions =
        (; top = (CO2 = soilco2_top_bc,), bottom = (CO2 = soilco2_bot_bc,))

    soilco2_drivers = Soil.Biogeochemistry.SoilDrivers(
        Soil.Biogeochemistry.PrognosticMet{FT}(),
        Soil.Biogeochemistry.PrescribedSOC{FT}(Csom),
        atmos,
    )

    soilco2_args = (;
        boundary_conditions = soilco2_boundary_conditions,
        sources = soilco2_sources,
        domain = soil_domain,
        parameters = soilco2_ps,
        drivers = soilco2_drivers,
    )

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
            ne = ne,
            ηsl = ηsl,
            σl = σl,
            μr = μr,
            μs = μs,
            f1 = f1,
            f2 = f2,
        )
    )
    # Set up radiative transfer
    radiative_transfer_args = (;
        parameters = TwoStreamParameters{FT}(;
            Ω = Ω,
            ld = ld,
            α_PAR_leaf = α_PAR_leaf,
            λ_γ_PAR = λ_γ_PAR,
            λ_γ_NIR = λ_γ_NIR,
            τ_PAR_leaf = τ_PAR_leaf,
            α_NIR_leaf = α_NIR_leaf,
            τ_NIR_leaf = τ_NIR_leaf,
            n_layers = n_layers,
            ϵ_canopy = ϵ_canopy,
        )
    )
    # Set up conductance
    conductance_args = (;
        parameters = MedlynConductanceParameters{FT}(;
            g1 = g1,
            Drel = Drel,
            g0 = g0,
        )
    )
    # Set up photosynthesis
    photosynthesis_args = (;
        parameters = FarquharParameters{FT}(
            Canopy.C3();
            oi = oi,
            ϕ = ϕ,
            θj = θj,
            f = f,
            sc = sc,
            pc = pc,
            Vcmax25 = Vcmax25,
            Γstar25 = Γstar25,
            Kc25 = Kc25,
            Ko25 = Ko25,
            To = To,
            ΔHkc = ΔHkc,
            ΔHko = ΔHko,
            ΔHVcmax = ΔHVcmax,
            ΔHΓstar = ΔHΓstar,
            ΔHJmax = ΔHJmax,
            ΔHRd = ΔHRd,
        )
    )
    # Set up plant hydraulics
    ai_parameterization = PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

    function root_distribution(z::T; rooting_depth = rooting_depth) where {T}
        return T(1.0 / rooting_depth) * exp(z / T(rooting_depth)) # 1/m
    end

    plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
        ai_parameterization = ai_parameterization,
        ν = plant_ν,
        S_s = plant_S_s,
        root_distribution = root_distribution,
        conductivity_model = conductivity_model,
        retention_model = retention_model,
    )
    plant_hydraulics_args = (
        parameters = plant_hydraulics_ps,
        n_stem = n_stem,
        n_leaf = n_leaf,
        compartment_midpoints = compartment_midpoints,
        compartment_surfaces = compartment_surfaces,
    )

    energy_args = (parameters = Canopy.BigLeafEnergyParameters{FT}(ac_canopy),)

    # Canopy component args
    canopy_component_args = (;
        autotrophic_respiration = autotrophic_respiration_args,
        radiative_transfer = radiative_transfer_args,
        photosynthesis = photosynthesis_args,
        conductance = conductance_args,
        hydraulics = plant_hydraulics_args,
        energy = energy_args,
    )

    # Other info needed
    shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
        z0_m,
        z0_b,
        earth_param_set,
    )

    canopy_model_args = (; parameters = shared_params, domain = canopy_domain)

    # Integrated plant hydraulics and soil model
    land_input = (atmos = atmos, radiation = radiation)
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

    #Initial conditions
    Y.soil.ϑ_l =
        drivers.SWC.status != absent ?
        drivers.SWC.values[1 + Int(round(t0 / DATA_DT))] : soil_ν / 2 # Get soil water content at t0
    # Both data and simulation are reference to 2005-01-01-00 (LOCAL)
    # or 2005-01-01-06 (UTC)
    Y.soil.θ_i = FT(0.0)
    T_0 =
        drivers.TS.status != absent ?
        drivers.TS.values[1 + Int(round(t0 / DATA_DT))] :
        drivers.TA.values[1 + Int(round(t0 / DATA_DT))] + 40# Get soil temperature at t0
    ρc_s =
        volumetric_heat_capacity.(Y.soil.ϑ_l, Y.soil.θ_i, Ref(land.soil.parameters))
    Y.soil.ρe_int =
        volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T_0,
            Ref(land.soil.parameters),
        )
    Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
    ψ_stem_0 = FT(-1e5 / 9800) # pressure in the leaf divided by rho_liquid*gravitational acceleration [m] 
    ψ_leaf_0 = FT(-2e5 / 9800)
    ψ_comps = n_stem > 0 ? [ψ_stem_0, ψ_leaf_0] : ψ_leaf_0

    S_l_ini =
        inverse_water_retention_curve.(retention_model, ψ_comps, plant_ν, plant_S_s)

    for i in 1:(n_stem + n_leaf)
        Y.canopy.hydraulics.ϑ_l.:($i) .=
            augmented_liquid_fraction.(plant_ν, S_l_ini[i])
    end

    Y.canopy.energy.T = drivers.TA.values[1 + Int(round(t0 / DATA_DT))] # Get atmos temperature at t0

    set_initial_cache! = make_set_initial_cache(land)
    set_initial_cache!(p, Y, t0);

    # Simulation
    sv = (;
        t = Array{FT}(undef, length(saveat)),
        saveval = Array{NamedTuple}(undef, length(saveat)),
    )
    cb = ClimaLSM.NonInterpSavingCallback(sv, saveat)

    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction((T_exp!) = exp_tendency!),
        Y,
        (t0, tf),
        p,
    );
    sol = SciMLBase.solve(
        prob,
        ode_algo;
        dt = dt,
        callback = cb,
        adaptive = false,
        saveat = saveat,
    )

    if isfile(
        joinpath(
            climalsm_dir,
            "experiments/integrated/fluxnet/$site_ID/Artifacts.toml",
        ),
    )
        rm(
            joinpath(
                climalsm_dir,
                "experiments/integrated/fluxnet/$site_ID/Artifacts.toml",
            ),
        )
    else
        nothing
    end

    return sv, sol, Y, p

end
