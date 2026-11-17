import ClimaTimeSteppers as CTS
using ClimaCore
import ClimaComms
ClimaComms.@import_required_backends
using CairoMakie
using Statistics
using Dates

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaParams
import ClimaUtilities.OutputPathGenerator: generate_output_path
using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaTimeSteppers

global climaland_dir = pkgdir(ClimaLand)
global site_ID = "US-MOz"

for float_type in (Float32, Float64)
    # Make these global so we can use them in other ozark files
    global FT = float_type
    global toml_dict = LP.create_toml_dict(FT)
    global earth_param_set = LP.LandParameters(toml_dict)

    # Create model, set initial conditions, and setup most simulation parameters
    # Convert site_ID string to a Val so we can dispatch on it
    site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)

    # Get the default values for this site's domain, location, and parameters
    # Use finer dz_top for conservation test
    dz_top = FT(0.025)
    (; dz_tuple, nelements, zmin, zmax) =
        FluxnetSimulations.get_domain_info(FT, Val(site_ID_val); dz_top)
    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(site_ID_val))
    (; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

    (;
        soil_ν,
        soil_K_sat,
        soil_S_s,
        soil_vg_n,
        soil_vg_α,
        θ_r,
        ν_ss_quartz,
        ν_ss_om,
        ν_ss_gravel,
        z_0m_soil,
        z_0b_soil,
        soil_ϵ,
        soil_α_PAR,
        soil_α_NIR,
        Ω,
        χl,
        G_Function,
        α_PAR_leaf,
        λ_γ_PAR,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        ϵ_canopy,
        ac_canopy,
        g1,
        Drel,
        g0,
        Vcmax25,
        SAI,
        f_root_to_shoot,
        K_sat_plant,
        ψ63,
        Weibull_param,
        a,
        conductivity_model,
        retention_model,
        plant_ν,
        plant_S_s,
        rooting_depth,
        h_canopy,
    ) = FluxnetSimulations.get_parameters(FT, Val(site_ID_val))

    # TIME STEPPING:
    t0 = Float64(120 * 3600 * 24)# start mid year
    dt = Float64(150)
    # Use smaller `tf` for Float32 simulation
    tf = (FT == Float64) ? t0 + 3600 * 24 * 10 : t0 + 2 * dt
    # Select conv. condition based on float type due to different precision
    err = (FT == Float64) ? 1e-8 : 1e-4
    norm_condition = CTS.MaximumError(err)
    conv_checker = CTS.ConvergenceChecker(; norm_condition)
    max_iterations = 20
    timestepper = CTS.ARS111()
    ode_algo = CTS.IMEXAlgorithm(
        timestepper,
        CTS.NewtonsMethod(
            max_iters = max_iterations,
            update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
            convergence_checker = conv_checker,
        ),
    )

    # Set up the domain for the simulation
    land_domain = Column(;
        zlim = (zmin, zmax),
        nelements = nelements,
        dz_tuple = dz_tuple,
        longlat = (long, lat),
    )
    canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

    prognostic_land_components = (:canopy, :soil, :soilco2)

    # Get the atmospheric and radiation forcing data
    (start_date, stop_date) =
        FluxnetSimulations.get_data_dates(site_ID, time_offset)
    (; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
        site_ID,
        lat,
        long,
        time_offset,
        atmos_h,
        start_date,
        toml_dict,
        FT,
    )

    # Soil model
    soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
        PAR_albedo = soil_α_PAR,
        NIR_albedo = soil_α_NIR,
    )
    retention_parameters = (;
        ν = soil_ν,
        θ_r,
        K_sat = soil_K_sat,
        hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
    )
    composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)
    soil_forcing = (; atmos, radiation)
    soil = Soil.EnergyHydrology{FT}(
        land_domain,
        soil_forcing,
        toml_dict;
        prognostic_land_components,
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
        albedo = soil_albedo,
        retention_parameters,
        composition_parameters,
        S_s = soil_S_s,
        z_0m = z_0m_soil,
        z_0b = z_0b_soil,
        emissivity = soil_ϵ,
    )

    # Soil microbes model
    co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
    drivers = Soil.Biogeochemistry.SoilDrivers(co2_prognostic_soil, atmos)
    # The litter of the canopy carbon pools makes SOC prognostic
    soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(
        land_domain,
        drivers,
        toml_dict;
        sources = (
            Soil.Biogeochemistry.MicrobeProduction{FT}(),
            ClimaLand.LitterInputMinusMicrobialDecomposition{FT}(),
        ),
    )

    # Canopy model
    # Set up radiative transfer
    radiation_parameters =
        (; Ω, G_Function, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf)
    radiative_transfer = Canopy.TwoStreamModel{FT}(
        canopy_domain,
        toml_dict;
        radiation_parameters,
        ϵ_canopy,
    )

    # Set up conductance
    conductance =
        Canopy.MedlynConductanceModel{FT}(canopy_domain, toml_dict; g1)

    # Set up photosynthesis
    photosynthesis_parameters = (; fractional_c3 = FT(1), Vcmax25)
    photosynthesis =
        FarquharModel{FT}(canopy_domain, toml_dict; photosynthesis_parameters)

    # Set up plant hydraulics
    # Read in LAI from MODIS data
    surface_space = land_domain.space.surface
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date + Second(t0),
        stop_date + Second(tf),
    )
    # Get the maximum LAI at this site over the first year of the simulation
    maxLAI = FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long)
    RAI = FT(maxLAI) * f_root_to_shoot # convert to float type of simulation

    hydraulics = Canopy.PlantHydraulicsModel{FT}(
        canopy_domain,
        toml_dict;
        ν = plant_ν,
        S_s = plant_S_s,
        conductivity_model,
        retention_model,
    )
    height = h_canopy
    lai_model = Canopy.PrescribedAreaIndexModel{FT}(;
        LAI,
        SAI,
        RAI,
        rooting_depth,
        height,
    )
    # Carbon pools around the prescribed LAI, for the carbon budget; they do not feed back on
    # the water and energy budgets
    biomass = Canopy.PrognosticCarbonModel{FT}(
        lai_model,
        toml_dict;
        tree_share = FT(1),
    )

    # Set up energy model
    energy = Canopy.BigLeafEnergyModel{FT}(toml_dict; ac_canopy)

    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)
    # Construct the canopy model using defaults for autotrophic respiration (from the pools)
    # and SIF models
    canopy = Canopy.CanopyModel{FT}(
        canopy_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components,
        radiative_transfer,
        photosynthesis,
        conductance,
        hydraulics,
        energy,
        biomass,
    )

    # Integrated plant hydraulics and soil model
    land = SoilCanopyModel{FT}(soilco2, soil, canopy)
    exp_tendency! = make_exp_tendency(land)
    imp_tendency! = make_compute_imp_tendency(land)
    jacobian! = make_compute_jacobian(land)
    set_initial_cache! = make_set_initial_cache(land)
    Y, p, cds = initialize(land)
    jac_kwargs =
        (; jac_prototype = ClimaLand.initialize_jacobian(Y), Wfact = jacobian!)
    CL_cache_imp! = ClimaLand.make_update_implicit_cache(land)

    set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
        site_ID,
        start_date,
        time_offset,
        land,
    )
    set_ic!(Y, p, t0, land)
    # Nonzero carbon stocks (kg C m^-2), so that the pools shed litter to the soil
    Y.canopy.biomass.C_sugar .= FT(0.05)
    Y.canopy.biomass.C_leaf .= FT(0.2)
    Y.canopy.biomass.C_stem .= FT(5)
    Y.canopy.biomass.C_root .= FT(0.5)
    set_initial_cache!(p, Y, t0)


    saveat = dt
    n_saves = length(t0:saveat:tf)
    sv = (;
        t = Array{Float64}(undef, n_saves),
        saveval = Array{NamedTuple}(undef, n_saves),
    )
    saving_affect! = ClimaLand.SavingAffect(sv, 0)
    saving_initialize = (_, _, _, x) -> saving_affect!(x)
    saving_cb = ClimaLand.IntervalBasedCallback(
        saveat,
        t0,
        dt,
        saving_affect!;
        callback_start = t0,
        initialize = saving_initialize,
    )

    drivers = ClimaLand.get_drivers(land)
    updatefunc = ClimaLand.make_update_drivers(drivers)
    driver_cb = ClimaLand.DriverUpdateCallback(updatefunc, dt, t0)
    prob = ClimaTimeSteppers.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = ClimaTimeSteppers.ODEFunction(
                imp_tendency!;
                jac_kwargs...,
            ),
            dss! = ClimaLand.dss!,
            cache_imp! = (Y, p, t) -> CL_cache_imp!(p, Y, t),
        ),
        Y,
        (t0, tf),
        p,
    )
    cb = ClimaTimeSteppers.CallbackSet(driver_cb, saving_cb)
    sol = ClimaTimeSteppers.solve(
        prob,
        ode_algo;
        dt = dt,
        callback = cb,
        adaptive = false,
        saveat = collect(t0:saveat:tf),
    )

    # Check that simulation still has correct float type
    @assert eltype(sol.u[end].soil) == FT
    @assert eltype(sol.u[end].soilco2) == FT
    @assert eltype(sol.u[end].canopy) == FT

    # Plotting for Float64 simulation
    if FT == Float64
        # Check that the driver variables stored in `p` are
        # updated correctly - just check one from radiation,
        # and one from atmos.
        cache_cosθs =
            [parent(sv.saveval[k].drivers.cosθs)[1] for k in 1:length(sv.t)]
        cache_Tair =
            [parent(sv.saveval[k].drivers.T)[1] for k in 1:length(sv.t)]
        @assert mean(
            abs.(
                first.(radiation.cosθs.(sv.t, radiation.start_date)) .-
                cache_cosθs,
            ),
        ) < eps(FT)
        T_mutable = Vector{FT}(undef, 1)
        atmos_T = map(sv.t) do time
            ClimaLand.evaluate!(T_mutable, atmos.T, time)
            return T_mutable[]
        end |> collect

        @assert mean(abs.(atmos_T .- cache_Tair)) < eps(FT)

        daily = sol.t[2:end] ./ 3600 ./ 24
        savedir = generate_output_path(
            "experiments/integrated/performance/conservation",
        )

        # The aux state in `sv` is an index off from the solution..
        # Aux state at index 2 corresponds to the solution at index 1
        # This is a bug that we will need to fix at some point
        # Internally during the integration, however, aux and state
        # are updated consistently.

        ##  Soil water balance ##

        # Evaporation
        E = [
            parent(
                sv.saveval[k].soil.turbulent_fluxes.vapor_flux_liq .+
                sv.saveval[k].soil.turbulent_fluxes.vapor_flux_ice,
            )[1] for k in 2:length(sol.t)
        ]
        # Root sink term: a positive root extraction is a sink term for soil; add minus sign
        root_sink =
            [sum(-1 .* sv.saveval[k].root_extraction) for k in 2:length(sol.t)]
        # Free Drainage BC, F_bot = -K_bot
        soil_bottom_flux =
            [-parent(sv.saveval[k].soil.K)[1] for k in 2:length(sol.t)]
        # Precip is not stored in the aux state, evaluated using sol.t
        p_mutable = Vector{FT}(undef, 1)

        precip =
            map(sol.t[k] for k in 1:(length(sol.t) - 1)) do time
                ClimaLand.evaluate!(p_mutable, atmos.liquid_precip, time)
                return p_mutable[]
            end |> collect

        # Water balance equation

        # d[∫(ϑ_l+ θ_i ρ_i/ρ_l)dz] = [-(F_sfc - F_bot) + ∫Sdz]dt = -ΔF dt + ∫Sdz dt
        # N.B. in ClimaCore, sum(field) -> integral
        rhs_soil = -(precip .+ E .- soil_bottom_flux) .+ root_sink

        ρ_liq = LP.ρ_cloud_liq(earth_param_set)
        ρ_ice = LP.ρ_cloud_ice(earth_param_set)
        net_soil_water_storage = [
            sum(sol.u[k].soil.ϑ_l .+ sol.u[k].soil.θ_i .* ρ_ice ./ ρ_liq)[1]
            for k in 1:length(sol.t)
        ]
        lhs_soil = net_soil_water_storage[2:end] .- net_soil_water_storage[1]
        soil_mass_change_actual = lhs_soil
        soil_mass_change_exp = cumsum(rhs_soil) .* dt

        ## Canopy water balance ##

        # Bottom flux from roots to plant
        root_flux =
            [sum(sv.saveval[k].root_extraction) for k in 2:length(sol.t)]

        # LAI (t)
        leaf_area_index = [
            parent(
                getproperty(sv.saveval[k].canopy.biomass.area_index, :leaf),
            )[1] for k in 2:length(sol.t)
        ]
        # Top boundary flux (transpiration)
        T = [
            parent(sv.saveval[k].canopy.turbulent_fluxes.vapor_flux)[1] for
            k in 2:length(sol.t)
        ]

        # Water balance equation
        # d[θ_leaf h_leaf] =  - [T - Root Flux]/LAI
        rhs_canopy = @. -T / leaf_area_index + root_flux / leaf_area_index
        net_plant_water_storage = [
            sum(parent(sol.u[k].canopy.hydraulics.ϑ_l) .* h_canopy) for
            k in 1:length(sol.t)
        ]
        lhs_canopy =
            net_plant_water_storage[2:end] .- net_plant_water_storage[1]

        canopy_mass_change_actual = lhs_canopy
        canopy_mass_change_exp = cumsum(rhs_canopy) .* dt

        fig = Figure(size = (1500, 400))
        ax = Axis(
            fig[1, 1],
            xlabel = "Days",
            ylabel = "Fractional Error",
            yscale = log10,
        )

        lines!(
            ax,
            daily,
            eps(FT) .+ abs.(
                (soil_mass_change_actual - soil_mass_change_exp) ./
                soil_mass_change_exp,
            ),
            label = "Soil Water Balance",
        )
        lines!(
            ax,
            daily,
            eps(FT) .+ abs.(
                (canopy_mass_change_actual - canopy_mass_change_exp) ./
                canopy_mass_change_exp,
            ),
            label = "Canopy Water Balance",
        )
        axislegend(position = :lb)
        CairoMakie.save(joinpath(savedir, "water_conservation.png"), fig)



        ##  Soil energy balance ##
        # Energy of liquid water infiltrating soil is ignored in our model.

        # Turbulent fluxes
        LHF = [
            parent(sv.saveval[k].soil.turbulent_fluxes.lhf)[1] for
            k in 2:length(sol.t)
        ]
        SHF = [
            parent(sv.saveval[k].soil.turbulent_fluxes.shf)[1] for
            k in 2:length(sol.t)
        ]
        # Radiation is computed in LW and SW components
        # with positive numbers indicating the soil gaining energy.
        soil_Rn =
            -1 .* [parent(sv.saveval[k].soil.R_n)[1] for k in 2:length(sol.t)]
        # Root sink term: a positive root extraction is a sink term for soil; add minus sign
        root_sink_energy = [
            sum(-1 .* sv.saveval[k].root_energy_extraction) for
            k in 2:length(sol.t)
        ]
        # Bottom energy BC
        soil_bottom_flux = FT(0)

        # Energy balance equation

        # d[∫Idz] = [-(F_sfc - F_bot) + ∫Sdz]dt = -ΔF dt + ∫Sdz dt
        # N.B. in ClimaCore, sum(field) -> integral
        rhs_soil_energy =
            -(LHF .+ SHF .+ soil_Rn .- soil_bottom_flux) .+ root_sink_energy

        net_soil_energy_storage =
            [sum(sol.u[k].soil.ρe_int)[1] for k in 1:length(sol.t)]
        lhs_soil_energy =
            net_soil_energy_storage[2:end] .- net_soil_energy_storage[1]
        soil_energy_change_actual = lhs_soil_energy
        soil_energy_change_exp = cumsum(rhs_soil_energy) .* dt

        fig = Figure(size = (1500, 400))
        ax = Axis(
            fig[1, 1],
            xlabel = "Days",
            ylabel = "Fractional Error",
            yscale = log10,
        )

        lines!(
            ax,
            daily,
            eps(FT) .+ abs.(
                (soil_energy_change_actual - soil_energy_change_exp) ./
                soil_energy_change_exp,
            ),
            label = "Soil Energy Balance",
        )

        axislegend(position = :lb)
        CairoMakie.save(joinpath(savedir, "energy_conservation.png"), fig)

        ## Carbon balance ##
        # Stocks in kg C m^-2: the vegetation pools, and the column integrals of SOC and of
        # soil CO2 (kg C per m^3 of soil, gas and dissolved).
        M_C = canopy.biomass.parameters.M_C
        surface_value(x) = parent(x)[1]
        # Per unit ground area; `sum` over this column carries a horizontal measure
        ∫dz = ClimaCore.Fields.zeros(surface_space)
        function column_integral(field)
            ClimaCore.Operators.column_integral_definite!(∫dz, field)
            return surface_value(∫dz)
        end
        pools = (:C_sugar, :C_leaf, :C_stem, :C_root)
        veg_carbon = [
            sum(
                pool ->
                    surface_value(getproperty(Y_k.canopy.biomass, pool)),
                pools,
            ) for Y_k in sol.u
        ]
        soc = [column_integral(Y_k.soilco2.SOC) for Y_k in sol.u]
        soil_co2 = [column_integral(Y_k.soilco2.CO2) for Y_k in sol.u]
        land_carbon = veg_carbon .+ soc .+ soil_co2

        # Fluxes in kg C m^-2 s^-1 of each step, from the aux state saved after it
        function carbon_fluxes(aux)
            carbon = aux.canopy.biomass.carbon
            return (;
                GPP = M_C *
                      surface_value(Canopy.get_GPP(aux, canopy.photosynthesis)),
                Ra = surface_value(carbon.Ra),
                litter = surface_value(carbon.L_leaf) +
                         surface_value(carbon.L_stem) +
                         surface_value(carbon.L_root),
                soil_litter = column_integral(aux.soil_litter_input),
                Sm = column_integral(aux.soilco2.Sm),
                F_top = surface_value(aux.soilco2.top_bc),
                F_bot = surface_value(aux.soilco2.bottom_bc),
            )
        end

        # Each budget: its stock, net flux into it and gross flux through it
        function carbon_budgets(f)
            return (
                (
                    "Vegetation",
                    veg_carbon,
                    f.GPP - f.Ra - f.litter,
                    f.GPP + f.Ra + f.litter,
                ),
                (
                    "Soil organic carbon",
                    soc,
                    f.soil_litter - f.Sm,
                    f.soil_litter + f.Sm,
                ),
                (
                    "Soil CO2",
                    soil_co2,
                    f.Sm - f.F_top + f.F_bot,
                    f.Sm + abs(f.F_top) + abs(f.F_bot),
                ),
                (
                    "Land",
                    land_carbon,
                    f.GPP - f.Ra - f.F_top + f.F_bot,
                    f.GPP + f.Ra + abs(f.F_top) + abs(f.F_bot),
                ),
            )
        end

        # Error in the change of each stock, relative to the carbon that went through it
        function carbon_errors(fluxes)
            budgets = carbon_budgets.(fluxes)
            return map(1:4) do i
                name, stock = budgets[1][i][1], budgets[1][i][2]
                net = [b[i][3] for b in budgets]
                gross = [b[i][4] for b in budgets]
                change = stock[2:end] .- stock[1]
                error =
                    abs.(change .- cumsum(net) .* dt) ./ (cumsum(gross) .* dt)
                (name, error)
            end
        end

        fluxes = [carbon_fluxes(sv.saveval[k]) for k in 2:length(sol.t)]
        carbon_error = carbon_errors(fluxes)

        fig = Figure(size = (1500, 400))
        ax = Axis(
            fig[1, 1],
            xlabel = "Days",
            ylabel = "Error relative to cumulative gross flux",
            yscale = log10,
        )
        for (name, error) in carbon_error
            lines!(ax, daily, eps(FT) .+ error, label = name)
        end
        axislegend(position = :lb)
        CairoMakie.save(joinpath(savedir, "carbon_conservation.png"), fig)
        for (name, error) in carbon_error
            @info "Carbon balance: $name" maximum(error)
            @assert maximum(error) < 1e-6 "Carbon is not conserved: $name"
        end
    end
end
