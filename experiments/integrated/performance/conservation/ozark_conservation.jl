import SciMLBase
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
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaParams
import ClimaUtilities.OutputPathGenerator: generate_output_path
using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations

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
        n_stem,
        n_leaf,
        h_leaf,
        h_stem,
        h_canopy,
        z_0m,
        z_0b,
    ) = FluxnetSimulations.get_parameters(FT, Val(site_ID_val))

    compartment_midpoints =
        n_stem > 0 ? [h_stem / 2, h_stem + h_leaf / 2] : [h_leaf / 2]
    compartment_surfaces =
        n_stem > 0 ? [zmax, h_stem, h_canopy] : [zmax, h_leaf]

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
    soil_organic_carbon =
        ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))
    co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
    drivers = Soil.Biogeochemistry.SoilDrivers(
        co2_prognostic_soil,
        soil_organic_carbon,
        atmos,
    )
    soilco2 =
        Soil.Biogeochemistry.SoilCO2Model{FT}(land_domain, drivers, toml_dict)

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
    photosynthesis_parameters = (; is_c3 = FT(1), Vcmax25)
    photosynthesis =
        FarquharModel{FT}(canopy_domain, toml_dict; photosynthesis_parameters)

    # Set up optimal LAI model
    lai_model =
        Canopy.OptimalLAIModel{FT}(Canopy.OptimalLAIParameters{FT}(toml_dict))

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
        n_stem,
        n_leaf,
        h_stem,
        h_leaf,
        ν = plant_ν,
        S_s = plant_S_s,
        conductivity_model,
        retention_model,
    )
    height = h_stem + h_leaf
    biomass = Canopy.PrescribedBiomassModel{FT}(;
        LAI,
        SAI,
        RAI,
        rooting_depth,
        height,
    )

    # Set up energy model
    energy = Canopy.BigLeafEnergyModel{FT}(toml_dict; ac_canopy)

    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)
    # Construct the canopy model using defaults for autotrophic respiration and SIF models
    canopy = Canopy.CanopyModel{FT}(
        canopy_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        z_0m,
        z_0b,
        prognostic_land_components,
        radiative_transfer,
        photosynthesis,
        conductance,
        lai_model,
        hydraulics,
        energy,
        biomass,
    )

    # Integrated plant hydraulics and soil model
    land = SoilCanopyModel{FT}(soilco2, soil, canopy)
    exp_tendency! = make_exp_tendency(land)
    imp_tendency! = make_imp_tendency(land)
    jacobian! = make_jacobian(land)
    set_initial_cache! = make_set_initial_cache(land)
    Y, p, cds = initialize(land)
    jac_kwargs = (;
        jac_prototype = ClimaLand.FieldMatrixWithSolver(Y),
        Wfact = jacobian!,
    )

    set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
        site_ID,
        start_date,
        time_offset,
        land,
    )
    set_ic!(Y, p, t0, land)
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
    cb = SciMLBase.CallbackSet(driver_cb, saving_cb)
    sol = SciMLBase.solve(
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
                cos.(radiation.θs.(sv.t, radiation.start_date)) .- cache_cosθs,
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

        # Bottom flux from roots to stem
        root_flux =
            [sum(sv.saveval[k].root_extraction) for k in 2:length(sol.t)]

        # Stem leaf flux
        stem_leaf_flux = [
            parent(sv.saveval[k].canopy.hydraulics.fa)[1] for
            k in 2:length(sol.t)
        ]
        # LAI (t)
        leaf_area_index = [
            parent(
                getproperty(sv.saveval[k].canopy.biomass.area_index, :leaf),
            )[1] for k in 2:length(sol.t)
        ]
        # Top boundary flux (transpiration)
        T = [
            parent(sv.saveval[k].canopy.turbulent_fluxes.transpiration)[1]
            for k in 2:length(sol.t)
        ]

        # Water balance equation
        # d[θ_leaf h_leaf+ θ_stem h_stem] = -[F_sl - Root Flux]/SAI - [T - F_sl]/LAI
        rhs_canopy = @. -T / leaf_area_index +
           root_flux / SAI +
           stem_leaf_flux * (1 / leaf_area_index - 1 / SAI)

        net_plant_water_storage = [
            sum(parent(sol.u[k].canopy.hydraulics.ϑ_l) .* [h_stem, h_leaf])
            for k in 1:length(sol.t)
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
            eps(FT) .+
            abs.(
                (soil_mass_change_actual - soil_mass_change_exp) ./
                soil_mass_change_exp,
            ),
            label = "Soil Water Balance",
        )
        lines!(
            ax,
            daily,
            eps(FT) .+
            abs.(
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
            eps(FT) .+
            abs.(
                (soil_energy_change_actual - soil_energy_change_exp) ./
                soil_energy_change_exp,
            ),
            label = "Soil Energy Balance",
        )

        axislegend(position = :lb)
        CairoMakie.save(joinpath(savedir, "energy_conservation.png"), fig)
    end
end
