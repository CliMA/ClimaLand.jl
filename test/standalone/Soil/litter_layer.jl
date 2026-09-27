using Test
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Dates
using Statistics
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaTimeSteppers as CTS

# A US-Var-like soil column under idealized diurnal forcing
function litter_test_model(FT, toml_dict, surface_layer; start_date)
    SW_d(t) = 800.0 * max(0.0, cos(2π * (float(t) / 86400 - 0.5)))
    T_atmos(t) = 293.0 + 6.0 * cos(2π * (float(t) / 86400 - 0.625))
    LW_d(t) = 0.8 * 5.67e-8 * T_atmos(t)^4
    atmos, radiation = ClimaLand.prescribed_analytic_forcing(
        FT;
        toml_dict,
        start_date,
        SW_d,
        LW_d,
        T_atmos,
        u_atmos = (t) -> 2.0,
        q_atmos = (t) -> 0.007,
        h_atmos = FT(2),
    )
    domain = Column(;
        zlim = (FT(-0.5), FT(0)),
        nelements = 24,
        dz_tuple = FT.((0.05, 0.02)),
    )
    params = Soil.EnergyHydrologyParameters(
        toml_dict;
        ν = FT(0.5),
        ν_ss_om = FT(0.02),
        ν_ss_quartz = FT(0.3),
        ν_ss_gravel = FT(0.0),
        hydrology_cm = vanGenuchten{FT}(; α = FT(2.0), n = FT(1.6)),
        K_sat = FT(0.45 / 3600 / 100),
        S_s = FT(1e-3),
        θ_r = FT(0.0),
        albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
            PAR_albedo = FT(0.35),
            NIR_albedo = FT(0.35),
        ),
        emissivity = FT(0.98),
        z_0m = FT(0.01),
        z_0b = FT(0.01),
    )
    boundary_conditions = (;
        top = Soil.AtmosDrivenFluxBC(atmos, radiation),
        bottom = WaterHeatBC(;
            water = WaterFluxBC((p, t) -> 0.0),
            heat = HeatFluxBC((p, t) -> 0.0),
        ),
    )
    return Soil.EnergyHydrology{FT}(;
        parameters = params,
        domain,
        boundary_conditions,
        sources = (),
        surface_layer,
    )
end

# Standalone soil has no canopy; seed the litter with a grassland-like PAI
const PAI_TEST = 1.5
function litter_test_ic!(Y, p, t0, model)
    FT = eltype(Y.soil.ϑ_l)
    params = model.parameters
    Y.soil.ϑ_l .= FT(0.2)
    Y.soil.θ_i .= FT(0)
    ρc_s = @. Soil.volumetric_heat_capacity(
        Y.soil.ϑ_l,
        Y.soil.θ_i,
        params.ρc_ds,
        params.earth_param_set,
    )
    Y.soil.ρe_int .= Soil.volumetric_internal_energy.(
        FT(0),
        ρc_s,
        FT(293),
        params.earth_param_set,
    )
    Soil.initialize_litter_temperature!(Y, model)
    Soil.initialize_litter_area_index!(Y, model, FT(PAI_TEST))
end

function run_litter_column(model, start_date, ndays, dt)
    stop_date = start_date + Day(ndays)
    saveat = Second(900)
    saving_cb = ClimaLand.NonInterpSavingCallback(start_date, stop_date, saveat)
    sv = saving_cb.affect!.saved_values
    simulation = LandSimulation(
        start_date,
        stop_date,
        dt,
        model;
        set_ic! = litter_test_ic!,
        updateat = Second(dt),
        solver_kwargs = (; saveat),
        user_callbacks = (saving_cb,),
        diagnostics = (),
    )
    sol = solve!(simulation)
    return sol, sv, simulation
end

for FT in (Float32, Float64)
    @testset "SlabLitter surface layer, FT = $FT" begin
        toml_dict = LP.create_toml_dict(FT)
        earth_param_set = LP.LandParameters(toml_dict)
        start_date = DateTime(2005, 7, 1)
        litter = SlabLitter{FT}(toml_dict, 450)
        @test litter.d_PAI == FT(toml_dict["litter_thickness_per_pai"])
        @test litter.c_vap == FT(toml_dict["litter_vapor_resistance_factor"])
        @test litter.Δt == FT(450)
        model = litter_test_model(FT, toml_dict, litter; start_date)
        @test model.surface_layer === litter
        @test :T_litter in ClimaLand.prognostic_vars(model)
        @test :PAI_mean in ClimaLand.prognostic_vars(model)
        @test :litter in ClimaLand.auxiliary_vars(model)
        @test :skin_solve in ClimaLand.auxiliary_vars(model)

        Y, p, coords = initialize(model)
        litter_test_ic!(Y, p, FT(0), model)
        set_initial_cache! = make_set_initial_cache(model)
        set_initial_cache!(p, Y, FT(0))
        T_top = ClimaLand.Domains.top_center_to_surface(p.soil.T)
        @test Y.soil.T_litter == T_top
        @test all(parent(Y.soil.PAI_mean) .== FT(PAI_TEST))
        d_l = Soil.litter_thickness(litter, Y)
        @test all(parent(@. d_l + 0) .≈ litter.d_PAI * FT(PAI_TEST))
        _D_vap = FT(LP.D_vapor(earth_param_set))
        r_vap_l = Soil.litter_vapor_resistance(litter, Y, _D_vap)
        @test all(
            parent(@. r_vap_l + 0) .≈
            litter.c_vap * (litter.d_PAI * FT(PAI_TEST) - litter.d_min) /
            _D_vap,
        )

        # The skin closes the energy balance against the litter node
        r_top = Soil.litter_half_resistance(litter, Y)
        F_atm = @. p.soil.R_n +
           p.soil.turbulent_fluxes.lhf +
           p.soil.turbulent_fluxes.shf
        F_scale = @. abs(p.soil.R_n) +
           abs(p.soil.turbulent_fluxes.lhf) +
           abs(p.soil.turbulent_fluxes.shf)
        G_skin = @. (Y.soil.T_litter - p.soil.turbulent_fluxes.T_sfc) / r_top
        @test all(
            abs.(parent(F_atm) .- parent(G_skin)) .<
            FT(0.01) .* parent(F_scale),
        )
        # The soil BC uses the litter temperature after the implicit substep;
        # the litter cache holds the atmospheric flux and a positive sensitivity
        r_bot = Soil.litter_soil_resistance(litter, model, Y, p)
        T_l = p.soil.litter.T
        @test all(
            parent(p.soil.top_bc.heat) .≈ parent(@. (T_top - T_l) / r_bot),
        )
        @test all(parent(T_l) .!= parent(Y.soil.T_litter))
        @test all(parent(p.soil.litter.F_atm) .≈ parent(F_atm))
        @test all(parent(p.soil.litter.∂F_atm∂T) .> 0)
        @test all(parent(p.soil.litter.T_n) .== parent(Y.soil.T_litter))
        # Total energy includes the litter
        total_energy = similar(p.soil.total_energy)
        ClimaLand.total_energy_per_area!(total_energy, model, Y, p, FT(0))
        soil_energy = similar(total_energy)
        ClimaCore.Operators.column_integral_definite!(
            soil_energy,
            Y.soil.ρe_int,
        )
        C_l = litter.ρc_l * litter.d_PAI * FT(PAI_TEST)
        _T_ref = FT(LP.T_0(earth_param_set))
        @test all(
            parent(total_energy) .≈
            parent(soil_energy) .+ C_l .* (parent(Y.soil.T_litter) .- _T_ref),
        )
        # Without a canopy the trailing PAI relaxes toward zero
        dY = similar(Y)
        dY .= 0
        exp_tendency! = make_compute_exp_tendency(model)
        exp_tendency!(dY, Y, p, FT(0))
        @test all(parent(dY.soil.PAI_mean) .≈ -FT(PAI_TEST) / litter.τ_PAI)
        @test all(parent(dY.soil.T_litter) .== 0)

        # Energy conservation over a two-day run: the change in soil plus
        # litter energy equals the accumulated boundary flux (including the
        # energy of the litter mass lost as the trailing PAI decays)
        sol, sv, sim = run_litter_column(model, start_date, 2, FT(450))
        p_end = sim._integrator.p
        Y_end = sim._integrator.u
        E_end = similar(p_end.soil.total_energy)
        ClimaLand.total_energy_per_area!(E_end, model, Y_end, p_end, FT(0))
        # The simulation was initialized with the same state as `total_energy`
        E_start = Array(parent(total_energy))[1]
        ΔE = Array(parent(E_end))[1] - E_start
        ∫F = Array(parent(Y_end.soil.∫F_e_dt))[1]
        @info "SlabLitter energy budget" FT ΔE ∫F ΔE - ∫F
        @test abs(ΔE - ∫F) < 1e-3 * abs(ΔE) + FT(10)
        @test all(isfinite, parent(Y_end.soil.T_litter))
        @test all(parent(Y_end.soil.PAI_mean) .< FT(PAI_TEST))

        # The litter damps the diurnal cycle of the top soil layer relative to
        # the skin scheme, and the skin swings less than without heat capacity
        nolitter = litter_test_model(FT, toml_dict, NoLitter{FT}(); start_date)
        sol0, sv0, _ = run_litter_column(nolitter, start_date, 2, FT(450))
        lastday(sv) =
            [Array(parent(s.soil.T))[end] for s in sv.saveval[(end - 96):end]]
        range_top_litter = maximum(lastday(sv)) - minimum(lastday(sv))
        range_top_nolitter = maximum(lastday(sv0)) - minimum(lastday(sv0))
        @test range_top_litter < 0.8 * range_top_nolitter
        skin(sv) = [
            Array(parent(s.soil.turbulent_fluxes.T_sfc))[1] for
            s in sv.saveval[(end - 96):end]
        ]
        range_skin_litter = maximum(skin(sv)) - minimum(skin(sv))
        range_skin_nolitter = maximum(skin(sv0)) - minimum(skin(sv0))
        @test range_skin_litter < 1.3 * range_skin_nolitter

        # Time step insensitivity of the implicit litter solve
        litter_s = SlabLitter{FT}(toml_dict, 90)
        model_s = litter_test_model(FT, toml_dict, litter_s; start_date)
        sol_s, sv_s, sim_s = run_litter_column(model_s, start_date, 2, FT(90))
        T_top_450 = Array(parent(sim._integrator.p.soil.T))[end]
        T_top_90 = Array(parent(sim_s._integrator.p.soil.T))[end]
        @test abs(T_top_450 - T_top_90) < FT(0.5)
        T_l_450 = Array(parent(sim._integrator.u.soil.T_litter))[1]
        T_l_90 = Array(parent(sim_s._integrator.u.soil.T_litter))[1]
        @test abs(T_l_450 - T_l_90) < FT(1)

        # The stored time step must match the simulation time step
        @test_throws ArgumentError run_litter_column(
            model,
            start_date,
            1,
            FT(900),
        )

        # The scheme is a backward Euler substep of the full step and lags the
        # soil by one Newton iteration: ARS111 with at least two iterations
        for timestepper in (
            CTS.IMEXAlgorithm(CTS.ARS343(), CTS.NewtonsMethod(max_iters = 3)),
            CTS.IMEXAlgorithm(CTS.ARS111(), CTS.NewtonsMethod(max_iters = 1)),
        )
            @test_throws ArgumentError LandSimulation(
                start_date,
                start_date + Day(1),
                FT(450),
                model;
                set_ic! = litter_test_ic!,
                updateat = Second(450),
                timestepper,
                diagnostics = (),
            )
        end

        # A SlabLitter requires an atmosphere-driven boundary condition, but
        # the atmosphere may be prescribed or coupled
        bc_heat = (;
            top = WaterHeatBC(;
                water = WaterFluxBC((p, t) -> 0.0),
                heat = HeatFluxBC((p, t) -> 0.0),
            ),
            bottom = model.boundary_conditions.bottom,
        )
        @test_throws ArgumentError Soil.EnergyHydrology{FT}(;
            parameters = model.parameters,
            domain = model.domain,
            boundary_conditions = bc_heat,
            sources = (),
            surface_layer = litter,
        )
        coupled_bc = (;
            top = Soil.AtmosDrivenFluxBC(
                ClimaLand.CoupledAtmosphere{FT, FT}(FT(2), FT(1)),
                ClimaLand.CoupledRadiativeFluxes{FT}(),
            ),
            bottom = model.boundary_conditions.bottom,
        )
        coupled_model = Soil.EnergyHydrology{FT}(;
            parameters = model.parameters,
            domain = model.domain,
            boundary_conditions = coupled_bc,
            sources = (),
            surface_layer = litter,
        )
        @test coupled_model.surface_layer === litter
        @test :T_litter in ClimaLand.prognostic_vars(coupled_model)
    end
end
