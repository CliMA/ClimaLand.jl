using Test
using Dates
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using Insolation
import ClimaParams as CP
using ClimaLand
using ClimaLand.Canopy
using ClimaLand.Domains: Point

import ClimaLand
import ClimaLand.Parameters as LP

@testset "Canopy interception pointwise functions" begin
    for FT in (Float32, Float64)
        p_liq = FT(1e-4)
        f_wet_max = FT(0.05)
        # Wetted fraction
        @test Canopy.wetted_fraction(FT(0), FT(3), p_liq, f_wet_max) == 0
        @test Canopy.wetted_fraction(FT(3e-4), FT(3), p_liq, FT(1)) ≈ 1
        @test Canopy.wetted_fraction(FT(3e-4), FT(3), p_liq, f_wet_max) ≈
              f_wet_max
        @test Canopy.wetted_fraction(FT(1e-6), FT(0.01), p_liq, FT(1)) == 0
        W = FT(0.5e-4)
        @test Canopy.wetted_fraction(W, FT(2), p_liq, FT(1)) ≈
              (W / (2p_liq))^(FT(2) / 3)

        # Conductances
        u_star = FT(0.4)
        Cd = FT(0.01)
        LAI = FT(3)
        PAI = FT(4)
        r_st = FT(50)
        ρ = FT(1.2)
        q_c = FT(0.015)
        q_a = FT(0.010)
        g_leaf = Cd * u_star * LAI
        g_dry = (1 / r_st) * g_leaf / (g_leaf + 1 / r_st)
        # Without interception, the dry-canopy conductance is recovered
        (g_wet, g_tr) = Canopy.canopy_vapor_conductances(
            u_star,
            Cd,
            LAI,
            PAI,
            r_st,
            FT(0),
            FT(0),
            q_c,
            q_a,
            ρ,
            false,
        )
        @test g_wet == 0
        @test g_tr == g_dry
        # A wet canopy with ample water evaporates at close to the boundary
        # layer conductance of the wetted area
        f_wet = FT(0.5)
        (g_wet, g_tr) = Canopy.canopy_vapor_conductances(
            u_star,
            Cd,
            LAI,
            PAI,
            r_st,
            f_wet,
            FT(1e3),
            q_c,
            q_a,
            ρ,
            true,
        )
        @test g_wet ≈ f_wet * Cd * u_star * PAI rtol = 1e-3
        @test g_tr ≈ (1 - f_wet) * g_dry
        # Wet-canopy evaporation is limited by the available water
        E_max = FT(1e-6)
        (g_wet, g_tr) = Canopy.canopy_vapor_conductances(
            u_star,
            Cd,
            LAI,
            PAI,
            r_st,
            f_wet,
            E_max,
            q_c,
            q_a,
            ρ,
            true,
        )
        @test ρ * g_wet * (q_c - q_a) < E_max
        # Condensation goes to the canopy store only
        (g_wet, g_tr) = Canopy.canopy_vapor_conductances(
            u_star,
            Cd,
            LAI,
            PAI,
            r_st,
            FT(0),
            FT(0),
            q_a,
            q_c,
            ρ,
            true,
        )
        @test g_tr == 0
        @test g_wet ≈ Cd * u_star * PAI
    end
end

@testset "Canopy interception water budget" begin
    FT = Float64
    longlat = (FT(-118.14), FT(34.15))
    domain = Point(; z_sfc = FT(0), longlat)
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    t0 = 0.0
    start_date = DateTime(2005)
    Δt = FT(180.0)

    SW_d = TimeVaryingInput((t) -> eltype(t)(20.0))
    LW_d = TimeVaryingInput((t) -> eltype(t)(300.0))
    radiation = PrescribedRadiativeFluxes(
        FT,
        SW_d,
        LW_d,
        start_date;
        cosθs = (t, s) -> default_cos_zenith_angle(
            t,
            s;
            insol_params = earth_param_set.insol_params,
            latitude = FT(40.0),
            longitude = FT(-120.0),
        ),
        toml_dict = toml_dict,
    )
    P_rain = FT(-2e-6) # m/s, 7.2 mm/h
    rain = TimeVaryingInput((t) -> eltype(t)(P_rain))
    snow = TimeVaryingInput((t) -> eltype(t)(0))
    T_atmos = TimeVaryingInput((t) -> eltype(t)(290.0))
    u_atmos = TimeVaryingInput((t) -> eltype(t)(2.0))
    q_atmos = TimeVaryingInput((t) -> eltype(t)(0.009))
    P_atmos = TimeVaryingInput((t) -> eltype(t)(101325))
    atmos = ClimaLand.PrescribedAtmosphere(
        rain,
        snow,
        T_atmos,
        u_atmos,
        q_atmos,
        P_atmos,
        start_date,
        FT(10),
        toml_dict,
    )
    ground = PrescribedGroundConditions{FT}()
    LAI_value = FT(3)
    SAI = FT(1)
    LAI = TimeVaryingInput(t -> LAI_value)
    biomass = Canopy.PrescribedBiomassModel{FT}(;
        LAI,
        SAI,
        RAI = FT(1),
        rooting_depth = FT(0.5),
        height = FT(2),
    )
    interception = Canopy.CLM5Interception{FT}(toml_dict, Δt; f_wet_max = 1)
    canopy = Canopy.CanopyModel{FT}(
        domain,
        (; radiation, atmos, ground),
        LAI,
        toml_dict;
        biomass,
        interception,
    )
    @test ClimaLand.prognostic_vars(canopy.interception) == (:W,)
    Y, p, cds = initialize(canopy)
    ϑ0 = canopy.hydraulics.parameters.ν / 2
    Y.canopy.hydraulics.ϑ_l .= ϑ0
    Y.canopy.energy.T .= FT(289)
    Y.canopy.interception.W .= 0
    set_initial_cache! = make_set_initial_cache(canopy)
    exp_tendency! = make_exp_tendency(canopy)
    dY = similar(Y)
    PAI = LAI_value + SAI
    W_max = interception.p_liq * PAI

    # With a dry canopy at the start, the interception is α tanh(PAI) P
    set_initial_cache!(p, Y, t0)
    I = interception.α_liq * tanh(PAI) * (-P_rain)
    @test all(parent(p.canopy.interception.intercepted_liq) .≈ I)
    @test all(parent(p.canopy.interception.drip) .== 0)
    # The wetted fraction is evaluated after interception in this step
    @test all(
        parent(p.canopy.interception.f_wet) .≈
        Canopy.wetted_fraction(I * Δt, PAI, interception.p_liq, FT(1)),
    )
    # Precipitation is partitioned into interception and throughfall
    @test all(
        parent(
            p.canopy.interception.throughfall_liq .-
            p.canopy.interception.intercepted_liq .+
            p.canopy.interception.drip,
        ) .≈ P_rain,
    )

    # Step the store forward: it stays within [0, W_max + I Δt] and fills up
    for step in 1:200
        set_initial_cache!(p, Y, t0 + step * Δt)
        exp_tendency!(dY, Y, p, t0 + step * Δt)
        # The canopy vapor flux is partitioned into transpiration and
        # evaporation of intercepted water
        E_wet =
            p.canopy.turbulent_fluxes.vapor_flux .-
            p.canopy.interception.transpiration
        @test all(
            parent(dY.canopy.interception.W) .≈ parent(
                p.canopy.interception.intercepted_liq .-
                p.canopy.interception.drip .- E_wet,
            ),
        )
        Y.canopy.interception.W .+= Δt .* dY.canopy.interception.W
        @test all(parent(Y.canopy.interception.W) .>= -sqrt(eps(FT)) * W_max)
        @test all(parent(Y.canopy.interception.W) .<= W_max + I * Δt)
    end
    @test all(parent(Y.canopy.interception.W) .> FT(0.5) * W_max)
    # The canopy total water includes the intercepted water
    total_water = ClimaCore.Fields.zeros(domain.space.surface)
    ClimaLand.total_liq_water_vol_per_area!(total_water, canopy, Y, p, t0)
    total_water_hydraulics = ClimaCore.Fields.zeros(domain.space.surface)
    ClimaLand.total_liq_water_vol_per_area!(
        total_water_hydraulics,
        canopy.hydraulics,
        canopy,
        Y,
        p,
        t0,
    )
    @test all(
        parent(total_water) .≈
        parent(total_water_hydraulics .+ Y.canopy.interception.W),
    )

    # Without rain the store evaporates and does not become negative, apart
    # from round-off and the small difference between the air density used
    # in the evaporation limit and in the surface flux computation
    rain_off = ClimaLand.PrescribedAtmosphere(
        snow,
        snow,
        T_atmos,
        u_atmos,
        q_atmos,
        P_atmos,
        start_date,
        FT(10),
        toml_dict,
    )
    canopy_dry = Canopy.CanopyModel{FT}(
        domain,
        (; radiation, atmos = rain_off, ground),
        LAI,
        toml_dict;
        biomass,
        interception,
    )
    set_initial_cache_dry! = make_set_initial_cache(canopy_dry)
    exp_tendency_dry! = make_exp_tendency(canopy_dry)
    W0 = copy(parent(Y.canopy.interception.W))
    for step in 1:200
        set_initial_cache_dry!(p, Y, t0 + step * Δt)
        exp_tendency_dry!(dY, Y, p, t0 + step * Δt)
        @test all(parent(p.canopy.interception.intercepted_liq) .== 0)
        Y.canopy.interception.W .+= Δt .* dY.canopy.interception.W
        @test all(parent(Y.canopy.interception.W) .>= -sqrt(eps(FT)) * W_max)
    end
    @test all(parent(Y.canopy.interception.W) .< W0)
end
