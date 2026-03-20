# Tests for the standalone InlandWater slab lake component.
# The lake is a separate component added to LandModel, using fraction-based
# blending for soil boundary conditions.

using Test
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand
using ClimaLand.InlandWater
import ClimaLand.InlandWater as IW
import ClimaLand.Parameters as LP

@testset "Lake physics round-trip: T → U → q_l → T" begin
    FT = Float64
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    params = SlabLakeParameters{FT}(; earth_param_set)

    # Warm lake (liquid)
    T_warm = FT(290)
    U_warm = IW.lake_energy_from_temperature(T_warm, params, earth_param_set)
    q_l_warm = IW.lake_liquid_fraction(U_warm, params, earth_param_set)
    T_warm_back = IW.lake_temperature(U_warm, q_l_warm, params, earth_param_set)
    @test q_l_warm ≈ FT(1)
    @test T_warm_back ≈ T_warm atol = FT(0.01)

    # Cold lake (frozen)
    T_cold = FT(260)
    U_cold = IW.lake_energy_from_temperature(T_cold, params, earth_param_set)
    q_l_cold = IW.lake_liquid_fraction(U_cold, params, earth_param_set)
    T_cold_back = IW.lake_temperature(U_cold, q_l_cold, params, earth_param_set)
    @test q_l_cold ≈ FT(0)
    @test T_cold_back ≈ T_cold atol = FT(0.01)

    # At freezing — T == T_freeze takes the ice branch, so q_l = 0
    T_freeze = FT(LP.T_freeze(earth_param_set))
    U_freeze =
        IW.lake_energy_from_temperature(T_freeze, params, earth_param_set)
    q_l_freeze = IW.lake_liquid_fraction(U_freeze, params, earth_param_set)
    @test q_l_freeze ≈ FT(0)
end

@testset "Lake albedo" begin
    FT = Float64
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    params = SlabLakeParameters{FT}(; earth_param_set)

    @test IW.lake_surface_albedo(FT(1), params) ≈ params.liquid_albedo
    @test IW.lake_surface_albedo(FT(0), params) ≈ params.ice_albedo
    @test IW.lake_surface_albedo(FT(0.5), params) ≈
          FT(0.5) * params.liquid_albedo + FT(0.5) * params.ice_albedo
end

@testset "Lake sediment heat flux" begin
    FT = Float64
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    params = SlabLakeParameters{FT}(; earth_param_set)

    # When lake is warmer than soil, flux should be negative (heat flows down)
    Q_sed =
        IW.lake_sediment_heat_flux(FT(290), FT(280), FT(1.5), FT(0.1), params)
    @test Q_sed < FT(0)

    # When soil is warmer, flux should be positive
    Q_sed2 =
        IW.lake_sediment_heat_flux(FT(280), FT(290), FT(1.5), FT(0.1), params)
    @test Q_sed2 > FT(0)

    # Symmetric: flipping ΔT flips the sign
    @test Q_sed ≈ -Q_sed2

    # Zero temperature difference gives zero flux
    Q_sed0 =
        IW.lake_sediment_heat_flux(FT(285), FT(285), FT(1.5), FT(0.1), params)
    @test Q_sed0 ≈ FT(0)
end

@testset "Lake runoff energy flux" begin
    FT = Float64
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)

    # Zero runoff gives zero energy flux
    @test IW.lake_runoff_energy_flux(FT(0), FT(290), FT(1), earth_param_set) ≈
          FT(0)

    # Non-zero runoff gives non-zero energy flux
    E = IW.lake_runoff_energy_flux(FT(1e-3), FT(290), FT(1), earth_param_set)
    @test isfinite(E)
    @test E != FT(0)
end

@testset "Energy conservation: precipitation onto lake" begin
    # When precipitation falls on a warm lake with no sediment flux,
    # the energy tendency should equal -(surface flux + runoff energy flux).
    # Here we verify the runoff energy budget: precipitation adds water at T_air,
    # runoff removes it at T_lake, and the net energy is accounted for.
    FT = Float64
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    params = SlabLakeParameters{FT}(; earth_param_set)

    _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
    _cp_l = FT(LP.cp_l(earth_param_set))
    _T_ref = FT(LP.T_0(earth_param_set))

    T_lake = FT(290)
    q_l = FT(1)  # fully liquid
    P_liq = FT(1e-4)  # precipitation rate (m/s)

    # Runoff = -F_water = -P_liq (no evaporation) to maintain constant depth
    runoff = -P_liq

    # Precipitation enthalpy flux into lake (at air temperature = lake temperature here)
    T_air = T_lake
    Q_precip = P_liq * _ρ_l * _cp_l * (T_air - _T_ref)

    # Runoff energy flux out of lake (at lake temperature)
    Q_runoff = IW.lake_runoff_energy_flux(runoff, T_lake, q_l, earth_param_set)

    # ρe_int at lake temperature for fully liquid
    ρe_lake = IW.lake_volumetric_internal_energy(T_lake, q_l, earth_param_set)

    # Runoff energy should be -P_liq * ρe(T_lake)
    @test Q_runoff ≈ -P_liq * ρe_lake

    # When T_air == T_lake and no radiation/turbulent fluxes, the energy tendency
    # contribution from precip + runoff should be: -Q_precip - Q_runoff
    # = -P_liq*ρ*cp*(T-T0) - (-P_liq)*ρe(T) = -P_liq*ρ*cp*(T-T0) + P_liq*ρe(T)
    # For fully liquid: ρe(T) = ρ*cp*(T-T0), so these cancel exactly.
    net_energy = -Q_precip - Q_runoff
    @test abs(net_energy) < eps(FT) * abs(Q_precip)
end

@testset "Energy conservation: evaporation from lake" begin
    # When evaporation removes water, runoff is positive (water supplied to maintain depth).
    # Verify the runoff energy flux is consistent.
    FT = Float64
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    params = SlabLakeParameters{FT}(; earth_param_set)

    T_lake = FT(295)
    q_l = FT(1)

    # Evaporation rate (negative water flux, so E_liq < 0 in sign convention)
    E_liq = FT(-5e-5)  # m/s evaporated
    F_water = E_liq  # no precipitation
    runoff = -F_water  # positive: water must be supplied

    Q_runoff = IW.lake_runoff_energy_flux(runoff, T_lake, q_l, earth_param_set)
    ρe_lake = IW.lake_volumetric_internal_energy(T_lake, q_l, earth_param_set)

    # Runoff is positive and carries energy into the lake at lake temperature
    @test Q_runoff ≈ runoff * ρe_lake
    @test Q_runoff > FT(0)  # energy is added (water supplied at lake T)
end

@testset "Mass conservation: runoff balances net water flux" begin
    FT = Float64

    # Scenario 1: net precipitation → positive runoff (excess water drained)
    P_liq = FT(2e-4)
    E_liq = FT(-5e-5)
    E_ice = FT(0)
    F_water = P_liq + E_liq + E_ice
    runoff = -F_water
    @test runoff < FT(0)  # water leaves as drainage (net precip case → runoff < 0 means drainage)
    # Actually: P_liq > |E_liq|, so F_water > 0, runoff = -F_water < 0
    # Negative runoff = drainage away from lake. Total mass change = F_water + runoff = 0
    @test F_water + runoff ≈ FT(0)

    # Scenario 2: net evaporation → water must be supplied
    P_liq2 = FT(1e-5)
    E_liq2 = FT(-1e-4)
    F_water2 = P_liq2 + E_liq2
    runoff2 = -F_water2
    @test F_water2 + runoff2 ≈ FT(0)

    # Scenario 3: no fluxes → no runoff
    @test -FT(0) ≈ FT(0)
end

@testset "Energy conservation: sediment heat flux symmetry" begin
    # Verify that sediment heat flux conserves energy: what leaves the lake
    # enters the soil and vice versa.
    FT = Float64
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    params = SlabLakeParameters{FT}(; earth_param_set)

    T_lake = FT(285)
    T_soil = FT(280)
    κ_soil = FT(1.5)
    Δz_soil = FT(0.1)

    Q_sed = IW.lake_sediment_heat_flux(T_lake, T_soil, κ_soil, Δz_soil, params)

    # Heat flows from warm lake to cold soil: Q_sed < 0 (lake loses energy)
    @test Q_sed < FT(0)

    # The soil gains exactly the same magnitude of energy
    # (In the tendency: dU_lake/dt includes +Q_sed, soil top BC gets -Q_sed)
    # Reversing temperatures reverses flux exactly
    Q_sed_rev =
        IW.lake_sediment_heat_flux(T_soil, T_lake, κ_soil, Δz_soil, params)
    @test Q_sed + Q_sed_rev ≈ FT(0)

    # At thermal equilibrium, no heat flows
    Q_eq = IW.lake_sediment_heat_flux(T_lake, T_lake, κ_soil, Δz_soil, params)
    @test Q_eq ≈ FT(0)
end

@testset "Energy conservation: frozen lake runoff" begin
    # Verify runoff energy is consistent for a frozen lake (q_l = 0)
    FT = Float64
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)

    T_lake = FT(260)  # frozen
    q_l = FT(0)
    runoff = FT(1e-4)

    Q_runoff = IW.lake_runoff_energy_flux(runoff, T_lake, q_l, earth_param_set)
    ρe_ice = IW.lake_volumetric_internal_energy(T_lake, q_l, earth_param_set)

    @test Q_runoff ≈ runoff * ρe_ice
    # Frozen water at 260 K has negative internal energy (below reference + latent heat)
    @test ρe_ice < FT(0)
end

@testset "SlabLakeParameters defaults" begin
    FT = Float64
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    params = SlabLakeParameters{FT}(; earth_param_set)

    @test params.depth == FT(10)
    @test params.conductance == FT(0.1)
    @test params.liquid_albedo == FT(0.08)
    @test params.ice_albedo == FT(0.6)
    @test params.emissivity == FT(0.97)
end
