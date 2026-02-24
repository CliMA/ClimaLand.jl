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
