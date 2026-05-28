"""
Test that the SoilCO2 model remains stable under saturated/near-saturated conditions.

The Henry's law air-water partitioning scheme should prevent blow-up when:
- θ_l → ν (fully saturated soil, no air)
- θ_l + θ_i → ν (saturated with ice)

Key invariants:
- CO2, O2, and their tendencies should never be NaN or Inf
- θ_eff and θ_eff_o2 should remain positive even when θ_a = 0
"""

using Test
import ClimaComms
ClimaComms.@import_required_backends
using Dates
using ClimaCore
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil.Biogeochemistry
import ClimaLand.Parameters as LP

FT = Float64

@testset "Henry's law parameterizations" begin
    # Test henry_constant
    K_H_298 = FT(3.4e-4)  # CO2 at 298K
    dln_K_H_dT = FT(2400)

    # At reference temperature, should return K_H_298
    T_ref = FT(298.15)
    @test henry_constant(K_H_298, dln_K_H_dT, T_ref, T_ref) ≈ K_H_298

    # At lower temperature, K_H should be higher (more soluble)
    T_cold = FT(273.15)
    K_H_cold = henry_constant(K_H_298, dln_K_H_dT, T_cold, T_ref)
    @test K_H_cold > K_H_298

    # At higher temperature, K_H should be lower
    T_warm = FT(313.15)
    K_H_warm = henry_constant(K_H_298, dln_K_H_dT, T_warm, T_ref)
    @test K_H_warm < K_H_298

    # Test beta_gas
    R = FT(8.314)
    T = FT(298.15)
    K_H = FT(3.4e-4)
    β = beta_gas(K_H, R, T)
    @test β ≈ K_H * R * T
    @test β > FT(0)

    # Test effective_porosity
    θ_a = FT(0.2)
    θ_l = FT(0.3)
    β_co2 = FT(0.8)  # typical value for CO2
    θ_eff = effective_porosity(θ_a, θ_l, β_co2)
    @test θ_eff ≈ θ_a + β_co2 * θ_l

    # Critical test: when θ_a = 0, θ_eff should still be positive
    θ_a_zero = FT(0)
    θ_l_sat = FT(0.5)
    θ_eff_sat = effective_porosity(θ_a_zero, θ_l_sat, β_co2)
    @test θ_eff_sat > FT(0)
    @test θ_eff_sat ≈ β_co2 * θ_l_sat
end

@testset "volumetric_air_content physical floor" begin
    ν = FT(0.5)

    # Normal case
    θ_w_normal = FT(0.3)
    θ_a = volumetric_air_content(θ_w_normal, ν)
    @test θ_a ≈ ν - θ_w_normal

    # Saturated case - should return 0, not 0.001
    θ_w_sat = ν
    θ_a_sat = volumetric_air_content(θ_w_sat, ν)
    @test θ_a_sat == FT(0)

    # Over-saturated case (shouldn't happen physically, but should be robust)
    θ_w_over = ν + FT(0.1)
    θ_a_over = volumetric_air_content(θ_w_over, ν)
    @test θ_a_over == FT(0)
end

@testset "SoilCO2 model stability under saturation" begin
    toml_dict = LP.create_toml_dict(FT)
    parameters = SoilCO2ModelParameters(toml_dict)

    # Verify Henry's law parameters are set
    @test parameters.K_H_co2_298 > FT(0)
    @test parameters.K_H_o2_298 > FT(0)
    @test parameters.dln_K_H_co2_dT > FT(0)
    @test parameters.dln_K_H_o2_dT > FT(0)

    # Set up a single column domain
    nelems = 10
    zmin = FT(-1)
    zmax = FT(0.0)
    domain = Column(; zlim = (zmin, zmax), nelements = nelems)

    sources = (MicrobeProduction{FT}(),)

    # Set up atmosphere
    precipitation_function = (t) -> 0.0
    snow_precip = (t) -> 0.0
    atmos_T = (t) -> 288.0
    atmos_u = (t) -> 2.0
    atmos_q = (t) -> 0.01
    atmos_p = (t) -> 101325.0
    UTC_DATETIME = Dates.now()
    atmos_h = FT(2)
    atmos_co2 = (t) -> 0.0004  # ~400 ppm as mass concentration

    atmos = ClimaLand.PrescribedAtmosphere(
        TimeVaryingInput(precipitation_function),
        TimeVaryingInput(snow_precip),
        TimeVaryingInput(atmos_T),
        TimeVaryingInput(atmos_u),
        TimeVaryingInput(atmos_q),
        TimeVaryingInput(atmos_p),
        UTC_DATETIME,
        atmos_h,
        toml_dict;
        c_co2 = TimeVaryingInput(atmos_co2),
    )

    ν = FT(0.5)
    θ_r = FT(0.05)
    α = FT(0.1)
    n = FT(2)
    hcm = ClimaLand.Soil.vanGenuchten{FT}(; α = α, n = n)

    # Test different saturation scenarios
    test_cases = [
        ("Normal conditions", (z, t) -> eltype(z)(0.25)),
        ("Near-saturated", (z, t) -> eltype(z)(0.49)),
        ("Fully saturated", (z, t) -> eltype(z)(0.5)),  # = ν
        ("Over-saturated (edge case)", (z, t) -> eltype(z)(0.51)),  # > ν
    ]

    for (name, θ_l_func) in test_cases
        @testset "$name" begin
            T_soil = (z, t) -> eltype(z)(288)
            prescribed_met = PrescribedMet{FT}(T_soil, θ_l_func, ν, θ_r, hcm)
            drivers = SoilDrivers(prescribed_met, atmos)

            model = SoilCO2Model{FT}(
                domain,
                drivers,
                toml_dict,
                FT(0);
                parameters,
                sources,
            )

            Y, p, coords = initialize(model)
            dY = similar(Y)
            exp_tendency! = make_exp_tendency(model)
            t = Float64(0)

            # Set initial conditions
            Y.soilco2.CO2 .= FT(0.001)  # Small positive CO2
            Y.soilco2.O2 .= FT(0.21)
            Y.soilco2.SOC .= FT(5.0)    # Soil organic carbon

            # Update cache
            set_initial_cache! = make_set_initial_cache(model)
            set_initial_cache!(p, Y, t)

            # Check auxiliary variables are finite
            @test all(isfinite.(parent(p.soilco2.θ_eff)))
            @test all(isfinite.(parent(p.soilco2.θ_eff_o2)))
            @test all(isfinite.(parent(p.soilco2.D)))
            @test all(isfinite.(parent(p.soilco2.Sm)))

            # θ_eff should be positive even when θ_a = 0
            @test all(parent(p.soilco2.θ_eff) .> FT(0))
            @test all(parent(p.soilco2.θ_eff_o2) .> FT(0))

            # Compute tendency
            exp_tendency!(dY, Y, p, t)

            # Check tendencies are finite (no blow-up)
            @test all(isfinite.(parent(dY.soilco2.CO2)))
            @test all(isfinite.(parent(dY.soilco2.O2)))
            @test all(isfinite.(parent(dY.soilco2.SOC)))

            # Check prognostic variables remain finite
            @test all(isfinite.(parent(Y.soilco2.CO2)))
            @test all(isfinite.(parent(Y.soilco2.O2)))
        end
    end
end

@testset "SoilCO2 stability with temperature variation" begin
    toml_dict = LP.create_toml_dict(FT)
    parameters = SoilCO2ModelParameters(toml_dict)

    nelems = 10
    zmin = FT(-1)
    zmax = FT(0.0)
    domain = Column(; zlim = (zmin, zmax), nelements = nelems)

    sources = (MicrobeProduction{FT}(),)

    # Set up atmosphere
    atmos = ClimaLand.PrescribedAtmosphere(
        TimeVaryingInput((t) -> 0.0),
        TimeVaryingInput((t) -> 0.0),
        TimeVaryingInput((t) -> 288.0),
        TimeVaryingInput((t) -> 2.0),
        TimeVaryingInput((t) -> 0.01),
        TimeVaryingInput((t) -> 101325.0),
        Dates.now(),
        FT(2),
        toml_dict;
        c_co2 = TimeVaryingInput((t) -> 0.0004),
    )

    ν = FT(0.5)
    θ_r = FT(0.05)
    hcm = ClimaLand.Soil.vanGenuchten{FT}(; α = FT(0.1), n = FT(2))

    # Test different temperatures
    test_temps = [FT(263), FT(273), FT(288), FT(303), FT(313)]

    for T_val in test_temps
        @testset "Temperature = $T_val K" begin
            T_soil = (z, t) -> eltype(z)(T_val)
            θ_l = (z, t) -> eltype(z)(0.45)  # Near-saturated
            prescribed_met = PrescribedMet{FT}(T_soil, θ_l, ν, θ_r, hcm)
            drivers = SoilDrivers(prescribed_met, atmos)

            model = SoilCO2Model{FT}(
                domain,
                drivers,
                toml_dict,
                FT(0);
                parameters,
                sources,
            )

            Y, p, coords = initialize(model)
            dY = similar(Y)
            exp_tendency! = make_exp_tendency(model)
            t = Float64(0)

            Y.soilco2.CO2 .= FT(0.001)
            Y.soilco2.O2 .= FT(0.21)
            Y.soilco2.SOC .= FT(5.0)

            set_initial_cache! = make_set_initial_cache(model)
            set_initial_cache!(p, Y, t)
            exp_tendency!(dY, Y, p, t)

            # All values should be finite
            @test all(isfinite.(parent(p.soilco2.θ_eff)))
            @test all(isfinite.(parent(dY.soilco2.CO2)))
            @test all(isfinite.(parent(dY.soilco2.O2)))
        end
    end
end

println("All saturation stability tests passed!")
