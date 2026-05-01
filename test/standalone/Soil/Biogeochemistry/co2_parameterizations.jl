using Test
import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams as CP
using ClimaLand.Soil.Biogeochemistry

import ClimaLand
import ClimaLand.Parameters as LP

for FT in (Float32, Float64)
    @testset "Soil CO2 production and transport, FT = $FT" begin
        toml_dict = LP.create_toml_dict(FT)
        earth_param_set = LP.LandParameters(toml_dict)
        # Parameters should be supplied in m/kg/s (Pa... etc)
        ν = FT(0.556)
        # Prognostic variables
        P_sfc = FT(101e3)
        T_soil = FT(303)
        θ_l = FT(0.3)
        θ_i = FT(0.0)
        θ_w = θ_l + θ_i
        Csom = FT(5.0)
        T_ref = FT(LP.T_0(earth_param_set))
        R = FT(LP.gas_constant(earth_param_set))

        parameters = SoilCO2ModelParameters(toml_dict)

        # Test volumetric_air_content
        θ_a = volumetric_air_content(θ_w, ν)
        @test θ_a == ν - θ_w
        @test typeof(θ_a) == FT

        # Test co2_diffusivity
        θ_a100 = FT(0.3)
        b = FT(0.45)
        D = co2_diffusivity(T_soil, θ_w, P_sfc, θ_a100, b, ν, parameters)
        @test D ==
              (
                  parameters.D_ref *
                  (T_soil / T_ref)^FT(1.75) *
                  (FT(LP.P_ref(parameters.earth_param_set)) / P_sfc)
              ) *
              (FT(2)θ_a100^FT(3) + FT(0.04) * θ_a100) *
              (θ_a / θ_a100)^(FT(2) + FT(3) / b)
        @test typeof(D) == FT

        # Test o2_diffusivity: same Ryan/Moldrup form as CO₂, only the reference
        # diffusivity differs, so D_o2 / D_co2 == D_ref_o2 / D_ref.
        D_o2 = o2_diffusivity(T_soil, θ_w, P_sfc, θ_a100, b, ν, parameters)
        @test typeof(D_o2) == FT
        @test D_o2 ≈ D * parameters.D_ref_o2 / parameters.D_ref

        # Test O2 concentration conversions
        O2_f = FT(0.21)  # Typical atmospheric value
        ρ_O2_air = o2_concentration(O2_f, T_soil, P_sfc, parameters)
        @test typeof(ρ_O2_air) == FT
        @test ρ_O2_air > FT(0)

        # Test inverse: converting back from concentration to fraction
        O2_f_reconstructed =
            o2_fraction_from_concentration(ρ_O2_air, T_soil, P_sfc, parameters)
        @test O2_f_reconstructed ≈ O2_f
        @test typeof(O2_f_reconstructed) == FT

        # Test O2 availability (dimensionless, accounts for tortuosity)
        (; D_oa) = parameters
        O2_avail = o2_availability(O2_f, θ_a, D_oa)
        @test O2_avail == D_oa * O2_f * θ_a^(FT(4 / 3))
        @test typeof(O2_avail) == FT

        # Test microbe_source with O2_avail
        ms = microbe_source(T_soil, θ_l, Csom, O2_avail, parameters)
        # Check that the value is correct
        (; p_sx, D_liq, kM_o2, kM_sx, V_ref_sx, T_ref_sx, Ea_sx) = parameters
        Sx = p_sx * Csom * D_liq * max(θ_l, FT(0))^3
        MM_sx = Sx / (kM_sx + Sx)
        MM_o2 = O2_avail / (kM_o2 + O2_avail)
        @test ms ==
              V_ref_sx *
              exp(-Ea_sx / R * (1 / T_soil - 1 / T_ref_sx)) *
              MM_sx *
              MM_o2
        @test typeof(ms) == FT
        @test ms >= FT(0)  # Respiration should be non-negative

        # Test henry_constant: at the reference temperature, returns K_H_ref
        T_ref_henry = parameters.T_ref_henry
        K_H_co2_ref = henry_constant(
            parameters.K_H_co2_298,
            parameters.dln_K_H_co2_dT,
            T_ref_henry,
            T_ref_henry,
        )
        @test K_H_co2_ref ≈ parameters.K_H_co2_298
        @test typeof(K_H_co2_ref) == FT
        # K_H decreases with temperature (positive dln_K_H_dT means more soluble when colder)
        K_H_co2_warm = henry_constant(
            parameters.K_H_co2_298,
            parameters.dln_K_H_co2_dT,
            T_ref_henry + FT(15),
            T_ref_henry,
        )
        @test K_H_co2_warm < parameters.K_H_co2_298

        # Test beta_gas: β = K_H * R * T
        β_co2 = beta_gas(K_H_co2_ref, R, T_ref_henry)
        @test β_co2 ≈ parameters.K_H_co2_298 * R * T_ref_henry
        @test typeof(β_co2) == FT
        @test β_co2 > FT(0)

        # Test effective_porosity = θ_a + β·θ_l (with floor θ_eff_min = 1e-4)
        θ_eff = effective_porosity(θ_a, θ_l, β_co2)
        @test θ_eff ≈ θ_a + β_co2 * θ_l
        @test typeof(θ_eff) == FT
        # When both θ_a and θ_l vanish, the floor kicks in
        θ_eff_floor = effective_porosity(FT(0), FT(0), β_co2)
        @test θ_eff_floor ≈ FT(1e-4)
        # When only θ_a vanishes (saturated soil), dissolved storage keeps θ_eff > 0
        θ_eff_sat = effective_porosity(FT(0), FT(0.5), β_co2)
        @test θ_eff_sat ≈ β_co2 * FT(0.5)
        @test θ_eff_sat > FT(0)

        # Test gas_diffusivity_in_soil exponent uses params.T_exp_diffusivity
        @test parameters.T_exp_diffusivity == FT(1.75)
    end
end
