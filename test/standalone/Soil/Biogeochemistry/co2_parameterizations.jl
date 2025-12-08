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
        (; p_sx, D_liq, kM_o2, kM_sx, α_sx, Ea_sx) = parameters
        Sx = p_sx * Csom * D_liq * max(θ_l, FT(0))^3
        MM_sx = Sx / (kM_sx + Sx)
        MM_o2 = O2_avail / (kM_o2 + O2_avail)
        @test ms == α_sx * exp(-Ea_sx / (R * T_soil)) * MM_sx * MM_o2
        @test typeof(ms) == FT
        @test ms >= FT(0)  # Respiration should be non-negative
    end
end
