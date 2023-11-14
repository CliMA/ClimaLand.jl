using Test
import CLIMAParameters as CP
using ClimaLSM.Soil.Biogeochemistry

import ClimaLSM
import ClimaLSM.Parameters as LSMP

include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

for FT in (Float32, Float64)
    @testset "Soil CO2 production and transport, FT = $FT" begin
        earth_param_set = create_lsm_parameters(FT)
        # Parameters should be supplied in m/kg/s (Pa... etc)
        D_liq = FT(3.17)
        ν = FT(0.556)
        θ_a100 = FT(0.1846)
        D_ref = FT(1.39e-5)
        b = FT(4.547)
        α_sx = FT(194e3)
        Ea_sx = FT(61e3)
        kM_sx = FT(5e-3)
        kM_o2 = FT(0.004)
        O2_a = FT(0.209)
        D_oa = FT(1.67)
        p_sx = FT(0.024)
        # Prognostic variables
        P_sfc = FT(101e3)
        T_soil = FT(303)
        θ_l = FT(0.3)
        θ_i = FT(0.0)
        θ_w = θ_l + θ_i
        Csom = FT(5.0)
        T_ref = FT(LSMP.T_0(earth_param_set))
        R = FT(LSMP.gas_constant(earth_param_set))

        parameters = SoilCO2ModelParameters{FT}(;
            D_liq = D_liq,
            ν = ν,
            θ_a100 = θ_a100,
            D_ref = D_ref,
            b = b,
            earth_param_set = earth_param_set,
        )

        # Test that parameterizations functions are working properly
        θ_a = volumetric_air_content(θ_w, parameters)
        @test θ_a == parameters.ν - θ_w
        @test typeof(θ_a) == FT

        D = co2_diffusivity(T_soil, θ_w, P_sfc, parameters)
        @test D ==
              (
                  parameters.D_ref *
                  (T_soil / T_ref)^FT(1.75) *
                  (FT(LSMP.P_ref(parameters.earth_param_set)) / P_sfc)
              ) *
              (FT(2)parameters.θ_a100^FT(3) + FT(0.04)parameters.θ_a100) *
              (θ_a / parameters.θ_a100)^(FT(2) + FT(3) / parameters.b)
        @test typeof(D) == FT

        ms = microbe_source(T_soil, θ_l, Csom, parameters)
        # check that the value is correct
        MM_sx =
            p_sx * Csom * D_liq * θ_l^3 / (kM_sx + p_sx * Csom * D_liq * θ_l^3)
        MM_o2 =
            D_oa * O2_a * ((ν - θ_l)^(FT(4 / 3))) /
            (kM_o2 + D_oa * O2_a * ((ν - θ_l)^(FT(4 / 3))))
        @test ms == α_sx * exp(-Ea_sx / (R * T_soil)) * MM_sx * MM_o2
        @test typeof(ms) == FT
    end
end
