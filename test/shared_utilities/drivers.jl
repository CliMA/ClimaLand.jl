using ClimaCore
using Test
using StaticArrays
using ClimaLSM

FT = Float32
@testset "Default model, FT = $FT" begin
    pa = ClimaLSM.PrescribedAtmosphere(
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        FT(1);
    )
    pr = ClimaLSM.PrescribedRadiativeFluxes(FT, nothing, nothing, nothing)
    coords = (; surface = [1, 2, 3])
    @test ClimaLSM.initialize_drivers(nothing, nothing, coords) == (;)
    pa_keys = (:P_liq, :P_snow, :T, :P, :u, :q, :c_co2)
    pa_vals = ([zeros(FT, 3) for k in pa_keys]...,)
    @test ClimaLSM.initialize_drivers(pa, nothing, coords) ==
          NamedTuple{pa_keys}(pa_vals)
    papr_keys = (:P_liq, :P_snow, :T, :P, :u, :q, :c_co2, :SW_d, :LW_d, :θs)
    papr_vals = ([zeros(FT, 3) for k in papr_keys]...,)
    @test ClimaLSM.initialize_drivers(pa, pr, coords) ==
          NamedTuple{papr_keys}(papr_vals)
end
