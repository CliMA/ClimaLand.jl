using ClimaLSM
using Test

@testset "Base runoff functionality" begin
    runoff = ClimaLSM.Soil.NoRunoff()
    precip = 5.0
    @test ClimaLSM.Soil.soil_surface_infiltration(runoff, precip) == precip
    @test ClimaLSM.Soil.subsurface_runoff_source(runoff) == nothing
    struct Foo{Float64} <: ClimaLSM.Soil.AbstractSoilSource{Float64} end
    srcs = (1, 2, 3)
    @test ClimaLSM.Soil.append_source(nothing, srcs) == srcs
    @test ClimaLSM.Soil.append_source(Foo{Float64}(), srcs) ==
          (srcs..., Foo{Float64}())
end
