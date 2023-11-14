using ClimaLSM
using Test

FT = Float32
@testset "Base runoff functionality, FT = $FT" begin
    runoff = ClimaLSM.Soil.NoRunoff()
    precip = 5.0
    @test ClimaLSM.Soil.soil_surface_infiltration(runoff, precip) == precip
    @test ClimaLSM.Soil.subsurface_runoff_source(runoff) == nothing
    struct Foo{FT} <: ClimaLSM.Soil.AbstractSoilSource{FT} end
    srcs = (1, 2, 3)
    @test ClimaLSM.Soil.append_source(nothing, srcs) == srcs
    @test ClimaLSM.Soil.append_source(Foo{FT}(), srcs) == (srcs..., Foo{FT}())
end
