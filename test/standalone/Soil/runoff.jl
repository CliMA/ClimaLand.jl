using ClimaLand
using Test

FT = Float32
@testset "Base runoff functionality, FT = $FT" begin
    runoff = ClimaLand.Soil.NoRunoff()
    precip = 5.0
    @test ClimaLand.Soil.soil_surface_infiltration(runoff, precip) == precip
    @test ClimaLand.Soil.subsurface_runoff_source(runoff) == nothing
    struct Foo{FT} <: ClimaLand.Soil.AbstractSoilSource{FT} end
    srcs = (1, 2, 3)
    @test ClimaLand.Soil.append_source(nothing, srcs) == srcs
    @test ClimaLand.Soil.append_source(Foo{FT}(), srcs) == (srcs..., Foo{FT}())
end
