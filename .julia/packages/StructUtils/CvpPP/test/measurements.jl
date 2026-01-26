using StructUtils, Measurements, Test

@testset "measurements.jl" begin
    m = 1.1 ± 3.3
    str = StructUtils.lower(m)
    @test str == "1.1 ± 3.3"
    m2 = StructUtils.lift(Measurements.Measurement{Float64}, str)
    @test m2 == m
end