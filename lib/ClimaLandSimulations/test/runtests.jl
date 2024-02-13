using Test
using ClimaLandSimulations
using ClimaLandSimulations.Fluxnet

@testset "Fluxnet single site" begin
    sv, sol, Y, p = run_fluxnet("US-MOz")
    @test typeof(sv) ==
          @NamedTuple{t::Vector{Float64}, saveval::Vector{NamedTuple}}
end
