using Test
using ClimaLandSimulations
using ClimaLandSimulations.Fluxnet

@testset "Fluxnet single site" begin
    sv, sol, Y, p = run_fluxnet("US-MOz")
    @test typeof(sv) ==
          @NamedTuple{t::Vector{Float64}, saveval::Vector{NamedTuple}}
end

@testset "Fluxnet single site, custom param" begin
    sv, sol, Y, p = run_fluxnet(
        "US-MOz";
        params = ozark_default_params(;
            hetero_resp = hetero_resp_ozark(; b = 2),
        ),
    )
    @test typeof(sv) ==
          @NamedTuple{t::Vector{Float64}, saveval::Vector{NamedTuple}}
end

@testset "Fluxnet single site, custom start time" begin
    sv, sol, Y, p = run_fluxnet(
        "US-MOz";
        setup = make_setup(;
            ozark_default_args(; t0 = Float64(125 * 3600 * 24))...,
        ),
    )
    @test typeof(sv) ==
          @NamedTuple{t::Vector{Float64}, saveval::Vector{NamedTuple}}
end
