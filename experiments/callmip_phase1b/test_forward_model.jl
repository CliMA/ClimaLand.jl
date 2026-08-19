# Phase 4 smoke test: one G(θ) evaluation (default params, years = [2010]).
# Checks length/layout, finiteness, and physical plausibility per flux block.
import ClimaComms
ClimaComms.@import_required_backends
using Test, Dates, Statistics

include(joinpath(@__DIR__, "forward_model_phase1b.jl"))

years = [2010]
t = @elapsed g = run_member(Float64[], years)
@info "run_member done" wall_min = round(t / 60; digits = 1) len = length(g)

@testset "G smoke (default params, 2010)" begin
    @test length(g) == 756
    @test count(isnan, g) == 0
    G = reshape(g, 36, N_SITES)      # per site: [nee12; lhf12; shf12]
    nee, lhf, shf = G[1:12, :], G[13:24, :], G[25:36, :]
    @info "block ranges" nee = extrema(nee) lhf = extrema(lhf) shf =
        extrema(shf)
    @test all(-30 .< nee .< 30)      # gC/m2/d, generous bounds
    @test all(-50 .< lhf .< 400)     # W/m2
    @test all(-150 .< shf .< 400)
    # coarse physics: growing-season NEE (uptake, negative) at DK-Sor (col 6)
    @info "DK-Sor monthly NEE (gC/m2/d)" round.(nee[:, 6]; digits = 2)
end
@info "test_forward_model.jl finished"
