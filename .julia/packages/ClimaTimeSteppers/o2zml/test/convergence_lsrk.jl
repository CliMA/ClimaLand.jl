#=
julia --project=test
using Revise; include("test/convergence_lsrk.jl")
=#
using ClimaTimeSteppers, LinearAlgebra, Test
import PrettyTables
import SciMLBase

if !@isdefined(IntegratorTestCase)
    include(joinpath(@__DIR__, "convergence_orders.jl"))
    include(joinpath(@__DIR__, "convergence_utils.jl"))
    include(joinpath(@__DIR__, "utils.jl"))
    include(joinpath(@__DIR__, "problems.jl"))
end

@testset "LSRK convergence" begin
    dts = 0.5 .^ (4:7)

    for (prob, sol, tscale) in [
        (linear_prob(), linear_sol, 1)
        (sincos_prob(), sincos_sol, 1)
    ]

        @test convergence_order(prob, sol, LSRKEulerMethod(), dts .* tscale) ≈ 1 atol = 0.1
        @test convergence_order(prob, sol, LSRK54CarpenterKennedy(), dts .* tscale) ≈ 4 atol = 0.05
        @test convergence_order(prob, sol, LSRK144NiegemannDiehlBusch(), dts .* tscale) ≈ 4 atol = 0.05
    end
end
