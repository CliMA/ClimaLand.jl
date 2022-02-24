using OrdinaryDiffEq: ODEProblem, solve, RK4
using Test
using DifferentialEquations
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Carbon

@testset "Carbon equlibrium" begin
    NPP_model =
        PrescribedNPPModel{Float64}((t) -> 60 + 5 * sin(2 * pi * t / 10)) # PgC/year
    carbon = Clima2PoolDefault(Float64, NPP_model)
    Y, p, coords = initialize(carbon)

    Y.carbon.pool[1] = 600.0
    Y.carbon.pool[2] = 1500.0
    ode_function! = make_ode_function(carbon)
    dY = similar(Y)
    ode_function!(dY, Y, p, 0.0)
    @test sum(abs.(dY.carbon.pool .- [0.0, 0.0])) < 1e-14
end
