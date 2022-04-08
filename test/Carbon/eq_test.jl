using OrdinaryDiffEq: ODEProblem, solve, RK4
using Test
using DifferentialEquations
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Carbon
using ClimaLSM.Domains

@testset "Carbon equlibrium" begin
    test_domain = Column(; zlim = (-1.0, 0.0), nelements = 3)
    NPP_model =
        PrescribedNPPModel{Float64}((t) -> 60 + 5 * sin(2 * pi * t / 10)) # PgC/year
    carbon = Clima2PoolDefault(Float64, test_domain, NPP_model)
    Y, p, coords = initialize(carbon)

    #= how do you intialize? why is this so hard? 
    function init_carbon(z::FT)::SVector{2, FT} where {FT<:AbstractFloat}
        return SVector{2, FT}([600.0, 1500.0])
    end

    Y.carbon.pool = init_carbon.(coords.z)
    =#
    ode_function! = make_ode_function(carbon)
    dY = similar(Y)
    ode_function!(dY, Y, p, 0.0)
    @test parent(dY.carbon.pool) == [60.0 0.0; 60.0 0.0; 60.0 0.0]
end
