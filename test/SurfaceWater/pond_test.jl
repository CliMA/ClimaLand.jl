using Test
using OrdinaryDiffEq: ODEProblem, solve, Euler
using ClimaCore

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Pond

FT = Float64
@testset "Pond integration tests" begin
    function precipitation(t::ft) where {ft}
        if t < ft(20)
            precip = -ft(1e-8)
        else
            precip = t < ft(100) ? -ft(5e-5) : ft(0.0)
        end
        return precip
    end

    pond_model = Pond.PondModel{FT}(;
        runoff = PrescribedRunoff{FT}(precipitation, (t) -> 0.0),
    )

    Y, p, coords = initialize(pond_model)
    Y.surface_water.η .= FT(0.0)

    pond_ode! = make_ode_function(pond_model)
    t0 = FT(0)
    tf = FT(200)
    dt = FT(1)

    prob = ODEProblem(pond_ode!, Y, (t0, tf), p)
    sol = solve(prob, Euler(), dt = dt)


    η = [sol.u[k].surface_water.η[1] for k in 1:1:length(sol.t)]
    t = sol.t[1:end]

    function expected_pond_height(t)
        if t < 20
            return 1e-8 * t
        elseif t < 100
            return (1e-8 * 20 + (t - 20.0) * (5e-5))
        else
            return (20 * (1e-8) + 80 * (5e-5))
        end
    end
    @test sum(abs.(expected_pond_height.(t) .- η)) < 1e-14
end
