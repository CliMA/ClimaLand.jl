using Test
using ClimaCore

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Pond

FT = Float64
@testset "Pond integration tests" begin
    function precipitation(t::T) where {T}
        if t < T(20)
            precip = -T(1e-8)
        else
            precip = t < T(100) ? T(-5e-5) : T(0.0)
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

    function expected_pond_height(t)
        if t < 20
            return 1e-8 * t
        elseif t < 100
            return (1e-8 * 20 + (t - 20.0) * (5e-5))
        else
            return (20 * (1e-8) + 80 * (5e-5))
        end
    end

    nsteps = Int((tf - t0) / dt)
    dY = similar(Y)
    t = t0
    for step in 1:nsteps
        pond_ode!(dY, Y, p, t)
        @. Y.surface_water.η = dY.surface_water.η * dt + Y.surface_water.η
        t = t + dt
        if t % 40 == 0
            @test abs.(
                expected_pond_height(t) .- parent(Y.surface_water.η)[1]
            ) < 1e-14
        end

    end

end
