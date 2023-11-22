using Test
using ClimaCore

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Pond

for FT in (Float32, Float64)
    @testset "Pond integration tests, FT = $FT" begin
        function precipitation(t)
            if t < 20
                precip = -1e-8
            else
                precip = t < 100 ? -5e-5 : 0.0
            end
            return precip
        end

        pond_model = Pond.PondModel{FT}(;
            runoff = PrescribedRunoff{FT}(precipitation, (t) -> 0.0),
        )

        Y, p, coords = initialize(pond_model)
        Y.surface_water.η .= FT(0.0)

        exp_tendency! = make_exp_tendency(pond_model)
        t0 = Float64(0)
        tf = Float64(200)
        dt = Float64(1)
        set_initial_aux_state! = make_set_initial_aux_state(pond_model)
        set_initial_aux_state!(p, Y, t0)


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
            exp_tendency!(dY, Y, p, t)
            @. Y.surface_water.η = dY.surface_water.η * dt + Y.surface_water.η
            t = t + dt
            if t % 40 == 0
                @test abs.(
                    FT(expected_pond_height(t)) .-
                    Array(parent(Y.surface_water.η))[1]
                ) < eps(FT)
            end
        end
    end
end
