using Test
using LinearAlgebra
using ClimaCore
import CLIMAParameters as CP
using ClimaLSM
using ClimaLSM.Domains: Column, HybridBox
using ClimaLSM.Soil

import ClimaLSM
import ClimaLSM.Parameters as LSMP

include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

@testset "Richards Jacobian entries, Moisture BC" begin
    FT = Float64
    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(1.43)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0.124)
    zmax = FT(0)
    zmin = FT(-1.5)
    nelems = 150
    soil_domains = [
        Column(; zlim = (zmin, zmax), nelements = nelems),
        HybridBox(;
            xlim = (0.0, 1.0),
            ylim = (0.0, 1.0),
            zlim = (zmin, zmax),
            nelements = (1, 1, nelems),
            npolynomial = 3,
        ),
    ]
    top_state_bc = MoistureStateBC((p, t) -> eltype(t)(ν - 1e-3))
    bot_flux_bc = FreeDrainage()
    sources = ()
    boundary_states =
        (; top = (water = top_state_bc,), bottom = (water = bot_flux_bc,))
    params = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r)

    for domain in soil_domains
        soil = Soil.RichardsModel{FT}(;
            parameters = params,
            domain = domain,
            boundary_conditions = boundary_states,
            sources = sources,
        )

        Y, p, coords = initialize(soil)
        Y.soil.ϑ_l .= FT(0.24)
        W = RichardsTridiagonalW(Y)
        Wfact! = make_tendency_jacobian(soil)
        dtγ = FT(1.0)
        Wfact!(W, Y, p, dtγ, FT(0.0))

        K_ic = hydraulic_conductivity(
            K_sat,
            vg_m,
            effective_saturation(ν, FT(0.24), θ_r),
        )
        dz = FT(0.01)
        dψdϑ_ic = dψdϑ(FT(0.24), ν, θ_r, vg_α, vg_n, vg_m, S_s)

        @test all(parent(W.temp2) .== FT(0.0))
        @test all(parent(W.temp2) .== FT(0.0))
        @test W.transform == false
        @test typeof(W.W_column_arrays) <:
              Vector{LinearAlgebra.Tridiagonal{Float64, Vector{Float64}}}
        @test length(W.W_column_arrays) == 1
        if typeof(domain) <: Column
            @test all(
                parent(W.∂ϑₜ∂ϑ.coefs.:1)[2:end] .≈
                parent(W.∂ϑₜ∂ϑ.coefs.:3)[1:(end - 1)],
            )
            @test parent(W.∂ϑₜ∂ϑ.coefs.:1)[1] == FT(0)
            @test parent(W.∂ϑₜ∂ϑ.coefs.:3)[end] == FT(0)
            @test all(parent(W.∂ϑₜ∂ϑ.coefs.:1)[2:end] .≈ K_ic / dz^2 * dψdϑ_ic)
            @test parent(W.∂ϑₜ∂ϑ.coefs.:2)[1] .≈ -K_ic / dz^2 * dψdϑ_ic
            @test all(
                parent(W.∂ϑₜ∂ϑ.coefs.:2)[2:(end - 1)] .≈
                -2 * K_ic / dz^2 * dψdϑ_ic,
            )
            @test parent(W.∂ϑₜ∂ϑ.coefs.:2)[end] .≈
                  -K_ic / dz^2 * dψdϑ_ic - K_ic / (dz * dz / 2) * dψdϑ_ic
        elseif typeof(domain) <: HybridBox
            @test all(
                parent(W.∂ϑₜ∂ϑ.coefs.:1)[2:end, :, 1, 1, 1] .≈
                parent(W.∂ϑₜ∂ϑ.coefs.:3)[1:(end - 1), :, 1, 1, 1],
            )
            @test all(parent(W.∂ϑₜ∂ϑ.coefs.:1)[1, :, 1, 1, 1] .== FT(0))
            @test all(parent(W.∂ϑₜ∂ϑ.coefs.:3)[end, :, 1, 1, 1] .== FT(0))
            @test all(
                parent(W.∂ϑₜ∂ϑ.coefs.:1)[2:end, :, 1, 1, 1] .≈
                K_ic / dz^2 * dψdϑ_ic,
            )
            @test all(
                parent(W.∂ϑₜ∂ϑ.coefs.:2)[1, :, 1, 1, 1] .≈
                -K_ic / dz^2 * dψdϑ_ic,
            )
            @test all(
                parent(W.∂ϑₜ∂ϑ.coefs.:2)[2:(end - 1), :, 1, 1, 1] .≈
                -2 * K_ic / dz^2 * dψdϑ_ic,
            )
            @test all(
                parent(W.∂ϑₜ∂ϑ.coefs.:2)[end, :, 1, 1, 1] .≈
                -K_ic / dz^2 * dψdϑ_ic - K_ic / (dz * dz / 2) * dψdϑ_ic,
            )
        end

    end
end


@testset "Richards Jacobian entries, Flux BC" begin
    FT = Float64
    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(1.43)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0.124)
    zmax = FT(0)
    zmin = FT(-1.5)
    nelems = 150
    soil_domains = [
        Column(; zlim = (zmin, zmax), nelements = nelems),
        HybridBox(;
            xlim = (0.0, 1.0),
            ylim = (0.0, 1.0),
            zlim = (zmin, zmax),
            nelements = (1, 1, nelems),
            npolynomial = 3,
        ),
    ]
    top_flux_bc = FluxBC((p, t) -> eltype(t)(-K_sat))
    bot_flux_bc = FreeDrainage()
    sources = ()
    boundary_states =
        (; top = (water = top_flux_bc,), bottom = (water = bot_flux_bc,))
    params = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r)
    for domain in soil_domains
        soil = Soil.RichardsModel{FT}(;
            parameters = params,
            domain = domain,
            boundary_conditions = boundary_states,
            sources = sources,
        )

        Y, p, coords = initialize(soil)
        Y.soil.ϑ_l .= FT(0.24)
        W = RichardsTridiagonalW(Y)
        Wfact! = make_tendency_jacobian(soil)
        dtγ = FT(1.0)
        Wfact!(W, Y, p, dtγ, FT(0.0))

        K_ic = hydraulic_conductivity(
            K_sat,
            vg_m,
            effective_saturation(ν, FT(0.24), θ_r),
        )
        dz = FT(0.01)
        dψdϑ_ic = dψdϑ(FT(0.24), ν, θ_r, vg_α, vg_n, vg_m, S_s)
        if typeof(domain) <: Column
            @test parent(W.∂ϑₜ∂ϑ.coefs.:2)[1] .≈ -K_ic / dz^2 * dψdϑ_ic
            @test all(
                parent(W.∂ϑₜ∂ϑ.coefs.:2)[2:(end - 1)] .≈
                -2 * K_ic / dz^2 * dψdϑ_ic,
            )
            @test parent(W.∂ϑₜ∂ϑ.coefs.:2)[end] .≈ -K_ic / dz^2 * dψdϑ_ic
        elseif typeof(domain) <: HybridBox
            @test all(
                parent(W.∂ϑₜ∂ϑ.coefs.:2)[1, :, 1, 1, 1] .≈
                -K_ic / dz^2 * dψdϑ_ic,
            )
            @test all(
                parent(W.∂ϑₜ∂ϑ.coefs.:2)[2:(end - 1), :, 1, 1, 1] .≈
                -2 * K_ic / dz^2 * dψdϑ_ic,
            )
            @test all(
                parent(W.∂ϑₜ∂ϑ.coefs.:2)[end, :, 1, 1, 1] .≈
                -K_ic / dz^2 * dψdϑ_ic,
            )
        end

    end
end
