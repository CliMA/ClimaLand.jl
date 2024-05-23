using Test
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
using LinearAlgebra
import ClimaCore: MatrixFields
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains: Column, HybridBox
using ClimaLand.Soil

import ClimaLand
import ClimaLand.Parameters as LP


for FT in (Float32, Float64)
    @testset "Richards Jacobian entries, Moisture BC, FT = $FT" begin
        ν = FT(0.495)
        K_sat = FT(0.0443 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(1.43)
        vg_α = FT(2.6) # inverse meters
        hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
        θ_r = FT(0.124)
        zmax = FT(0)
        zmin = FT(-1.5)
        nelems = 150
        soil_domains = [
            Column(; zlim = (zmin, zmax), nelements = nelems),
            HybridBox(;
                xlim = FT.((0, 1)),
                ylim = FT.((0, 1)),
                zlim = (zmin, zmax),
                nelements = (1, 1, nelems),
                npolynomial = 3,
            ),
        ]
        top_state_bc = MoistureStateBC((p, t) -> ν - 1e-3)
        bot_flux_bc = FreeDrainage()
        sources = ()
        boundary_states = (; top = top_state_bc, bottom = bot_flux_bc)
        params = Soil.RichardsParameters(ν, hcm, K_sat, S_s, θ_r)

        for domain in soil_domains
            soil = Soil.RichardsModel{FT}(;
                parameters = params,
                domain = domain,
                boundary_conditions = boundary_states,
                sources = sources,
            )

            Y, p, coords = initialize(soil)
            Y.soil.ϑ_l .= FT(0.24)
            # We do not set the initial aux state here because
            # we want to test that it is updated correctly in
            # the jacobian correctly.
            jacobian = ImplicitEquationJacobian(Y)
            jac_tendency! = make_tendency_jacobian(soil)
            dtγ = FT(1.0)
            jac_tendency!(jacobian, Y, p, dtγ, FT(0.0))

            K_ic = hydraulic_conductivity(
                hcm,
                K_sat,
                effective_saturation(ν, FT(0.24), θ_r),
            )
            dz = FT(0.01)
            dψdϑ_ic = dψdϑ(hcm, FT(0.24), ν, θ_r, S_s)

            @test jacobian.solver isa MatrixFields.FieldMatrixSolver
            @test jacobian.solver.alg isa MatrixFields.BlockDiagonalSolve
            @test jacobian.matrix.keys.values ==
                  ((MatrixFields.@name(soil.ϑ_l), MatrixFields.@name(soil.ϑ_l)),)

            jac_ϑ_l = jacobian.matrix[
                MatrixFields.@name(soil.ϑ_l),
                MatrixFields.@name(soil.ϑ_l)
            ]
            if typeof(domain) <: Column
                # Check diagonals on either side of the main diagonal
                @test all(
                    parent(jac_ϑ_l.entries.:1)[2:end] .≈
                    parent(jac_ϑ_l.entries.:3)[1:(end - 1)],
                )
                @test Array(parent(jac_ϑ_l.entries.:1))[1] == FT(0)
                @test Array(parent(jac_ϑ_l.entries.:3))[end] == FT(0)
                @test all(
                    Array(parent(jac_ϑ_l.entries.:1))[2:end] .≈
                    dtγ * (K_ic / dz^2 * dψdϑ_ic),
                )
                # Check values on main diagonal: note jac_ϑ_l is I-γJ
                @test Array(parent(jac_ϑ_l.entries.:2))[1] .≈
                      dtγ * (-K_ic / dz^2 * dψdϑ_ic) - I
                @test all(
                    Array(parent(jac_ϑ_l.entries.:2))[2:(end - 1)] .≈
                    dtγ * (-2 * K_ic / dz^2 * dψdϑ_ic) - I,
                )
                @test Array(parent(jac_ϑ_l.entries.:2))[end] .≈
                      dtγ * (
                    -K_ic / dz^2 * dψdϑ_ic - K_ic / (dz * dz / 2) * dψdϑ_ic
                ) - I
            elseif typeof(domain) <: HybridBox
                # Check diagonals on either side of the main diagonal
                @test all(
                    Array(parent(jac_ϑ_l.entries.:1))[2:end, :, 1, 1, 1] .≈
                    Array(parent(jac_ϑ_l.entries.:3))[1:(end - 1), :, 1, 1, 1],
                )
                @test all(
                    Array(parent(jac_ϑ_l.entries.:1))[1, :, 1, 1, 1] .== FT(0),
                )
                @test all(
                    Array(parent(jac_ϑ_l.entries.:3))[end, :, 1, 1, 1] .==
                    FT(0),
                )
                @test all(
                    Array(parent(jac_ϑ_l.entries.:1))[2:end, :, 1, 1, 1] .≈
                    dtγ * (K_ic / dz^2 * dψdϑ_ic),
                )
                # Check values on main diagonal: note jac_ϑ_l is I-γJ
                @test all(
                    Array(parent(jac_ϑ_l.entries.:2))[1, :, 1, 1, 1] .≈
                    dtγ * (-K_ic / dz^2 * dψdϑ_ic) - I,
                )
                @test all(
                    Array(parent(jac_ϑ_l.entries.:2))[
                        2:(end - 1),
                        :,
                        1,
                        1,
                        1,
                    ] .≈ dtγ * (-2 * K_ic / dz^2 * dψdϑ_ic) - I,
                )
                @test all(
                    Array(parent(jac_ϑ_l.entries.:2))[end, :, 1, 1, 1] .≈
                    dtγ *
                    (-K_ic / dz^2 * dψdϑ_ic - K_ic / (dz * dz / 2) * dψdϑ_ic) -
                    I,
                )
            end

        end
    end

    @testset "Richards Jacobian entries, Flux BC, FT = $FT" begin
        ν = FT(0.495)
        K_sat = FT(0.0443 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(1.43)
        vg_α = FT(2.6) # inverse meters
        hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
        θ_r = FT(0.124)
        zmax = FT(0)
        zmin = FT(-1.5)
        nelems = 150
        soil_domains = [
            Column(; zlim = (zmin, zmax), nelements = nelems),
            HybridBox(;
                xlim = FT.((0, 1)),
                ylim = FT.((0, 1)),
                zlim = (zmin, zmax),
                nelements = (1, 1, nelems),
                npolynomial = 3,
            ),
        ]
        top_flux_bc = WaterFluxBC((p, t) -> -K_sat)
        bot_flux_bc = FreeDrainage()
        sources = ()
        boundary_states = (; top = top_flux_bc, bottom = bot_flux_bc)
        params = Soil.RichardsParameters(ν, hcm, K_sat, S_s, θ_r)
        for domain in soil_domains
            soil = Soil.RichardsModel{FT}(;
                parameters = params,
                domain = domain,
                boundary_conditions = boundary_states,
                sources = sources,
            )

            Y, p, coords = initialize(soil)
            Y.soil.ϑ_l .= FT(0.24)
            # We do not set the initial aux state here because
            # we want to test that it is updated correctly in
            # the jacobian correctly.
            jacobian = ImplicitEquationJacobian(Y)
            jacobian_tendency! = make_tendency_jacobian(soil)
            dtγ = FT(1.0)
            jacobian_tendency!(jacobian, Y, p, dtγ, FT(0.0))

            K_ic = hydraulic_conductivity(
                hcm,
                K_sat,
                effective_saturation(ν, FT(0.24), θ_r),
            )
            dz = FT(0.01)
            dψdϑ_ic = dψdϑ(hcm, FT(0.24), ν, θ_r, S_s)
            jac_ϑ_l = jacobian.matrix[
                MatrixFields.@name(soil.ϑ_l),
                MatrixFields.@name(soil.ϑ_l)
            ]
            if typeof(domain) <: Column
                @test Array(parent(jac_ϑ_l.entries.:2))[1] .≈
                      dtγ * (-K_ic / dz^2 * dψdϑ_ic) - I
                @test all(
                    Array(parent(jac_ϑ_l.entries.:2))[2:(end - 1)] .≈
                    dtγ * (-2 * K_ic / dz^2 * dψdϑ_ic) - I,
                )
                @test Array(parent(jac_ϑ_l.entries.:2))[end] .≈
                      dtγ * (-K_ic / dz^2 * dψdϑ_ic) - I
            elseif typeof(domain) <: HybridBox
                @test all(
                    Array(parent(jac_ϑ_l.entries.:2))[1, :, 1, 1, 1] .≈
                    dtγ * (-K_ic / dz^2 * dψdϑ_ic) - I,
                )
                @test all(
                    Array(parent(jac_ϑ_l.entries.:2))[
                        2:(end - 1),
                        :,
                        1,
                        1,
                        1,
                    ] .≈ dtγ * (-2 * K_ic / dz^2 * dψdϑ_ic) - I,
                )
                @test all(
                    Array(parent(jac_ϑ_l.entries.:2))[end, :, 1, 1, 1] .≈
                    dtγ * (-K_ic / dz^2 * dψdϑ_ic) - I,
                )
            end
        end
    end
end
