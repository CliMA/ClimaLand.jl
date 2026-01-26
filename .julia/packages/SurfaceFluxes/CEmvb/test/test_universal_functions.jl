using Test

import QuadGK

import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP

FT = Float32
param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
thermo_params = param_set.thermo_params

# TODO: Right now, we test these functions for
# type stability and correctness in the asymptotic
# limit. We may want to extend correctness tests.

@testset "UniversalFunctions" begin
    @testset "Type stability" begin
        FT = Float32
        ζ = FT(-2):FT(0.01):FT(200)
        for L in (-FT(10), FT(10))
            for ufp in (
                UF.GryanikParams(FT),
                UF.GrachevParams(FT),
                UF.BusingerParams(FT),
            )
                uft = UF.universal_func_type(typeof(ufp))
                uf = UF.universal_func(uft, L, ufp)
                for transport in (UF.MomentumTransport(), UF.HeatTransport())
                    ϕ = UF.phi.(uf, ζ, transport)
                    @test eltype(ϕ) == FT
                    ψ = UF.psi.(uf, ζ, transport)
                    @test eltype(ψ) == FT
                end
            end
        end

        # More type stability (phi/psi):
        FT = Float32
        ζ = (-FT(1), FT(0.5) * eps(FT), 2 * eps(FT))
        for L in (-FT(10), FT(10))
            for ufp in (
                UF.GryanikParams(FT),
                UF.GrachevParams(FT),
                UF.BusingerParams(FT),
            )
                uft = UF.universal_func_type(typeof(ufp))
                uf = UF.universal_func(uft, L, ufp)
                for transport in (UF.MomentumTransport(), UF.HeatTransport())
                    ϕ = UF.phi.(uf, ζ, transport)
                    @test eltype(ϕ) == FT
                    ψ = UF.psi.(uf, ζ, transport)
                    @test eltype(ψ) == FT
                end
            end
        end

        # More type stability (Psi):
        FT = Float32
        ζ = (-FT(1), -FT(0.5) * eps(FT), FT(0.5) * eps(FT), 2 * eps(FT))
        for L in (-FT(10), FT(10))
            for ufp in (UF.GryanikParams(FT), UF.BusingerParams(FT))
                uft = UF.universal_func_type(typeof(ufp))
                uf = UF.universal_func(uft, L, ufp)
                for transport in (UF.MomentumTransport(), UF.HeatTransport())
                    Ψ = UF.Psi.(uf, ζ, transport)
                    @test eltype(Ψ) == FT
                end
            end
        end

    end
    @testset "Asymptotic range" begin
        FT = Float32

        ϕ_h_ζ∞(uf::UF.Grachev, ζ) = 1 + FT(UF.b_h(uf))
        ϕ_m_ζ∞(uf::UF.Grachev, ζ) =
            FT(UF.a_m(uf)) / FT(UF.b_m(uf)) * ζ^FT(1 / 3)

        ϕ_h_ζ∞(uf::UF.Gryanik, ζ) =
            FT(1) +
            (FT(ζ) * FT(UF.Pr_0(uf)) * FT(UF.a_h(uf))) /
            (1 + FT(UF.b_h(uf)) * FT(ζ))
        ϕ_m_ζ∞(uf::UF.Gryanik, ζ) =
            FT(UF.a_m(uf) / UF.b_m(uf)^FT(2 / 3)) * ζ^FT(1 / 3)


        for L in (-FT(10), FT(10))
            for ufp in (UF.GryanikParams(FT), UF.GrachevParams(FT))
                uft = UF.universal_func_type(typeof(ufp))
                uf = UF.universal_func(uft, L, ufp)
                for ζ in FT(10) .^ (4, 6, 8, 10)
                    ϕ_h = UF.phi(uf, ζ, UF.HeatTransport())
                    @test isapprox(ϕ_h, ϕ_h_ζ∞(uf, ζ))
                end
                for ζ in FT(10) .^ (8, 9, 10)
                    ϕ_m = UF.phi(uf, ζ, UF.MomentumTransport())
                    @test isapprox(ϕ_m, ϕ_m_ζ∞(uf, ζ))
                end
            end
        end

    end

    # Test for Gryanik2021 Eq. 2 & 3; ensures Ψ(0) = 0
    @testset "Vanishes at Zero" begin
        FT = Float32
        for L in (-FT(10), FT(10))
            for ufp in (UF.GryanikParams(FT), UF.GrachevParams(FT))
                for transport in (UF.HeatTransport(), UF.MomentumTransport())
                    uft = UF.universal_func_type(typeof(ufp))
                    uf = UF.universal_func(uft, L, ufp)
                    Ψ_0 = UF.psi(uf, FT(0), transport)
                    @test isapprox(Ψ_0, FT(0))
                end
            end
        end
    end

    # Test that the integrated forms of the stability correction functions ψ(ζ) are consistent
    # when directly evaluated via numerical integration from ϕ(ζ).
    @testset "Test Correctness: ψ(ζ) = ∫ 𝒻(ϕ(ζ′)) dζ′" begin
        FloatType = (Float32, Float64)
        for FT in FloatType
            ζ_array = (
                FT(-20),
                FT(-10),
                FT(-1),
                -sqrt(eps(FT)),
                sqrt(eps(FT)),
                FT(1),
                FT(10),
                FT(20),
            )
            for L in FT(10) .* sign.(ζ_array)
                for ζ in ζ_array
                    for ufp in (
                        UF.GryanikParams(FT),
                        UF.GrachevParams(FT),
                        UF.BusingerParams(FT),
                    )
                        uft = UF.universal_func_type(typeof(ufp))
                        uf = UF.universal_func(uft, L, ufp)
                        for transport in
                            (UF.MomentumTransport(), UF.HeatTransport())
                            # Compute ψ via numerical integration of 𝒻(ϕ(ζ))
                            ψ_int = QuadGK.quadgk(
                                ζ′ ->
                                    (FT(1) - UF.phi(uf, ζ′, transport)) / ζ′,
                                eps(FT),
                                ζ,
                            )
                            # Compute ψ using function definitions of ψ(ζ)
                            ψ = UF.psi(uf, ζ, transport)
                            @test isapprox(
                                ψ_int[1] - ψ,
                                FT(0),
                                atol = 10sqrt(eps(FT)),
                            )
                        end
                    end
                end
            end
        end
    end

end
