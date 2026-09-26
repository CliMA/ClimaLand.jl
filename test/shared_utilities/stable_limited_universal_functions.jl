using Test
import ClimaLand
import ClimaLand.Parameters as LP
import SurfaceFluxes
import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics

@testset "Stable-limited universal functions" begin
    for FT in (Float32, Float64)
        toml_dict = LP.create_toml_dict(FT)
        earth_param_set = LP.LandParameters(toml_dict)
        sfp = LP.surface_fluxes_parameters(earth_param_set)
        uf = sfp.ufp
        @test uf isa LP.StableLimitedUniversalFunctionParams{FT}
        ζc = uf.ζ_max_stable
        @test ζc == toml_dict["zeta_max_stable"]
        base = uf.base
        @test UF.Pr_0(uf) == UF.Pr_0(base)
        Δz = FT(30)
        z0 = FT(0.5)
        for scheme in (UF.PointValueScheme(), UF.LayerAverageScheme()),
            tt in (UF.MomentumTransport(), UF.HeatTransport())

            # Unchanged for unstable and moderately stable conditions
            for ζ in FT.((-10, -1, -0.1, 0, 0.1, ζc))
                @test UF.dimensionless_profile(uf, Δz, ζ, z0, tt, scheme) ==
                      UF.dimensionless_profile(base, Δz, ζ, z0, tt, scheme)
            end
            # Independent of ζ beyond the limit
            F_c = UF.dimensionless_profile(base, Δz, ζc, z0, tt, scheme)
            for ζ in FT.((2ζc, 10, 100))
                @test UF.dimensionless_profile(uf, Δz, ζ, z0, tt, scheme) ==
                      F_c
            end
            @test UF.phi(uf, FT(3), tt) == UF.phi(base, FT(3), tt)
            @test UF.psi(uf, FT(3), tt) == UF.psi(base, FT(3), tt)
        end
        # Bulk Richardson number increases monotonically with ζ, so the
        # Monin-Obukhov solve has a root for any stable Ri_b
        ζs = FT.(range(0, 100; length = 400))
        Ri = [UF.bulk_richardson_number(uf, Δz, ζ, z0, z0 / 10) for ζ in ζs]
        @test all(diff(Ri) .> 0)

        # Strongly stable case: exchange no longer collapses
        thermo_params = LP.thermodynamic_parameters(earth_param_set)
        T_a = FT(290)
        P = FT(1e5)
        q = FT(0.008)
        ρ = Thermodynamics.air_density(thermo_params, T_a, P, q)
        rough = SurfaceFluxes.ConstantRoughnessParams{FT}(FT(1), FT(0.1))
        config = SurfaceFluxes.SurfaceFluxConfig(
            rough,
            SurfaceFluxes.ConstantGustinessSpec(FT(1)),
        )
        out = SurfaceFluxes.surface_fluxes(
            sfp,
            T_a,
            q,
            FT(0),
            FT(0),
            ρ,
            T_a - FT(8),
            q,
            FT(0),
            FT(30),
            FT(0),
            (FT(1.5), FT(0)),
            (FT(0), FT(0)),
            nothing,
            config,
        )
        @test out.ζ > ζc
        @test out.ζ < 100 # a root is found below SurfaceFluxes' ζ limit
        # u* equals the neutral-profile value evaluated at ζ = ζc
        κ = sfp.von_karman_const
        U = FT(1.5) # effective wind speed max(|u|, gustiness)
        F_m = UF.dimensionless_profile(
            base,
            FT(30),
            ζc,
            FT(1),
            UF.MomentumTransport(),
        )
        @test out.ustar ≈ κ * U / F_m rtol = FT(0.05)
        @test out.shf < 0
    end
end
