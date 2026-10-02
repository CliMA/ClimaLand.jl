using Test
import ClimaLand
import ClimaLand.Parameters as LP
import SurfaceFluxes
import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics

@testset "Stability cap in land surface fluxes" begin
    for FT in (Float32, Float64)
        toml_dict = LP.create_toml_dict(FT)
        earth_param_set = LP.LandParameters(toml_dict)
        sfp = LP.surface_fluxes_parameters(earth_param_set)
        @test sfp.ufp isa UF.GryanikParams
        rough = SurfaceFluxes.ConstantRoughnessParams{FT}(FT(1), FT(0.1))
        gust = SurfaceFluxes.ConstantGustinessSpec(FT(1))
        config = ClimaLand.surface_flux_config(rough, gust)
        @test config.stability_cap isa SurfaceFluxes.MaxHeatFluxStabilityCap

        # Strongly stable case (supercritical bulk Richardson number), in
        # which the capped exchange persists
        thermo_params = LP.thermodynamic_parameters(earth_param_set)
        T_a = FT(290)
        P = FT(1e5)
        q = FT(0.008)
        ρ = Thermodynamics.air_density(thermo_params, T_a, P, q)
        Δz = FT(30)
        solve(config, T_sfc) = SurfaceFluxes.surface_fluxes(
            sfp,
            T_a,
            q,
            FT(0),
            FT(0),
            ρ,
            T_sfc,
            q,
            FT(0),
            Δz,
            FT(0),
            (FT(1.5), FT(0)),
            (FT(0), FT(0)),
            nothing,
            config,
        )
        out = solve(config, T_a - FT(8))
        ζ_p = SurfaceFluxes.max_heat_flux_stability(sfp, Δz, FT(1))
        @test out.converged
        @test ζ_p < out.ζ < 100
        # u* equals the value at the cap
        κ = sfp.von_karman_const
        U = FT(1.5) # effective wind speed max(|u|, gustiness)
        F_m = UF.dimensionless_profile(
            sfp.ufp,
            Δz,
            ζ_p,
            FT(1),
            UF.MomentumTransport(),
        )
        @test out.ustar ≈ κ * U / F_m rtol = FT(0.01)
        @test out.shf < 0
        # Uncapped Monin-Obukhov exchange collapses
        out_uncapped =
            solve(SurfaceFluxes.SurfaceFluxConfig(rough, gust), T_a - FT(8))
        @test abs(out_uncapped.shf) < abs(out.shf) / 5
        # The downward heat flux increases with the surface-air temperature
        # difference
        shf = [solve(config, T_a - FT(ΔT)).shf for ΔT in (1, 2, 4, 8, 16)]
        @test all(diff(shf) .< 0)
    end
end
