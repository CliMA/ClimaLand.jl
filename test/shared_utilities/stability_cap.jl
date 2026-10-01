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

        # Strongly stable case (supercritical bulk Richardson number):
        # the exchange no longer collapses
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
        out_uncapped = solve(SurfaceFluxes.SurfaceFluxConfig(rough, gust), T_a - FT(8))
        @test abs(out_uncapped.shf) < abs(out.shf) / 5
        # Downward heat flux increases with the surface-air temperature difference
        shf = [solve(config, T_a - FT(ΔT)).shf for ΔT in (1, 2, 4, 8, 16)]
        @test all(diff(shf) .< 0)

        @testset "Undercanopy resistance and above-canopy stability, FT = $FT" begin
            PAI = FT(4.0)
            h_c = FT(18.0)
            z_0m_c = FT(1.5)
            z_0b_c = FT(1.35)
            displ_c = FT(12.0)
            leaf_Cd = FT(0.07)
            Δz_ref = FT(30.0)
            u_air = FT(3.0)
            gustiness = FT(1.0)
            C_s = FT(0.004)
            grav = LP.grav(earth_param_set)
            cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
            T_air = FT(290.0)
            T_a_uc = T_air + grav * Δz_ref / cp_d

            # Neutral above and under canopy (T_canopy == T_a_uc == T_ground)
            r_neutral = ClimaLand.undercanopy_resistance_at_a_point(
                PAI,
                h_c,
                z_0m_c,
                z_0b_c,
                displ_c,
                leaf_Cd,
                T_a_uc,
                T_a_uc,
                T_air,
                u_air,
                Δz_ref,
                gustiness,
                C_s,
                earth_param_set,
            )
            z_eff = Δz_ref - displ_c
            u_star_0 = κ * u_air / log(z_eff / z_0m_c)
            @test r_neutral ≈ 1 / (C_s * u_star_0)

            # Stable above canopy (T_canopy < T_a_uc, T_ground == T_a_uc so under-canopy Ri == 0):
            # u*_c < u*_0 => r_stable > r_neutral, bounded by max_heat_flux_stability cap
            r_stable = ClimaLand.undercanopy_resistance_at_a_point(
                PAI,
                h_c,
                z_0m_c,
                z_0b_c,
                displ_c,
                leaf_Cd,
                T_a_uc - FT(3),
                T_a_uc,
                T_air,
                u_air,
                Δz_ref,
                gustiness,
                C_s,
                earth_param_set,
            )
            r_very_stable = ClimaLand.undercanopy_resistance_at_a_point(
                PAI,
                h_c,
                z_0m_c,
                z_0b_c,
                displ_c,
                leaf_Cd,
                T_a_uc - FT(20),
                T_a_uc,
                T_air,
                u_air,
                Δz_ref,
                gustiness,
                C_s,
                earth_param_set,
            )
            @test r_stable > r_neutral
            @test r_very_stable ≈ r_stable rtol = FT(0.05)
            @test FT(1.4) * r_neutral < r_very_stable < FT(1.6) * r_neutral

            # Unstable above canopy (T_canopy > T_a_uc, T_ground == T_canopy so under-canopy Ri <= 0):
            # u*_c > u*_0 => r_unstable < r_neutral
            r_unstable = ClimaLand.undercanopy_resistance_at_a_point(
                PAI,
                h_c,
                z_0m_c,
                z_0b_c,
                displ_c,
                leaf_Cd,
                T_a_uc + FT(4),
                T_a_uc + FT(4),
                T_air,
                u_air,
                Δz_ref,
                gustiness,
                C_s,
                earth_param_set,
            )
            @test r_unstable < r_neutral

            # Zero PAI: T_af == T_a_uc regardless of T_canopy, recovering neutral u*_c
            r_zero_pai = ClimaLand.undercanopy_resistance_at_a_point(
                FT(0),
                h_c,
                z_0m_c,
                z_0b_c,
                displ_c,
                leaf_Cd,
                T_a_uc - FT(10),
                T_a_uc,
                T_air,
                u_air,
                Δz_ref,
                gustiness,
                C_s,
                earth_param_set,
            )
            @test r_zero_pai ≈ r_neutral
        end
    end
end

