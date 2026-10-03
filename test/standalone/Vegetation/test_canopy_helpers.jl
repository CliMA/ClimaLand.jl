using Test
import ClimaComms
ClimaComms.@import_required_backends
import StaticArrays
import SurfaceFluxes
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP

for FT in (Float32, Float64)
    @testset "Plant area index weighting, FT = $FT" begin
        LAI = FT(2)
        SAI = FT(1)
        @test Canopy.stem_area_fraction(LAI, SAI) ≈ FT(1 / 3)
        @test Canopy.stem_area_fraction(LAI, FT(0)) == FT(0)
        @test Canopy.stem_area_fraction(FT(0), FT(0)) == FT(0)
        @test Canopy.stem_area_fraction(FT(0), SAI) == FT(1)

        x_leaf = FT(0.1)
        x_stem = FT(0.4)
        @test Canopy.plant_area_weighted(x_leaf, x_stem, LAI, SAI) ≈
              (LAI * x_leaf + SAI * x_stem) / (LAI + SAI)
        @test Canopy.plant_area_weighted(x_leaf, x_stem, LAI, FT(0)) == x_leaf
        @test Canopy.plant_area_weighted(x_leaf, x_stem, FT(0), SAI) == x_stem
        @test Canopy.plant_area_weighted(x_leaf, x_stem, FT(0), FT(0)) == x_leaf

        @test Canopy.leaf_absorption_fraction(LAI, SAI) ≈ FT(2 / 3)
        @test Canopy.leaf_absorption_fraction(LAI, FT(0)) == FT(1)
        @test Canopy.leaf_absorption_fraction(FT(0), FT(0)) == FT(1)
        @test Canopy.leaf_absorption_fraction(FT(0), SAI) == FT(0)
        @test typeof(Canopy.leaf_absorption_fraction(LAI, SAI)) == FT
    end

    @testset "Canopy nitrogen scaling, FT = $FT" begin
        kn = FT(0.5)
        @test Canopy.canopy_nitrogen_scaling(kn, FT(0)) == FT(1)
        LAI = FT(4)
        expected = (1 - exp(-kn * LAI)) / (kn * LAI)
        @test Canopy.canopy_nitrogen_scaling(kn, LAI) ≈ expected
        @test 0 < Canopy.canopy_nitrogen_scaling(kn, LAI) < 1
        # Decreasing in LAI, and the canopy total LAI * scaling is increasing
        @test Canopy.canopy_nitrogen_scaling(kn, FT(1)) >
              Canopy.canopy_nitrogen_scaling(kn, FT(2))
        @test FT(2) * Canopy.canopy_nitrogen_scaling(kn, FT(2)) >
              FT(1) * Canopy.canopy_nitrogen_scaling(kn, FT(1))
        @test typeof(Canopy.canopy_nitrogen_scaling(kn, LAI)) == FT
    end

    @testset "Sub-canopy reference height and wind, FT = $FT" begin
        toml_dict = LP.create_toml_dict(FT)
        α = toml_dict["canopy_subcanopy_wind_extinction_coefficient"]
        z_min = toml_dict["canopy_subcanopy_min_reference_height"]
        @test α == FT(0.5)
        @test z_min == FT(2)
        h_atmos = FT(30)
        u = FT(3)
        g = FT(1)

        # Tall canopy: reference at the apparent sink height d + z_0m
        h = FT(20)
        d = FT(0.67) * h
        z_0m = FT(0.0823) * h + FT(0.01)
        z_ref = subcanopy_reference_height(FT(0), d, z_0m, z_min, h_atmos)
        @test z_ref == d + z_0m
        # Deep snow raises the floor; the forcing height is the ceiling
        @test subcanopy_reference_height(FT(14), d, z_0m, z_min, h_atmos) ==
              FT(16)
        @test subcanopy_reference_height(FT(29), d, z_0m, z_min, h_atmos) ==
              h_atmos
        # Short canopy: the floor binds
        @test subcanopy_reference_height(
            FT(0),
            FT(0.67) * FT(0.5),
            FT(0.0823) * FT(0.5) + FT(0.01),
            z_min,
            h_atmos,
        ) == z_min

        # At the forcing height with zero PAI, the wind is the effective wind
        # of the Monin-Obukhov solve; with positive PAI it is attenuated by
        # exp(-α * PAI)
        PAI = FT(4)
        @test subcanopy_wind(u, g, h_atmos, h_atmos, h, d, z_0m, FT(0), α) == u
        @test subcanopy_wind(
            FT(0.5),
            g,
            h_atmos,
            h_atmos,
            h,
            d,
            z_0m,
            FT(0),
            α,
        ) == g
        @test subcanopy_wind(u, g, h_atmos, h_atmos, h, d, z_0m, PAI, α) ≈
              u * exp(-α * PAI)
        # At or below the canopy top: the neutral log profile to max(z, h),
        # attenuated across the full plant area index PAI
        u_h = u * log((h - d) / z_0m) / log((h_atmos - d) / z_0m)
        @test subcanopy_wind(u, g, h_atmos, h, h, d, z_0m, FT(0), α) ≈ u_h
        @test subcanopy_wind(u, g, h_atmos, h, h, d, z_0m, PAI, α) ≈
              u_h * exp(-α * PAI)
        u_ref = subcanopy_wind(u, g, h_atmos, z_ref, h, d, z_0m, PAI, α)
        @test u_ref ≈ u_h * exp(-α * PAI)
        @test u_ref < subcanopy_wind(u, g, h_atmos, z_ref, h, d, z_0m, FT(1), α)
        @test typeof(u_ref) == FT
        # Short canopy (height < z_min): the log profile is evaluated at z_min
        # and still attenuated across the full PAI of the canopy below it
        h_s = FT(0.5)
        d_s = FT(0.67) * h_s
        z_0m_s = FT(0.0823) * h_s + FT(0.01)
        @test subcanopy_wind(u, g, h_atmos, z_min, h_s, d_s, z_0m_s, PAI, α) ≈
              u * log((z_min - d_s) / z_0m_s) / log((h_atmos - d_s) / z_0m_s) *
              exp(-α * PAI)
        # No canopy at all
        @test subcanopy_wind(
            u,
            g,
            h_atmos,
            z_min,
            FT(0),
            FT(0),
            FT(0.01),
            FT(0),
            α,
        ) ≈ u * log(z_min / FT(0.01)) / log(h_atmos / FT(0.01))

        # A wind vector keeps its direction
        uv = StaticArrays.SVector{2, FT}(3, 4)
        uv_ground = subcanopy_wind(uv, g, h_atmos, z_ref, h, d, z_0m, PAI, α)
        @test uv_ground isa StaticArrays.SVector{2, FT}
        @test hypot(uv_ground...) ≈
              subcanopy_wind(FT(5), g, h_atmos, z_ref, h, d, z_0m, PAI, α)
        @test uv_ground[1] / uv_ground[2] ≈ FT(3 / 4)
        # Calm air with gustiness still gives the attenuated gust speed
        calm = StaticArrays.SVector{2, FT}(0, 0)
        calm_ground =
            subcanopy_wind(calm, g, h_atmos, z_ref, h, d, z_0m, PAI, α)
        @test calm_ground[1] ≈
              subcanopy_wind(FT(0), g, h_atmos, z_ref, h, d, z_0m, PAI, α)
        @test calm_ground[2] == 0

        # The parameterization reads both parameters from the TOML
        sf = MoninObukhovCanopyFluxes(toml_dict, FT(10))
        @test sf.subcanopy_wind_extinction == α
        @test sf.subcanopy_min_reference_height == z_min
    end

    @testset "Ground gustiness below a canopy, FT = $FT" begin
        c = SurfaceFluxes.ConstantGustinessSpec(FT(2))
        # Below plants the floor is folded into the sub-canopy wind; bare
        # ground keeps the model of the forcing
        @test Base.materialize(ground_gustiness(c, true)) ==
              SurfaceFluxes.ConstantGustinessSpec(FT(0))
        @test Base.materialize(ground_gustiness(c, false)) == c
        @test Base.materialize(ground_gustiness(FT(2), false)) == c
        d = SurfaceFluxes.DeardorffGustinessSpec()
        @test ground_gustiness(d, true) === d
    end
end
