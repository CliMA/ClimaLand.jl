using Test
import ClimaComms
ClimaComms.@import_required_backends
import StaticArrays
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

    @testset "Sub-canopy wind, FT = $FT" begin
        toml_dict = LP.create_toml_dict(FT)
        α = toml_dict["canopy_subcanopy_wind_extinction_coefficient"]
        @test α == FT(0.5)
        u = FT(3)
        g = FT(1)
        # No canopy: the effective wind speed of the Monin-Obukhov solve
        @test subcanopy_wind(u, g, FT(0), α) == u
        @test subcanopy_wind(FT(0.5), g, FT(0), α) == g
        @test subcanopy_wind(u, FT(0), FT(0), α) == u
        PAI = FT(4)
        @test subcanopy_wind(u, g, PAI, α) ≈ exp(-α * PAI) * u
        @test subcanopy_wind(u, g, PAI, α) < subcanopy_wind(u, g, FT(1), α)
        @test typeof(subcanopy_wind(u, g, PAI, α)) == FT

        # A wind vector keeps its direction
        uv = StaticArrays.SVector{2, FT}(3, 4)
        uv_ground = subcanopy_wind(uv, g, PAI, α)
        @test uv_ground isa StaticArrays.SVector{2, FT}
        @test hypot(uv_ground...) ≈ subcanopy_wind(FT(5), g, PAI, α)
        @test uv_ground[1] / uv_ground[2] ≈ FT(3 / 4)
        # Calm air with gustiness still gives the attenuated gust speed
        calm = StaticArrays.SVector{2, FT}(0, 0)
        @test subcanopy_wind(calm, g, PAI, α) ==
              StaticArrays.SVector{2, FT}(exp(-α * PAI) * g, 0)

        # The parameterization reads the extinction coefficient from the TOML
        sf = MoninObukhovCanopyFluxes(toml_dict, FT(10))
        @test sf.subcanopy_wind_extinction == α
    end
end
