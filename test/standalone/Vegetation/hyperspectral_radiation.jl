using Test
import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams as CP
using ClimaCore
using ClimaLand.Canopy

import ClimaLand

import ClimaLand.Parameters as LP

for FT in (Float32, Float64)
    @testset "Hyperspectral radiative transfer model, FT = $FT" begin

        LAI = FT(5.0) # m2 (leaf) m-2 (ground)

        # Set up radiative transfer using 16 band discretization
        spectral_discretization = HyperspectralDiscretization{FT}()
        ρ_leaf = FT.(ntuple(_ -> 0.0, 16))
        RTparams = BeerLambertParameters(FT;
            spectral_discretization,
            ρ_leaf,
        )

        # Drivers
        θs = FT.(Array(0:0.1:(π / 2)))
        G = compute_G(RTparams.G_Function, θs)
        K_c = extinction_coeff.(G, θs)
        α_soil = FT.(ntuple(_ -> 0.0, 16))

        # Test that radiative transfer works
        output = canopy_sw_rt_beer_lambert.(
            RTparams.Ω,
            (RTparams.ρ_leaf,),
            LAI,
            K_c,
            (α_soil,),
        )

        @test length(output) == length(θs)
        @test length(output[1]) == 16
        get_PAR_sum =
            (rt, sym, PAR_coeff) -> sum(map(x -> getproperty(x, sym), rt) .* PAR_coeff)
        @test all(get_PAR_sum.(output, :refl, RTparams.spectral_discretization.PAR_proportions) .< 1e-6)

    end
end
