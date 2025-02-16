# This script tests the output of the ClimaLand TwoStream implementation
# against the T. Quaife pySellersTwoStream implementation by providing
# the same setups to each model and checking that the outputs are equal.
# The output of the python module for a variety of inputs is given in the
# linked file which then gets read and compared to the Clima version, ensuring
# the FAPAR of each model is within 1/2 of a percentage point.

using Test
using ClimaLand
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand.Canopy
using DelimitedFiles
using ClimaLand.Domains: Point

include("../../Artifacts.jl")

# Read the test data from the ClimaArtifact
datapath = twostr_test_data_path()

data = joinpath(datapath, "twostr_test.csv")
test_set = readdlm(data, ',')

# Floating point precision to use
for FT in (Float32, Float64)
    @testset "Two-Stream Model Correctness, FT = $FT" begin
        # Use a single band for the test
        spectral_discretization = TwoBandSpectralDiscretization{FT}()
        ε_canopy = FT(0.99)

        # Read the conditions for each setup parameter from the test file
        column_names = test_set[1, :]
        θs = acos.(FT.(test_set[2:end, column_names .== "mu"]))
        LAI = FT.(test_set[2:end, column_names .== "LAI"])
        a_soil_vals = FT.(test_set[2:end, column_names .== "a_soil"])
        n_layers = UInt64.(test_set[2:end, column_names .== "n_layers"])
        PropDif = FT.(test_set[2:end, column_names .== "prop_diffuse"])

        # setup spatially varying params as both float and spatially varying
        domain = Point(; z_sfc = FT(0.0))
        lds = FT.(test_set[2:end, column_names .== "ld"])
        lds_field = map(x -> fill(x, domain.space.surface), lds)
        a_soil_scalars = map(x -> FT.((x, x)), a_soil_vals)
        a_soil_fields =
            map(x -> fill(FT.((x, x)), domain.space.surface), a_soil_vals)
        a_soil_cases = (a_soil_scalars, a_soil_fields)
        ρ_PAR_leaf_scalar_vals = FT.(test_set[2:end, column_names .== "rho"])
        ρ_PAR_leaf_scalars = map(x -> (x, x), ρ_PAR_leaf_scalar_vals)
        ρ_PAR_leaf_fields =
            map(x -> fill((x, x), domain.space.surface), ρ_PAR_leaf_scalar_vals)
        τ_scalar_vals = FT.(test_set[2:end, column_names .== "tau"])
        τ_scalars = map(x -> (x, x), τ_scalar_vals)
        τ_fields = map(x -> fill((x, x), domain.space.surface), τ_scalar_vals)
        # loop through once with params as floats, then with params as fields
        Ω_cases = (FT(1), fill(FT(1), domain.space.surface))
        ρ_PAR_leaf_cases = (ρ_PAR_leaf_scalars, ρ_PAR_leaf_fields)
        τ_PAR_leaf_cases = (τ_scalars, τ_fields)
        lds_cases = (lds, lds_field)
        zipped_params = zip(
            Ω_cases,
            ρ_PAR_leaf_cases,
            τ_PAR_leaf_cases,
            lds_cases,
            a_soil_cases,
        )
        for (Ω, ρ_PAR_leaf, τ_PAR_leaf, lds, a_soil) in zipped_params
            # Read the result for each setup from the Python output
            py_FAPAR = FT.(test_set[2:end, column_names .== "FAPAR"])

            # Python code does not use clumping index, and λ_γ does not impact FAPAR
            # Test over all rows in the stored output from the Python module
            for i in 2:(size(test_set, 1) - 1)

                # Set the parameters based on the setup read from the file
                RT_params = TwoStreamParameters(;
                    spectral_discretization = spectral_discretization,
                    Ω = Ω,
                    G_Function = ConstantGFunction(FT.(lds[i])),
                    ρ_leaf = ρ_PAR_leaf[i],
                    τ_leaf = τ_PAR_leaf[i],
                    ϵ_canopy = FT(0.0), # Unused
                    λ_γ_PAR = FT(0.0),  # Unusedsp
                    n_layers = n_layers[i],
                )

                # Initialize the TwoStream model
                RT = TwoStreamModel(RT_params)

                # Compute the predicted FAPAR using the ClimaLand TwoStream implementation
                G = compute_G(RT_params.G_Function, θs)
                K = extinction_coeff.(G, θs[i])
                if typeof(RT_params.ρ_leaf) <: Tuple
                    output = canopy_sw_rt_two_stream(
                        G,
                        RT_params.Ω,
                        RT_params.n_layers,
                        RT_params.ρ_leaf,
                        RT_params.τ_leaf,
                        LAI[i],
                        K,
                        θs[i],
                        a_soil[i],
                        PropDif[i],
                    )
                else
                    output = fill(
                        (
                            (abs = FT(0.0), refl = FT(0.0), trans = FT(0.0)),
                            (abs = FT(0.0), refl = FT(0.0), trans = FT(0.0)),
                        ),
                        domain.space.surface,
                    )
                    output .=
                        canopy_sw_rt_two_stream.(
                            G,
                            RT_params.Ω,
                            RT_params.n_layers,
                            RT_params.ρ_leaf,
                            RT_params.τ_leaf,
                            LAI[i],
                            K,
                            θs[i],
                            a_soil[i],
                            PropDif[i],
                        )
                end
                if typeof(output) <: Tuple
                    FAPAR = output[1].abs
                else
                    get_abs = (x) -> x[1].abs
                    FAPAR = get_abs.(output)
                end
                # Check that the predictions are app. equivalent to the Python model
                # Create a field of the expect value because isapprox cannot be broadcast
                # over a field of floats. The domain is a point, so it makes no difference
                # to the error when FAPAR is a float
                expected_output = fill(py_FAPAR[i], domain.space.surface)
                @test isapprox(0, sum(FAPAR .- expected_output), atol = 0.005)
            end
        end
    end
end
