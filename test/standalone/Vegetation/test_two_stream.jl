# This script tests the canopy absorbed output of the ClimaLand TwoStream implementation
# against the T. Quaife pySellersTwoStream implementation by providing
# the same setups to each model and checking that the outputs are equal.
# The output of the python module for a variety of inputs is given in the
# linked file which then gets read and compared to the Clima version, ensuring
# the FAPAR of each model is within 1/2 of a percentage point.

# It also tests that the sum of absorbed and reflected radiation is 1.

using Test
using ClimaLand
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand.Canopy
using DelimitedFiles
using ClimaLand.Domains: Point
import ClimaLand.Parameters as LP
import ClimaParams
using ClimaCore

@testset "Comparison to pySellersTwoStream" begin
    include("../../Artifacts.jl")

    # Read the test data from the ClimaArtifact
    datapath = twostr_test_data_path()

    data = joinpath(datapath, "twostr_test.csv")
    test_set = readdlm(data, ',')

    # Floating point precision to use
    for FT in (Float32, Float64)
        @testset "Two-Stream Model Correctness, FT = $FT" begin
            # Read the conditions for each setup parameter from the test file
            column_names = test_set[1, :]
            cosθs = FT.(test_set[2:end, column_names .== "mu"])
            LAI = FT.(test_set[2:end, column_names .== "LAI"])
            a_soil = FT.(test_set[2:end, column_names .== "a_soil"])
            n_layers = UInt64.(test_set[2:end, column_names .== "n_layers"])
            PropDif = FT.(test_set[2:end, column_names .== "prop_diffuse"])

            # setup spatially varying params as both float and spatially varying
            domain = Point(; z_sfc = FT(0.0))
            lds = FT.(test_set[2:end, column_names .== "ld"])
            lds_field = map(x -> fill(x, domain.space.surface), lds)
            α_PAR_leaf_scalars = FT.(test_set[2:end, column_names .== "rho"])
            α_PAR_leaf_fields =
                map(x -> fill(x, domain.space.surface), α_PAR_leaf_scalars)
            τ_scalars = FT.(test_set[2:end, column_names .== "tau"])
            τ_fields = map(x -> fill(x, domain.space.surface), τ_scalars)
            # loop through once with params as floats, then with params as fields
            Ω_cases = (FT(1), fill(FT(1), domain.space.surface))
            α_PAR_leaf_cases = (α_PAR_leaf_scalars, α_PAR_leaf_fields)
            τ_PAR_leaf_cases = (τ_scalars, τ_fields)
            α_NIR_leaf_cases = (FT(0.4), fill(FT(0.4), domain.space.surface))
            τ_NIR_leaf_cases = (FT(0.25), fill(FT(0.24), domain.space.surface))
            lds_cases = (lds, lds_field)
            zipped_params = zip(
                Ω_cases,
                α_PAR_leaf_cases,
                τ_PAR_leaf_cases,
                α_NIR_leaf_cases,
                τ_NIR_leaf_cases,
                lds_cases,
            )
            output = ClimaCore.Fields.zeros(
                @NamedTuple{abs::FT, refl::FT, trans::FT},
                domain.space.surface,
            )
            for (Ω, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf, lds) in
                zipped_params
                # Read the result for each setup from the Python output
                py_FAPAR = FT.(test_set[2:end, column_names .== "FAPAR"])

                # Python code does not use clumping index, and λ_γ does not impact FAPAR
                # Test over all rows in the stored output from the Python module
                for i in 2:(size(test_set, 1) - 1)

                    # Set the parameters based on the setup read from the file
                    RT_params = TwoStreamParameters(
                        FT;
                        Ω = Ω,
                        G_Function = ConstantGFunction.(FT.(lds[i])),
                        α_PAR_leaf = α_PAR_leaf[i],
                        τ_PAR_leaf = τ_PAR_leaf[i],
                        α_NIR_leaf = α_NIR_leaf,
                        τ_NIR_leaf = τ_NIR_leaf,
                        n_layers = n_layers[i],
                    )

                    # Initialize the TwoStream model
                    RT = TwoStreamModel(RT_params)

                    # Compute the predicted FAPAR using the ClimaLand TwoStream implementation
                    @. output = canopy_sw_rt_two_stream(
                        RT_params.G_Function,
                        RT_params.Ω,
                        RT_params.n_layers,
                        RT_params.α_PAR_leaf,
                        RT_params.τ_PAR_leaf,
                        LAI[i],
                        cosθs[i],
                        a_soil[i],
                        PropDif[i],
                    )
                    FAPAR = output.abs
                    # Check that the predictions are app. equivalent to the Python model
                    # Create a field of the expect value because isapprox cannot be broadcast
                    # over a field of floats. The domain is a point, so it makes no difference
                    # to the error when FAPAR is a float
                    expected_output = fill(py_FAPAR[i], domain.space.surface)
                    @test isapprox(
                        0,
                        sum(FAPAR .- expected_output),
                        atol = 0.005,
                    )
                end
            end
        end
    end
end

@testset "Test physicality" begin
    N = 100000
    θs = [rand(N - 1) * 2π..., π / 2]
    cosθs = cos.(θs)
    α_leaf = [rand(N - 4)..., 0.0, 0.0, 1.0, 1.0]
    τ_leaf = (1.0 .- α_leaf) .* rand(N)
    α_soil = [rand(N - 6)..., 0.0, 1.0, 0.2, 0.2, 0.2, 0.2]
    G_Function = ConstantGFunction(0.5)
    K = ClimaLand.Canopy.extinction_coeff.(G_Function, cosθs)
    frac_diff = rand(N)
    n_layers = UInt64(20)
    Ω = rand(N)
    LAI = round.(rand(N))
    output =
        ClimaLand.Canopy.canopy_sw_rt_two_stream.(
            G_Function,
            Ω,
            n_layers,
            α_leaf,
            τ_leaf,
            LAI,
            cosθs,
            α_soil,
            frac_diff,
        )
    expected = zeros(N)
    for i in 1:N
        expected[i] =
            output[i].trans * (1 - α_soil[i]) + output[i].abs + output[i].refl
    end
    @assert all(expected .≈ 1)
end
