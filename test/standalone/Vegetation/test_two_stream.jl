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
using ArtifactWrappers
using ClimaLand.Domains: Point
import ClimaLand.Parameters as LP
import ClimaParams
# Get the test data
data_file = ArtifactWrapper(
    @__DIR__,
    "PySellersTwoStream Data",
    ArtifactFile[ArtifactFile(
        url = "https://caltech.box.com/shared/static/e7angzdnw18tmf8gctse5flsrkjsrlhx.csv",
        filename = "2_str_test_data.csv",
    ),],
)
datapath = get_data_folder(data_file)
data = joinpath(datapath, "2_str_test_data.csv")
test_set = readdlm(data, ',')

# Floating point precision to use
for FT in (Float32, Float64)
    @testset "Two-Stream Model Correctness, FT = $FT" begin
        # Read the conditions for each setup parameter from the test file
        column_names = test_set[1, :]
        θs = acos.(FT.(test_set[2:end, column_names .== "mu"]))
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
        α_PAR_leaf_cases = (α_PAR_leaf_scalars, α_PAR_leaf_fields)
        τ_PAR_leaf_cases = (τ_scalars, τ_fields)
        α_NIR_leaf_cases = (FT(0.4), fill(FT(0.4), domain.space.surface))
        τ_NIR_leaf_cases = (FT(0.25), fill(FT(0.24), domain.space.surface))
        lds_cases = (lds, lds_field)
        zipped_params = zip(
            α_PAR_leaf_cases,
            τ_PAR_leaf_cases,
            α_NIR_leaf_cases,
            τ_NIR_leaf_cases,
            lds_cases,
        )
        for (α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf, lds) in
            zipped_params
            # Read the result for each setup from the Python output
            py_FAPAR = FT.(test_set[2:end, column_names .== "FAPAR"])

            # Python code does not use clumping index, and λ_γ does not impact FAPAR
            # Test over all rows in the stored output from the Python module
            for i in 2:(size(test_set, 1) - 1)

                # Set the parameters based on the setup read from the file
                RT_params = TwoStreamParameters(
                    FT;
                    G_Function = ConstantGFunction(FT.(lds[i])),
                    α_PAR_leaf = α_PAR_leaf[i],
                    τ_PAR_leaf = τ_PAR_leaf[i],
                    α_NIR_leaf = α_NIR_leaf,
                    τ_NIR_leaf = τ_NIR_leaf,
                    Ω = FT(1),
                    n_layers = n_layers[i],
                )

                # Initialize the TwoStream model
                RT = TwoStreamModel(RT_params)

                # Compute the predicted FAPAR using the ClimaLand TwoStream implementation
                G = compute_G(RT_params.G_Function, θs)
                K = extinction_coeff.(G, θs[i])
                output =
                    canopy_sw_rt_two_stream.(
                        G,
                        RT_params.Ω,
                        RT_params.n_layers,
                        RT_params.α_PAR_leaf,
                        RT_params.τ_PAR_leaf,
                        LAI[i],
                        K,
                        θs[i],
                        a_soil[i],
                        PropDif[i],
                    )
                FAPAR = output.abs
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
