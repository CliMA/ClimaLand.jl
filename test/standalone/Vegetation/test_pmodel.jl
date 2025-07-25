"""
# P-Model Regression Tests

This module tests the ClimaLand P-Model implementation against reference outputs 
from the R-based P-Model package. The P-Model predicts photosynthetic carbon 
assimilation and related variables based on environmental conditions.

## Test Structure

The tests load input and expected output data from CSV files provided by the
`pmodel_unittests` artifact. Each test case includes:

### Input Variables:
- `tc`: Temperature in Celsius
- `patm`: Atmospheric pressure in Pa
- `vpd`: Vapor pressure deficit
- `co2`: CO2 concentration in ppm
- `fapar`: Fraction of absorbed photosynthetically active radiation
- `ppfd`: Photosynthetic photon flux density
- `soilm`: Soil moisture (when soil moisture stress is enabled)
- `beta`: Beta parameter for CO2 limitation
- `kphio`: intrinsic quantum yield parameter
- `do_ftemp_kphio`: Flag for temperature-dependent intrinsic quantum yield
- `do_soilmstress`: Flag for soil moisture stress

### Expected Outputs:
- `gpp`: Gross primary productivity
- Additional photosynthetic variables as computed by the P-Model

## Validation Criteria

Tests use relative tolerance of 1e-3 and absolute tolerance of 1e-3 to account
for numerical differences between implementations while ensuring scientific accuracy.
"""

using Test
using ClimaLand
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand.Canopy
using DelimitedFiles
using ClimaLand.Domains: Point
import ClimaLand.Parameters as LP
import ClimaParams


"""
    percent_difference(a, b)

Calculate the percentage difference between `a` and `b` (reference). 
"""
function percent_difference(a, b)
    return abs(a - b) / abs(b) * 100
end

"""
    create_pmodel_parameters(inputs::Dict{String, Any}, FT)

Constructs the parameter structure for the P-Model
"""
function PModelParameters(inputs::Dict{String, Any}, FT)
    # these are default values used in Stocker 2020
    β = FT(inputs["beta"])
    cstar = FT(0.41)

    # handle temperature-dependent quantum yield 
    if Bool(inputs["do_ftemp_kphio"])
        ϕ0 = FT(NaN)
        ϕc = FT(inputs["kphio"])
        # Eqn 20 from Stocker 2020
        ϕa0 = FT(0.352)
        ϕa1 = FT(0.022)
        ϕa2 = FT(-0.00034)
    else
        ϕ0 = FT(inputs["kphio"])
        ϕc, ϕa0, ϕa1, ϕa2 = FT(NaN), FT(NaN), FT(NaN), FT(NaN)
    end

    return ClimaLand.Canopy.PModelParameters(
        cstar = cstar,
        β = β,
        ϕc = ϕc,
        ϕ0 = ϕ0,
        ϕa0 = ϕa0,
        ϕa1 = ϕa1,
        ϕa2 = ϕa2,
    )
end


"""
    csv_to_dict(data, header, row_idx)

Convert a CSV row to a dictionary, skipping missing values and the testcase column.
"""
function csv_to_dict(data, header, row_idx)
    dict = Dict{String, Any}()
    for (col_idx, col_name) in enumerate(header)
        col_name_str = string(col_name)
        value = data[row_idx, col_idx]

        # Skip testcase column and missing values
        if col_name_str != "testcase" && !ismissing(value)
            dict[col_name_str] = value
        end
    end
    return dict
end

@testset "P-model regression tests against R output" begin
    # verbose output?
    verbose = false

    # Directory containing CSV test cases
    include("../../Artifacts.jl")
    datadir = pmodel_unittests_path()
    inputs_file = joinpath(datadir, "inputs.csv")
    outputs_file = joinpath(datadir, "outputs.csv")

    # tolerances
    atol = 1e-3
    rtol = 1e-3

    # read inputs and outputs from CSV files
    inputs_data, inputs_header = readdlm(inputs_file, ',', header = true)
    outputs_data, outputs_header = readdlm(outputs_file, ',', header = true)
    inputs_header = vec(inputs_header)
    outputs_header = vec(outputs_header)

    # Get testcase column indices
    inputs_testcase_idx = findfirst(h -> string(h) == "testcase", inputs_header)
    outputs_testcase_idx =
        findfirst(h -> string(h) == "testcase", outputs_header)

    # Process each test case
    for row_idx in 1:size(inputs_data, 1)
        testcase_name = inputs_data[row_idx, inputs_testcase_idx]

        # Find corresponding output row
        output_row_idx = findfirst(
            r -> outputs_data[r, outputs_testcase_idx] == testcase_name,
            1:size(outputs_data, 1),
        )

        if output_row_idx === nothing
            @warn "No output data found for testcase: $testcase_name"
            continue
        end

        # Create inputs and outputs dictionaries using helper function
        inputs = csv_to_dict(inputs_data, inputs_header, row_idx)
        ref_outputs = csv_to_dict(outputs_data, outputs_header, output_row_idx)

        for FT in (Float32, Float64)
            # Convert ref_outputs to the correct FT
            ref_outputs_typed =
                Dict{String, FT}(k => FT(v) for (k, v) in ref_outputs)

            @testset "Test $testcase_name, FT = $FT" begin
                verbose &&
                    println("Running test case: $testcase_name with FT = $FT")

                # prepare constants, parameters, and drivers for the current FT
                constants = PModelConstants{FT}()
                parameters = PModelParameters(inputs, FT)

                T_canopy = FT(inputs["tc"] + 273.15)  # Convert from Celsius to Kelvin
                VPD = FT(inputs["vpd"])
                ca = FT(inputs["co2"]) * FT(1e-6) * FT(101325.0)  # Convert ppm to Pa
                P_air = FT(get(inputs, "patm", 101325.0))
                I_abs = FT(inputs["fapar"]) * FT(inputs["ppfd"])
                βm =
                    Bool(inputs["do_soilmstress"]) ?
                    quadratic_soil_moisture_stress(FT(inputs["soilm"])) :
                    FT(1.0)

                # Run the model
                outputs = compute_full_pmodel_outputs(
                    parameters,
                    constants,
                    T_canopy,
                    I_abs,
                    ca,
                    P_air,
                    VPD,
                    βm,
                )

                # Compare each output field
                for key in keys(ref_outputs_typed)
                    # Convert string key to symbol for named tuple access
                    key_symbol = Symbol(key)

                    if haskey(outputs, key_symbol)
                        r_out = ref_outputs_typed[key]
                        # convert gpp to kg/m^2/s   
                        r_out = key == "gpp" ? r_out / FT(1000.0) : r_out

                        j_out = getproperty(outputs, key_symbol)
                        diff = percent_difference(j_out, r_out)

                        # Verbose output
                        if verbose
                            println("Output: $key")
                            println("  Expected: $r_out")
                            println("  Computed: $j_out")
                            println("  Percent Difference: $diff%")
                        end
                        # Test for approximate equality
                        @test isapprox(j_out, r_out, rtol = rtol, atol = atol)
                    else
                        verbose && @warn "Missing key $key in Julia outputs"
                    end
                end
            end
        end
    end
end
