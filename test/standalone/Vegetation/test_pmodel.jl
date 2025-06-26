using Test
using ClimaLand
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand.Canopy
using DelimitedFiles
using ClimaLand.Domains: Point
import ClimaLand.Parameters as LP
import ClimaParams

# Directory containing CSV test cases
datadir = "/Users/yuchenli/Documents/CliMA Land (local)/testcases"
inputs_file = joinpath(datadir, "inputs.csv")
outputs_file = joinpath(datadir, "outputs.csv")

# Allowed relative error (can tighten for Float64)
atol = 1e-3
rtol = 1e-2

function percent_difference(a, b)
    return abs(a - b) / abs(b) * 100
end

function create_pmodel_drivers(inputs::Dict{String, Any}, FT)
    T_canopy = FT(inputs["tc"] + 273.15)  # Convert from Celsius to Kelvin
    VPD = FT(inputs["vpd"])
    ca = FT(inputs["co2"]) * FT(1e-6) * FT(101325.0)  # Convert ppm to Pa
    P_air = FT(get(inputs, "patm", 101325.0))
    
    # Calculate I_abs directly from fapar and ppfd
    I_abs = FT(inputs["fapar"]) *  FT(inputs["ppfd"])
    
    return PModelDrivers(
        T_canopy = T_canopy,
        I_abs = I_abs,
        ca = ca,
        P_air = P_air,
        VPD = VPD
    )
end

function create_pmodel_parameters(inputs::Dict{String, Any}, FT)
    β = FT(inputs["beta"])
    cstar = FT(0.41) 
    sc = FT(0.0)  # placeholder since this isn't used
    pc = FT(0.0)  # placeholder since this isn't used

    if Bool(inputs["do_ftemp_kphio"])
        ϕ0 = FT(NaN) 
        ϕc = FT(1.0) 
    else
        ϕ0 = FT(inputs["kphio"])
        ϕc = FT(NaN)
    end

    return PModelParameters(
        sc = sc,
        pc = pc,
        cstar = cstar,
        β = β,
        ϕc = ϕc,
        ϕ0 = ϕ0
    )
end


@testset "P-model regression tests against R output" begin
    # Read inputs and outputs CSV files
    inputs_data, inputs_header = readdlm(inputs_file, ',', header=true)
    outputs_data, outputs_header = readdlm(outputs_file, ',', header=true)
    
    inputs_header = vec(inputs_header)
    outputs_header = vec(outputs_header)
    
    # Get testcase column indices
    inputs_testcase_idx = findfirst(h -> string(h) == "testcase", inputs_header)
    outputs_testcase_idx = findfirst(h -> string(h) == "testcase", outputs_header)
    
    # Process each test case
    for row_idx in 1:size(inputs_data, 1)
        testcase_name = inputs_data[row_idx, inputs_testcase_idx]
        
        # Find corresponding output row
        output_row_idx = findfirst(r -> outputs_data[r, outputs_testcase_idx] == testcase_name, 1:size(outputs_data, 1))
        
        if output_row_idx === nothing
            @warn "No output data found for testcase: $testcase_name"
            continue
        end
        
        # Create inputs dictionary
        inputs = Dict{String, Any}()
        for (col_idx, col_name) in enumerate(inputs_header)
            col_name_str = string(col_name)
            value = inputs_data[row_idx, col_idx]
            
            # Skip testcase column and missing values
            if col_name_str == "testcase" || ismissing(value)
                continue
            end
            
            inputs[col_name_str] = value
        end
        
        # Create expected outputs dictionary
        ref_outputs = Dict{String, Any}()
        for (col_idx, col_name) in enumerate(outputs_header)
            col_name_str = string(col_name)
            value = outputs_data[output_row_idx, col_idx]
            
            # Skip testcase column and missing values
            if col_name_str == "testcase" || ismissing(value)
                continue
            end
            
            ref_outputs[col_name_str] = value
        end

        for FT in (Float32, )
            # Convert ref_outputs to the correct FT
            ref_outputs_typed = Dict{String, FT}(k => FT(v) for (k, v) in ref_outputs)

            @testset "Test $testcase_name, FT = $FT" begin
                # Create constants, drivers, and parameters for the current FT
                constants = create_pmodel_constants(FT)
                drivers = create_pmodel_drivers(inputs, FT)
                parameters = create_pmodel_parameters(inputs, FT)
                
                # Run the model
                outputs = compute_pmodel_outputs(parameters, drivers, constants)
                
                # Convert outputs dict keys to strings for comparison
                outputs = Dict(string(k) => v for (k, v) in outputs)

                # Compare each output field
                for key in keys(ref_outputs_typed)
                    if haskey(outputs, key)
                        r_out = ref_outputs_typed[key]
                        # convert gpp to kg/m^2/s   
                        if key == "gpp"
                            r_out = r_out / 1000.0 
                        end

                        j_out = outputs[key]
                        diff = percent_difference(j_out, r_out)

                        # Verbose output
                        println("Output: $key")
                        println("  Expected: $r_out")
                        println("  Computed: $j_out")
                        println("  Percent Difference: $diff%")

                        # Test for approximate equality
                        @test isapprox(j_out, r_out, rtol=rtol, atol=atol)
                    else
                        @warn "Missing key $key in Julia outputs"
                    end
                end
            end
        end
    end
end
