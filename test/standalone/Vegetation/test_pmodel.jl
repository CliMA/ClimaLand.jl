using Test
using ClimaLand
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand.Canopy
using DelimitedFiles
using ClimaLand.Domains: Point
import ClimaLand.Parameters as LP
import ClimaParams
using JSON3
using Glob

# Directory containing R-generated JSON test cases
datadir = "/Users/yuchenli/Documents/CliMA Land (local)/testcases"
json_files = Glob.glob("*.json", datadir)

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
    f_abs = FT(inputs["fapar"])
    par_d = FT(inputs["ppfd"])
    P_air = FT(get(inputs, "patm", 101325.0))

    return PModelDrivers(
        T_canopy = T_canopy,
        f_abs = f_abs,
        ca = ca,
        P_air = P_air,
        par_d = par_d,
        VPD = VPD
    )
end

function create_pmodel_parameters(inputs::Dict{String, Any}, FT)
    β = FT(inputs["beta"])
    cstar = FT(0.41) 
    sc = FT(0.0)  # placeholder since this isn't used
    pc = FT(0.0)  # placeholder since this isn't used

    if inputs["do_ftemp_kphio"]
        ϕ0 = NaN 
        ϕc = FT(1.0) 
    else
        ϕ0 = FT(inputs["kphio"])
        ϕc = NaN
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
    for file in json_files
        testcase = JSON3.read(read(file, String))
        inputs = Dict{String, Any}(string(k) => v for (k, v) in testcase["inputs"])

        for FT in (Float32, )
            ref_outputs = Dict{String, FT}(string(k) => FT(v) for (k, v) in testcase["outputs"])  # Ensure ref_outputs is String -> FT

            @testset "Test $(basename(file)), FT = $FT" begin
                # Create constants, drivers, and parameters for the current FT
                constants = create_pmodel_constants(FT)
                drivers = create_pmodel_drivers(inputs, FT)
                parameters = create_pmodel_parameters(inputs, FT)

                # Run the model
                outputs = compute_pmodel_outputs(parameters, drivers, constants)
                outputs = Dict(string(k) => v for (k, v) in outputs)  # Ensure outputs is String -> FT

                # Compare each output field
                for key in keys(ref_outputs)
                    if haskey(outputs, key)
                        r_out = ref_outputs[key]
                        j_out = outputs[key]
                        diff = percent_difference(j_out, r_out)

                        # Verbose output
                        println("Key: $key")
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
