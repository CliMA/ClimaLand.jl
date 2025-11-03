using Test
using ClimaLand
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand.Canopy

@testset "GPU-compatible lambertw0 function tests" begin
    @testset "CPU tests for FT = $FT" for FT in (Float32, Float64)
        # Define relative tolerance based on precision
        rtol = FT == Float32 ? FT(1e-6) : FT(1e-12)

        # Test values in the valid domain
        test_values = [
            (-1 / exp(1) + FT(1e-7), -FT(1)),  # Near branch point
            (-FT(0.1), -FT(0.11183255915896293)),
            (FT(0.0), FT(0.0)),
            (FT(0.1), FT(0.09127652716086226)),
            (FT(1.0), FT(0.5671432904097838)),
            (FT(10.0), FT(1.7455280027406994)),
        ]

        @testset "lambertw0 accuracy for x = $x" for (x, expected) in test_values
            result = Canopy.lambertw0(FT(x))
            @test result isa FT
            @test isapprox(result, FT(expected), rtol = rtol)
        end

        # Test invalid inputs return NaN
        @testset "Invalid inputs return NaN" begin
            @test isnan(Canopy.lambertw0(FT(-1.0)))  # x < -1/e
            @test isnan(Canopy.lambertw0(FT(NaN)))    # NaN input
            @test isnan(Canopy.lambertw0(FT(Inf)))    # Inf input (should handle gracefully)
        end

        # Test broadcastability on CPU arrays
        @testset "Broadcasting on CPU arrays" begin
            x_vals = FT[-0.3, -0.1, 0.0, 0.1, 1.0, 10.0]
            results = Canopy.lambertw0.(x_vals)
            @test results isa Vector{FT}
            @test length(results) == length(x_vals)
            @test all(isfinite.(results))
        end
    end

    # GPU tests - only run if CUDA is available
    @testset "GPU tests" begin
        device = ClimaComms.device()
        
        if device isa ClimaComms.CUDADevice
            @testset "GPU broadcasting for Float32" begin
                FT = Float32
                ArrayType = ClimaComms.array_type(device)
                
                # Create test data on CPU
                x_cpu = FT[-0.3, -0.1, 0.0, 0.1, 1.0, 10.0]
                expected_cpu = Canopy.lambertw0.(x_cpu)
                
                # Transfer to GPU
                x_gpu = ArrayType(x_cpu)
                
                # Compute on GPU
                results_gpu = Canopy.lambertw0.(x_gpu)
                
                # Transfer back to CPU for comparison
                results_cpu = Array(results_gpu)
                
                # Compare with CPU results
                @test results_cpu isa Vector{FT}
                @test isapprox(results_cpu, expected_cpu, rtol = FT(1e-5))
            end
        else
            @info "Skipping GPU tests: CUDA not available (device: $device)"
        end
    end
end
