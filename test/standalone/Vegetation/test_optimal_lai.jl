"""
# Optimal LAI Model Tests

This module tests the optimal LAI model implementation from Zhou et al. (2025).
Tests verify the core LAI computation functions and their mathematical properties.

Reference:
Zhou, B., Cai, W., Zhu, Z., Wang, H., Harrison, S. P., & Prentice, C. (2025).
A General Model for the Seasonal to Decadal Dynamics of Leaf Area.
Global Change Biology, 31(1), e70125. https://doi.org/10.1111/gcb.70125
"""

using Test
using ClimaLand
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand.Canopy
import ClimaLand.Parameters as LP

@testset "Optimal LAI Model Functions" begin
    for FT in (Float32, Float64)
        @testset "Test compute_Ao for FT = $FT" begin
            # Test that compute_Ao correctly inverts Beer-Lambert law
            A = FT(10.0)  # Canopy-integrated assimilation
            k = FT(0.5)   # Extinction coefficient
            L = FT(3.0)   # LAI
            
            Ao = Canopy.compute_Ao(A, k, L)
            
            # Verify: A = Ao * (1 - exp(-k*L))
            A_reconstructed = Ao * (1 - exp(-k * L))
            @test isapprox(A, A_reconstructed, rtol = FT(1e-5))
            
            # Test that Ao > A (top of canopy should have higher assimilation)
            @test Ao > A
            
            # Test edge case: very small L should give Ao ≈ A
            L_small = FT(0.01)
            Ao_small = Canopy.compute_Ao(A, k, L_small)
            @test isapprox(Ao_small, A, rtol = FT(0.1))
        end
        
        @testset "Test compute_L_max for FT = $FT" begin
            # Test maximum LAI constraint
            k = FT(0.5)
            z = FT(12.227)  # Minimum assimilation threshold
            Ao = FT(50.0)   # Top-of-canopy assimilation
            
            L_max = Canopy.compute_L_max(k, z, Ao)
            
            # Verify: z = Ao * exp(-k * L_max)
            z_reconstructed = Ao * exp(-k * L_max)
            @test isapprox(z, z_reconstructed, rtol = FT(1e-5))
            
            # L_max should be positive
            @test L_max > 0
            
            # Higher Ao should allow larger L_max
            Ao_high = FT(100.0)
            L_max_high = Canopy.compute_L_max(k, z, Ao_high)
            @test L_max_high > L_max
        end
        
        @testset "Test optimality condition functions for FT = $FT" begin
            # Test g and dgdL functions
            μ = FT(15.0)  # Scaled cost parameter
            k = FT(0.5)
            L = FT(2.0)
            
            g_val = Canopy.g(μ, k, L)
            dgdL_val = Canopy.dgdL(μ, k, L)
            
            # Test derivative using finite differences
            dL = FT(1e-6)
            g_plus = Canopy.g(μ, k, L + dL)
            g_minus = Canopy.g(μ, k, L - dL)
            dgdL_numerical = (g_plus - g_minus) / (2 * dL)
            
            @test isapprox(dgdL_val, dgdL_numerical, rtol = FT(1e-4))
            
            # Test that g is continuous
            @test isfinite(g_val)
            @test isfinite(dgdL_val)
        end
        
        @testset "Test compute_L_opt for FT = $FT" begin
            # Test Newton-Raphson solver for optimal LAI
            μ = FT(15.0)
            k = FT(0.5)
            L_init = FT(1.0)
            
            L_opt = Canopy.compute_L_opt(μ, k, L_init)
            
            # Verify that L_opt satisfies the optimality condition g ≈ 0
            g_at_opt = Canopy.g(μ, k, L_opt)
            @test abs(g_at_opt) < FT(1e-3)
            
            # L_opt should be positive
            @test L_opt > 0
            
            # Test that different initial guesses converge to same solution
            L_init2 = FT(5.0)
            L_opt2 = Canopy.compute_L_opt(μ, k, L_init2)
            @test isapprox(L_opt, L_opt2, rtol = FT(1e-2))
        end
        
        @testset "Test compute_L (EMA update) for FT = $FT" begin
            # Test EMA update at local noon
            L_current = FT(2.0)
            Ao = FT(50.0)
            m = FT(0.3)
            k = FT(0.5)
            α = FT(0.933)
            z = FT(12.227)
            
            # Test at local noon (mask = 1)
            local_noon_mask = FT(1.0)
            L_new = Canopy.compute_L(L_current, Ao, m, k, α, z, local_noon_mask)
            
            # L should be updated
            @test L_new != L_current
            @test L_new > 0
            
            # Test that L stays within reasonable bounds
            @test L_new < FT(20.0)  # Unrealistically high LAI
            
            # Test outside local noon (mask = 0)
            local_noon_mask_off = FT(0.0)
            L_unchanged = Canopy.compute_L(L_current, Ao, m, k, α, z, local_noon_mask_off)
            
            # L should not be updated when mask is 0
            @test L_unchanged == L_current
            
            # Test EMA property: new value is weighted average
            # When α is close to 1, new value should be close to current
            α_high = FT(0.99)
            L_ema = Canopy.compute_L(L_current, Ao, m, k, α_high, z, FT(1.0))
            @test abs(L_ema - L_current) < abs(L_new - L_current)
        end
        
        @testset "Test OptimalLAIParameters for FT = $FT" begin
            # Test parameter structure creation
            toml_dict = LP.create_toml_dict(FT)
            
            params = Canopy.OptimalLAIParameters(toml_dict)
            
            @test params.m ≈ FT(0.3)
            @test params.z ≈ FT(12.227)
            @test params.α ≈ FT(0.933)
            
            # Test with custom parameters
            custom_params = Canopy.OptimalLAIParameters{FT}(
                m = FT(0.25),
                z = FT(10.0),
                α = FT(0.9),
            )
            
            @test custom_params.m == FT(0.25)
            @test custom_params.z == FT(10.0)
            @test custom_params.α == FT(0.9)
        end
        
        @testset "Test physical constraints for FT = $FT" begin
            # Test that functions handle edge cases appropriately
            k = FT(0.5)
            
            # Test with very high assimilation
            Ao_high = FT(1000.0)
            z = FT(12.227)
            L_max_high = Canopy.compute_L_max(k, z, Ao_high)
            @test isfinite(L_max_high)
            @test L_max_high > 0
            
            # Test with very small assimilation (near threshold)
            A_small = FT(0.1)
            L_test = FT(1.0)
            Ao_small = Canopy.compute_Ao(A_small, k, L_test)
            @test isfinite(Ao_small)
            @test Ao_small > 0
            
            # Test compute_L with extreme α values
            L_current = FT(2.0)
            Ao = FT(50.0)
            m = FT(0.3)
            z = FT(12.227)
            
            # α = 0 means no memory (all new)
            α_zero = FT(0.0)
            L_no_memory = Canopy.compute_L(L_current, Ao, m, k, α_zero, z, FT(1.0))
            @test isfinite(L_no_memory)
            
            # α = 1 means all memory (no update)
            α_one = FT(1.0)
            L_all_memory = Canopy.compute_L(L_current, Ao, m, k, α_one, z, FT(1.0))
            @test L_all_memory == L_current
        end
    end
end
