using Test
using ClimaLand
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand.Canopy
using ClimaLand.Domains: Point
import ClimaLand.Parameters as LP
import ClimaParams

@testset "Optimal LAI Model Tests" begin
    for FT in (Float32, Float64)
        toml_dict = LP.create_toml_dict(FT)

        @testset "OptimalLAIParameters construction for FT = $FT" begin
            # Test parameter construction from TOML
            params = Canopy.OptimalLAIParameters{FT}(toml_dict)

            @test params.k isa FT
            @test params.z isa FT
            @test params.chi isa FT
            @test params.f0 isa FT
            @test params.sigma isa FT
            @test params.alpha isa FT

            # Check expected values from default_parameters.toml
            @test params.k ≈ FT(0.5)
            @test params.z ≈ FT(12.227)
            @test params.chi ≈ FT(0.7)
            @test params.f0 ≈ FT(0.62)
            @test params.sigma ≈ FT(0.771)
            @test params.alpha ≈ FT(0.067)

            @test eltype(params) == FT
        end

        @testset "OptimalLAIModel construction for FT = $FT" begin
            params = Canopy.OptimalLAIParameters{FT}(toml_dict)
            model = Canopy.OptimalLAIModel{FT}(params)

            @test model.parameters === params
            @test eltype(model) == FT

            # Test auxiliary variables
            @test Canopy.auxiliary_vars(model) == (:LAI,)
            @test Canopy.auxiliary_types(model) == (FT,)
            @test Canopy.auxiliary_domain_names(model) == (:surface,)
        end

        @testset "compute_L_max function for FT = $FT" begin
            # Test with typical temperate forest conditions
            Ao_annual = FT(100.0)   # mol m⁻² yr⁻¹
            P_annual = FT(60000.0)  # mol m⁻² yr⁻¹ (~1080 mm)
            D_growing = FT(1000.0)  # Pa
            k = FT(0.5)
            z = FT(12.227)
            ca = FT(40.0)           # Pa (~400 ppm)
            chi = FT(0.7)
            f0 = FT(0.62)

            LAI_max = Canopy.compute_L_max(
                Ao_annual,
                P_annual,
                D_growing,
                k,
                z,
                ca,
                chi,
                f0,
            )

            @test LAI_max isa FT
            @test LAI_max >= FT(0.0)  # LAI should be non-negative
            @test LAI_max < FT(20.0)  # LAI should be reasonable (< 20)

            # Test water-limited vs energy-limited scenarios
            # Very wet conditions - should be energy limited
            P_wet = FT(200000.0)  # Very high precipitation
            LAI_wet = Canopy.compute_L_max(
                Ao_annual,
                P_wet,
                D_growing,
                k,
                z,
                ca,
                chi,
                f0,
            )

            # Dry conditions - should be water limited
            P_dry = FT(20000.0)  # Low precipitation
            LAI_dry = Canopy.compute_L_max(
                Ao_annual,
                P_dry,
                D_growing,
                k,
                z,
                ca,
                chi,
                f0,
            )

            @test LAI_wet > LAI_dry  # Wetter conditions should support higher LAI
        end

        @testset "compute_m function for FT = $FT" begin
            GSL = FT(180.0)         # days
            LAI_max = FT(3.0)       # m² m⁻²
            Ao_annual = FT(100.0)   # mol m⁻² yr⁻¹
            sigma = FT(0.771)
            k = FT(0.5)

            m = Canopy.compute_m(GSL, LAI_max, Ao_annual, sigma, k)

            @test m isa FT
            @test m > FT(0.0)  # m should be positive

            # Test that m scales with GSL
            m_short = Canopy.compute_m(FT(90.0), LAI_max, Ao_annual, sigma, k)
            m_long = Canopy.compute_m(FT(270.0), LAI_max, Ao_annual, sigma, k)
            @test m_long > m_short  # Longer GSL should give larger m
        end

        @testset "lambertw0 function for FT = $FT" begin
            # Test known values of Lambert W function
            @test Canopy.lambertw0(FT(0.0)) ≈ FT(0.0) atol = FT(1e-6)
            @test Canopy.lambertw0(FT(1.0)) ≈ FT(0.5671432904097838) atol =
                FT(1e-6)
            @test Canopy.lambertw0(FT(ℯ)) ≈ FT(1.0) atol = FT(1e-6)

            # Test near branch point
            @test Canopy.lambertw0(-FT(1.0) / FT(ℯ)) ≈ -FT(1.0) atol = FT(1e-6)

            # Test domain error for invalid input
            @test_throws DomainError Canopy.lambertw0(-FT(1.0))
        end

        @testset "compute_steady_state_LAI function for FT = $FT" begin
            Ao_daily = FT(0.4)      # mol m⁻² day⁻¹
            m = FT(7.0)
            k = FT(0.5)
            LAI_max = FT(3.0)

            L_steady = Canopy.compute_steady_state_LAI(Ao_daily, m, k, LAI_max)

            @test L_steady isa FT
            @test L_steady >= FT(0.0)
            @test L_steady <= LAI_max  # Should not exceed LAI_max

            # Test with zero GPP
            L_zero = Canopy.compute_steady_state_LAI(FT(0.0), m, k, LAI_max)
            @test L_zero ≈ FT(0.0)

            # Test that higher GPP gives higher steady-state LAI
            L_low = Canopy.compute_steady_state_LAI(FT(0.1), m, k, LAI_max)
            L_high = Canopy.compute_steady_state_LAI(FT(0.8), m, k, LAI_max)
            @test L_high > L_low
        end

        @testset "update_LAI! function for FT = $FT" begin
            LAI_prev = FT(2.0)
            L_steady = FT(3.0)
            alpha = FT(0.067)

            # Test with local noon mask = 1 (update)
            LAI_new = Canopy.update_LAI!(LAI_prev, L_steady, alpha, FT(1.0))
            expected = alpha * L_steady + (1 - alpha) * LAI_prev
            @test LAI_new ≈ expected
            @test LAI_new > LAI_prev  # Should move toward higher steady state
            @test LAI_new < L_steady  # But not reach it in one step

            # Test with local noon mask = 0 (no update)
            LAI_no_update =
                Canopy.update_LAI!(LAI_prev, L_steady, alpha, FT(0.0))
            @test LAI_no_update == LAI_prev

            # Test that LAI is non-negative
            LAI_negative_test =
                Canopy.update_LAI!(-FT(1.0), FT(0.0), alpha, FT(1.0))
            @test LAI_negative_test >= FT(0.0)
        end

        @testset "compute_Ao_daily function for FT = $FT" begin
            A = FT(10.0)  # Some instantaneous GPP
            k = FT(0.5)
            L = FT(2.0)

            Ao_daily = Canopy.compute_Ao_daily(A, k, L)

            @test Ao_daily isa FT
            @test Ao_daily > FT(0.0)

            # Test with zero LAI - should handle gracefully
            Ao_zero_lai = Canopy.compute_Ao_daily(A, k, FT(0.0))
            @test isfinite(Ao_zero_lai)
        end

        @testset "get_local_noon_mask function for FT = $FT" begin
            dt = FT(3600.0)  # 1 hour timestep
            local_noon = FT(43200.0)  # 12:00 noon in seconds

            # Test at local noon
            mask_noon = Canopy.get_local_noon_mask(43200.0, dt, local_noon)
            @test mask_noon == FT(1.0)

            # Test within window
            mask_before =
                Canopy.get_local_noon_mask(43200.0 - dt / 4, dt, local_noon)
            @test mask_before == FT(1.0)

            mask_after =
                Canopy.get_local_noon_mask(43200.0 + dt / 4, dt, local_noon)
            @test mask_after == FT(1.0)

            # Test outside window
            mask_outside =
                Canopy.get_local_noon_mask(43200.0 + dt, dt, local_noon)
            @test mask_outside == FT(0.0)
        end

        @testset "update_optimal_LAI integration test for FT = $FT" begin
            local_noon_mask = FT(1.0)
            A = FT(10.0)  # mol m⁻² s⁻¹
            L = FT(2.0)   # m² m⁻²

            L_new = Canopy.update_optimal_LAI(
                local_noon_mask,
                A,
                L;
                k = FT(0.5),
                Ao_annual = FT(100.0),
                P_annual = FT(60000.0),
                D_growing = FT(1000.0),
                z = FT(12.227),
                ca = FT(40.0),
                chi = FT(0.7),
                f0 = FT(0.62),
                GSL = FT(180.0),
                sigma = FT(0.771),
                alpha = FT(0.067),
            )

            @test L_new isa FT
            @test L_new >= FT(0.0)
            @test L_new < FT(20.0)  # Reasonable LAI range

            # Test with no update (mask = 0)
            L_no_update = Canopy.update_optimal_LAI(FT(0.0), A, L; k = FT(0.5))
            @test L_no_update == L
        end
    end
end
