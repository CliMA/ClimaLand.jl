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
            @test params.sigma isa FT
            @test params.alpha isa FT

            # Check expected values from default_parameters.toml
            @test params.k ≈ FT(0.5)
            @test params.z ≈ FT(12.227)
            @test params.sigma ≈ FT(1.1)
            @test params.alpha ≈ FT(0.067)  # ~15 days of memory

            @test eltype(params) == FT
        end

        @testset "OptimalLAIModel construction for FT = $FT" begin
            params = Canopy.OptimalLAIParameters{FT}(toml_dict)
            model = Canopy.OptimalLAIModel{FT}(params)

            @test model.parameters === params
            @test eltype(model) == FT

            # Test auxiliary variables (LAI + A0 tracking variables)
            @test Canopy.auxiliary_vars(model) == (
                :LAI,
                :A0_daily,
                :A0_annual,
                :A0_daily_acc,
                :A0_annual_acc,
                :A0_annual_daily_acc,
                :last_day_of_year,
            )
            @test Canopy.auxiliary_types(model) == (FT, FT, FT, FT, FT, FT, FT)
            @test Canopy.auxiliary_domain_names(model) ==
                  (:surface, :surface, :surface, :surface, :surface, :surface, :surface)
        end

        @testset "compute_L_max function (energy-limited only) for FT = $FT" begin
            # Test with typical conditions
            Ao_annual = FT(100.0)   # mol m⁻² yr⁻¹
            k = FT(0.5)
            z = FT(12.227)

            LAI_max = Canopy.compute_L_max(Ao_annual, k, z)

            @test LAI_max isa FT
            @test LAI_max >= FT(0.0)  # LAI should be non-negative
            @test LAI_max < FT(20.0)  # LAI should be reasonable (< 20)

            # Test that higher A0_annual gives higher LAI_max
            LAI_high_gpp = Canopy.compute_L_max(FT(300.0), k, z)
            LAI_low_gpp = Canopy.compute_L_max(FT(50.0), k, z)
            @test LAI_high_gpp > LAI_low_gpp

            # Test energy limitation formula
            fAPAR_energy = FT(1) - z / (k * Ao_annual)
            fAPAR_max = max(FT(0), min(FT(1), fAPAR_energy))
            LAI_max_manual = -(FT(1) / k) * log(FT(1) - fAPAR_max)
            @test LAI_max ≈ LAI_max_manual
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

            # Test near branch point - at x = -1/e + 1e-8, W(x) ≈ -1 + sqrt(2*1e-8*e)
            # For Float64: W(-1/e + 1e-8) ≈ -0.9997668
            # For Float32: -1/e + 1e-8 rounds to exactly -1/e, so W(-1/e) = -1
            x_near_branch = -FT(1.0) / FT(ℯ) + FT(1e-8)
            w_near_branch = Canopy.lambertw0(x_near_branch)
            @test w_near_branch ≈ -FT(1.0) atol = FT(1e-3)  # Looser tolerance near branch point

            # Test invalid input returns NaN (GPU-friendly behavior)
            @test isnan(Canopy.lambertw0(-FT(1.0)))
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

        @testset "compute_PPFD function for FT = $FT" begin
            # Test PPFD computation from PAR
            par_d = FT(500.0)  # W m⁻² (typical midday)
            λ_γ_PAR = FT(5e-7)  # 500 nm
            lightspeed = FT(3e8)  # m s⁻¹
            planck_h = FT(6.626e-34)  # J s
            N_a = FT(6.022e23)  # mol⁻¹

            PPFD = Canopy.compute_PPFD(par_d, λ_γ_PAR, lightspeed, planck_h, N_a)

            @test PPFD isa FT
            @test PPFD > FT(0.0)
            @test isfinite(PPFD)
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
            A0_daily = FT(0.5)  # mol CO₂ m⁻² day⁻¹ (with actual β)
            L = FT(2.0)   # m² m⁻²

            L_new = Canopy.update_optimal_LAI(
                local_noon_mask,
                A0_daily,
                L,
                FT(0.5),       # k
                FT(100.0),     # A0_annual (with β=1)
                FT(12.227),    # z
                FT(180.0),     # GSL
                FT(0.771),     # sigma
                FT(0.017),     # alpha
            )

            @test L_new isa FT
            @test L_new >= FT(0.0)
            @test L_new < FT(20.0)  # Reasonable LAI range

            # Test with no update (mask = 0)
            L_no_update = Canopy.update_optimal_LAI(
                FT(0.0),       # local_noon_mask = 0
                A0_daily,
                L,
                FT(0.5),       # k
                FT(100.0),     # A0_annual
                FT(12.227),    # z
                FT(180.0),     # GSL
                FT(0.771),     # sigma
                FT(0.017),     # alpha
            )
            @test L_no_update == L
        end

        @testset "update_A0_and_LAI_at_noon function for FT = $FT" begin
            # Test accumulation at local noon
            A0_daily = FT(0.5)
            A0_annual = FT(100.0)
            A0_daily_acc = FT(0.6)  # Accumulated during the day (with actual β)
            A0_annual_acc = FT(50.0)  # Accumulated during the year (with β=1)
            A0_annual_daily_acc = FT(0.7)  # Daily accumulator for β=1 GPP
            last_doy = FT(100.0)
            current_doy = FT(101.0)
            L = FT(2.0)

            # Optimal LAI parameters
            k = FT(0.5)
            z = FT(12.227)
            GSL = FT(180.0)
            sigma = FT(0.771)
            alpha = FT(0.067)

            # At local noon, should finalize daily A0 and accumulate to annual
            new_daily, new_annual, new_daily_acc, new_annual_acc, new_annual_daily_acc, new_last_doy, new_L =
                Canopy.update_A0_and_LAI_at_noon(
                    FT(1.0),  # local noon mask
                    A0_daily,
                    A0_annual,
                    A0_daily_acc,
                    A0_annual_acc,
                    A0_annual_daily_acc,
                    last_doy,
                    current_doy,
                    L,
                    k,
                    z,
                    GSL,
                    sigma,
                    alpha,
                )

            @test new_daily ≈ A0_daily_acc  # Daily A0 (with β) is finalized from accumulator
            @test new_daily_acc ≈ FT(0.0)   # Daily accumulator reset
            @test new_annual_daily_acc ≈ FT(0.0)  # Annual daily accumulator reset
            @test new_annual_acc ≈ A0_annual_acc + A0_annual_daily_acc  # Annual accumulator updated with β=1 daily
            @test new_last_doy ≈ current_doy
            @test new_L isa FT  # LAI is updated

            # Test year change detection (day of year decreased)
            _, new_annual_year_change, _, _, _, _, _ =
                Canopy.update_A0_and_LAI_at_noon(
                    FT(1.0),
                    A0_daily,
                    A0_annual,
                    A0_daily_acc,
                    A0_annual_acc,
                    A0_annual_daily_acc,
                    FT(365.0),  # last day of year
                    FT(1.0),    # first day of new year
                    L,
                    k,
                    z,
                    GSL,
                    sigma,
                    alpha,
                )

            @test new_annual_year_change ≈ A0_annual_acc  # Annual finalized before new year

            # Test no update when not at noon
            same_daily, same_annual, same_daily_acc, same_annual_acc, same_annual_daily_acc, same_last_doy, same_L =
                Canopy.update_A0_and_LAI_at_noon(
                    FT(0.0),  # not at noon
                    A0_daily,
                    A0_annual,
                    A0_daily_acc,
                    A0_annual_acc,
                    A0_annual_daily_acc,
                    last_doy,
                    current_doy,
                    L,
                    k,
                    z,
                    GSL,
                    sigma,
                    alpha,
                )

            @test same_daily == A0_daily
            @test same_annual == A0_annual
            @test same_daily_acc == A0_daily_acc
            @test same_annual_acc == A0_annual_acc
            @test same_annual_daily_acc == A0_annual_daily_acc
            @test same_L == L  # LAI unchanged when not at noon
        end

        @testset "Full LAI computation trace with Ozark-like parameters for FT = $FT" begin
            # Parameters matching Ozark site conditions
            A0_annual = FT(258.0)   # mol CO₂ m⁻² yr⁻¹ (with β=1)
            k = FT(0.5)
            z = FT(12.227)          # mol m⁻² yr⁻¹
            GSL = FT(180.0)         # days
            sigma = FT(0.771)
            alpha = FT(0.067)       # ~15 days memory

            # Step 1: Compute LAI_max (energy-limited only)
            fAPAR_energy = FT(1) - z / (k * A0_annual)
            fAPAR_max = max(FT(0), min(FT(1), fAPAR_energy))
            LAI_max_manual = -(FT(1) / k) * log(FT(1) - fAPAR_max)

            LAI_max_func = Canopy.compute_L_max(A0_annual, k, z)

            @test LAI_max_manual ≈ LAI_max_func
            @test LAI_max_func > FT(4.0)  # Ozark should support LAI > 4
            @test LAI_max_func < FT(6.0)  # But not unreasonably high

            # Step 2: Compute m
            fAPAR_max_from_LAI = FT(1) - exp(-k * LAI_max_func)
            m_manual = (sigma * GSL * LAI_max_func) / (A0_annual * fAPAR_max_from_LAI)

            m_func = Canopy.compute_m(GSL, LAI_max_func, A0_annual, sigma, k)

            @test m_manual ≈ m_func
            @test m_func > FT(2.0)  # Should be reasonable positive value
            @test m_func < FT(5.0)

            # Step 3: Compute L_steady for a summer day with A0_daily = 1.5 (with actual β)
            A0_daily_summer = FT(1.5)  # mol CO₂ m⁻² day⁻¹
            mu = m_func * A0_daily_summer
            arg = -k * mu * exp(-k * mu)
            w_val = Canopy.lambertw0(arg)
            L_steady_manual = mu + (FT(1) / k) * w_val
            L_steady_manual = min(L_steady_manual, LAI_max_func)
            L_steady_manual = max(FT(0), L_steady_manual)

            L_steady_func = Canopy.compute_steady_state_LAI(A0_daily_summer, m_func, k, LAI_max_func)

            @test L_steady_manual ≈ L_steady_func
            @test L_steady_func > FT(2.5)  # Summer should give reasonable LAI
            @test L_steady_func < FT(5.0)

            # Step 4: Test update_optimal_LAI integrates correctly
            L_prev = FT(1.0)
            L_new = Canopy.update_optimal_LAI(
                FT(1.0),       # local_noon_mask
                A0_daily_summer,
                L_prev,
                k,
                A0_annual,
                z,
                GSL,
                sigma,
                alpha,
            )

            # Manual calculation of expected new LAI
            L_expected = alpha * L_steady_func + (FT(1) - alpha) * L_prev
            @test L_new ≈ L_expected
            @test L_new > L_prev  # LAI should increase toward summer steady state

            # Step 5: Test that after many iterations, LAI converges to L_steady
            # With alpha = 0.067 (~15 days), converges in ~100 days
            L_current = FT(1.0)
            for _ in 1:100  # 100 days of updates for convergence
                L_current = Canopy.update_optimal_LAI(
                    FT(1.0),
                    A0_daily_summer,
                    L_current,
                    k,
                    A0_annual,
                    z,
                    GSL,
                    sigma,
                    alpha,
                )
            end
            @test L_current ≈ L_steady_func atol = FT(0.1)  # Should converge to steady state

            # Step 6: Test with winter-like A0_daily = 0.2 (low β due to cold/dry)
            A0_daily_winter = FT(0.2)
            L_steady_winter = Canopy.compute_steady_state_LAI(A0_daily_winter, m_func, k, LAI_max_func)
            @test L_steady_winter < FT(1.0)  # Winter should have low steady-state LAI
            @test L_steady_winter >= FT(0.0)

            # Print diagnostic values for debugging
            @info "Ozark-like LAI computation trace" LAI_max=LAI_max_func m=m_func L_steady_summer=L_steady_func L_steady_winter=L_steady_winter fAPAR_energy=fAPAR_energy
        end

        @testset "Water limitation through daily β test for FT = $FT" begin
            # Test that varying β in daily A0 affects LAI dynamics
            # This demonstrates the new water limitation mechanism

            A0_annual = FT(200.0)  # High potential GPP (with β=1)
            k = FT(0.5)
            z = FT(12.227)
            GSL = FT(180.0)
            sigma = FT(0.771)
            alpha = FT(0.067)

            LAI_max = Canopy.compute_L_max(A0_annual, k, z)
            m = Canopy.compute_m(GSL, LAI_max, A0_annual, sigma, k)

            # Wet conditions: high β → high daily A0
            A0_daily_wet = FT(1.0)  # mol CO₂ m⁻² day⁻¹
            L_steady_wet = Canopy.compute_steady_state_LAI(A0_daily_wet, m, k, LAI_max)

            # Dry conditions: low β → low daily A0
            A0_daily_dry = FT(0.3)  # mol CO₂ m⁻² day⁻¹
            L_steady_dry = Canopy.compute_steady_state_LAI(A0_daily_dry, m, k, LAI_max)

            # Wet conditions should support higher steady-state LAI
            @test L_steady_wet > L_steady_dry

            # Simulate a "flushing event" - dry conditions followed by wet
            L_current = FT(0.5)  # Low initial LAI from dry period

            # 30 days of dry conditions
            for _ in 1:30
                L_current = Canopy.update_optimal_LAI(
                    FT(1.0), A0_daily_dry, L_current, k, A0_annual, z, GSL, sigma, alpha)
            end
            L_after_dry = L_current

            # Then 60 days of wet conditions (flushing event)
            for _ in 1:60
                L_current = Canopy.update_optimal_LAI(
                    FT(1.0), A0_daily_wet, L_current, k, A0_annual, z, GSL, sigma, alpha)
            end
            L_after_wet = L_current

            # LAI should increase significantly after wet period
            @test L_after_wet > L_after_dry

            @info "Flushing event test" L_steady_dry=L_steady_dry L_steady_wet=L_steady_wet L_after_dry=L_after_dry L_after_wet=L_after_wet
        end
    end
end
