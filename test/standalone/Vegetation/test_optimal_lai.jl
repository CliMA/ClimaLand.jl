using Test
using ClimaLand
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand.Canopy
using ClimaLand.Domains: Point
import ClimaLand.Parameters as LP
import ClimaParams
using ClimaCore

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

        @testset "ZhouOptimalLAIModel construction for FT = $FT" begin
            params = Canopy.OptimalLAIParameters{FT}(toml_dict)
            # For unit tests, use scalar values for initial conditions
            optimal_lai_inputs = (;
                GSL = FT(240.0),
                A0_annual = FT(258.0),
                precip_annual = FT(1000.0),
                vpd_gs = FT(1000.0),
                lai_init = FT(2.0),
                f0 = FT(0.65),
            )
            model = Canopy.ZhouOptimalLAIModel{FT}(
                params,
                optimal_lai_inputs;
                SAI = FT(0.0),
                RAI = FT(1.0),
                rooting_depth = FT(1.0),
                height = FT(10.0),
            )

            @test model.parameters === params
            @test model.optimal_lai_inputs === optimal_lai_inputs
            @test eltype(model) == FT
            @test model.SAI == FT(0.0)
            @test model.RAI == FT(1.0)

            # Test auxiliary variables
            aux_vars = Canopy.auxiliary_vars(model)
            @test :area_index in aux_vars
            @test :A0_daily in aux_vars
            @test :A0_annual in aux_vars
            @test :A0_daily_acc in aux_vars
            @test :A0_annual_acc in aux_vars
            @test :days_since_reset in aux_vars
            @test :GSL in aux_vars
            @test :precip_annual in aux_vars
            @test :vpd_gs in aux_vars
            @test :f0 in aux_vars
        end

        @testset "compute_L_max function (energy-limited only) for FT = $FT" begin
            # Test with typical conditions
            Ao_annual = FT(100.0)   # mol m^-2 yr^-1
            k = FT(0.5)
            z = FT(12.227)
            # Use high precip and low VPD to ensure energy-limited
            precip_annual = FT(100000.0)  # mol H2O m^-2 yr^-1 (very high)
            f0 = FT(0.65)
            ca_pa = FT(40.0)  # Pa
            chi = FT(0.77)  # typical tropical value
            vpd_gs = FT(1000.0)  # Pa

            LAI_max = Canopy.compute_L_max(
                Ao_annual,
                k,
                z,
                precip_annual,
                f0,
                ca_pa,
                chi,
                vpd_gs,
            )

            @test LAI_max isa FT
            @test LAI_max >= FT(0.0)  # LAI should be non-negative
            @test LAI_max < FT(20.0)  # LAI should be reasonable (< 20)

            # Test that higher A0_annual gives higher LAI_max
            LAI_high_gpp = Canopy.compute_L_max(
                FT(300.0),
                k,
                z,
                precip_annual,
                f0,
                ca_pa,
                chi,
                vpd_gs,
            )
            LAI_low_gpp = Canopy.compute_L_max(
                FT(50.0),
                k,
                z,
                precip_annual,
                f0,
                ca_pa,
                chi,
                vpd_gs,
            )
            @test LAI_high_gpp > LAI_low_gpp

            # Test energy limitation formula (with high precip, should be energy-limited)
            fAPAR_energy = FT(1) - z / (k * Ao_annual)
            fAPAR_max = max(FT(0), min(FT(1), fAPAR_energy))
            LAI_max_manual = -(FT(1) / k) * log(FT(1) - fAPAR_max)
            @test LAI_max ≈ LAI_max_manual
        end

        @testset "compute_m function for FT = $FT" begin
            GSL = FT(180.0)         # days
            LAI_max = FT(3.0)       # m^2 m^-2
            Ao_annual = FT(100.0)   # mol m^-2 yr^-1
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

            # Test near branch point - at x = -1/e + 1e-8, W(x) ~ -1 + sqrt(2*1e-8*e)
            # For Float64: W(-1/e + 1e-8) ~ -0.9997668
            # For Float32: -1/e + 1e-8 rounds to exactly -1/e, so W(-1/e) = -1
            x_near_branch = -FT(1.0) / FT(ℯ) + FT(1e-8)
            w_near_branch = Canopy.lambertw0(x_near_branch)
            @test w_near_branch ≈ -FT(1.0) atol = FT(1e-3)  # Looser tolerance near branch point

            # Test invalid input returns NaN (GPU-friendly behavior)
            @test isnan(Canopy.lambertw0(-FT(1.0)))
        end

        @testset "compute_steady_state_LAI function for FT = $FT" begin
            Ao_daily = FT(0.4)      # mol m^-2 day^-1
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

        @testset "compute_LAI function for FT = $FT" begin
            LAI_prev = FT(2.0)
            L_steady = FT(3.0)
            alpha = FT(0.067)

            # Test with local noon mask = 1 (update)
            LAI_new = Canopy.compute_LAI(LAI_prev, L_steady, alpha, FT(1.0))
            expected = alpha * L_steady + (1 - alpha) * LAI_prev
            @test LAI_new ≈ expected
            @test LAI_new > LAI_prev  # Should move toward higher steady state
            @test LAI_new < L_steady  # But not reach it in one step

            # Test with local noon mask = 0 (no update)
            LAI_no_update =
                Canopy.compute_LAI(LAI_prev, L_steady, alpha, FT(0.0))
            @test LAI_no_update == LAI_prev

            # Test that LAI is non-negative
            LAI_negative_test =
                Canopy.compute_LAI(-FT(1.0), FT(0.0), alpha, FT(1.0))
            @test LAI_negative_test >= FT(0.0)
        end

        @testset "compute_PPFD function for FT = $FT" begin
            # Test PPFD computation from PAR
            par_d = FT(500.0)  # W m^-2 (typical midday)
            λ_γ_PAR = FT(5e-7)  # 500 nm
            lightspeed = FT(3e8)  # m s^-1
            planck_h = FT(6.626e-34)  # J s
            N_a = FT(6.022e23)  # mol^-1

            PPFD =
                Canopy.compute_PPFD(par_d, λ_γ_PAR, lightspeed, planck_h, N_a)

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
            mask_morning = Canopy.get_local_noon_mask(21600.0, dt, local_noon)  # 6 AM
            @test mask_morning == FT(0.0)

            mask_evening = Canopy.get_local_noon_mask(64800.0, dt, local_noon)  # 6 PM
            @test mask_evening == FT(0.0)
        end

        @testset "update_optimal_LAI function for FT = $FT" begin
            # Test the full LAI update function
            A0_daily = FT(0.5)       # mol m^-2 day^-1
            L = FT(2.0)              # current LAI
            k = FT(0.5)
            A0_annual = FT(100.0)    # mol m^-2 yr^-1
            z = FT(12.227)
            GSL = FT(180.0)          # days
            sigma = FT(0.771)
            alpha = FT(0.067)
            precip_annual = FT(100000.0)  # mol H2O m^-2 yr^-1 (high, energy-limited)
            f0 = FT(0.65)
            ca_pa = FT(40.0)         # Pa
            chi = FT(0.77)
            vpd_gs = FT(1000.0)      # Pa

            # Test update at noon
            L_new = Canopy.update_optimal_LAI(
                FT(1.0),  # local noon mask
                A0_daily,
                L,
                k,
                A0_annual,
                z,
                GSL,
                sigma,
                alpha,
                precip_annual,
                f0,
                ca_pa,
                chi,
                vpd_gs,
            )

            @test L_new isa FT
            @test L_new >= FT(0.0)
            @test isfinite(L_new)

            # Test no update when not at noon
            L_no_update = Canopy.update_optimal_LAI(
                FT(0.0),  # not local noon
                A0_daily,
                L,
                k,
                A0_annual,
                z,
                GSL,
                sigma,
                alpha,
                precip_annual,
                f0,
                ca_pa,
                chi,
                vpd_gs,
            )
            @test L_no_update == L
        end

        @testset "optimal_lai_initial_conditions for single-point domains for FT = $FT" begin
            # Test that optimal_lai_initial_conditions returns reasonable values
            # for single-point domains at various locations (Fluxnet sites)

            test_sites = [
                # (name, longitude, latitude)
                ("US-MOz (Ozark)", FT(-92.2000), FT(38.7441)),   # Missouri, USA - Deciduous forest
                ("US-Ha1 (Harvard)", FT(-72.1715), FT(42.5378)), # Massachusetts, USA - Mixed forest
                ("Amazon", FT(-60.0), FT(-3.0)),                 # Amazon rainforest
                ("Sahel", FT(0.0), FT(15.0)),                    # Semi-arid Africa
            ]

            for (site_name, long, lat) in test_sites
                # Create a point domain at this location
                domain = Point(; z_sfc = FT(0.0), longlat = (long, lat))
                surface_space = domain.space.surface

                # Load initial conditions from global data file
                optimal_lai_inputs =
                    Canopy.optimal_lai_initial_conditions(surface_space)

                # Extract scalar values from Fields
                GSL_val = Array(parent(optimal_lai_inputs.GSL))[1]
                A0_annual_val = Array(parent(optimal_lai_inputs.A0_annual))[1]
                precip_annual_val =
                    Array(parent(optimal_lai_inputs.precip_annual))[1]
                vpd_gs_val = Array(parent(optimal_lai_inputs.vpd_gs))[1]
                lai_init_val = Array(parent(optimal_lai_inputs.lai_init))[1]
                f0_val = Array(parent(optimal_lai_inputs.f0))[1]

                @testset "$site_name" begin
                    # GSL should be positive and reasonable (0-365 days)
                    @test GSL_val >= FT(0) "GSL should be non-negative at $site_name"
                    @test GSL_val <= FT(365) "GSL should be <= 365 days at $site_name"
                    @test GSL_val > FT(0) "GSL should be positive (non-zero) at $site_name, got $GSL_val"

                    # A0_annual should be positive (mol CO2 m^-2 yr^-1)
                    # Typical values range from ~50 (arid) to ~500 (tropical rainforest)
                    @test A0_annual_val >= FT(0) "A0_annual should be non-negative at $site_name"
                    @test A0_annual_val > FT(0) "A0_annual should be positive at $site_name, got $A0_annual_val"
                    @test A0_annual_val < FT(1000) "A0_annual should be < 1000 at $site_name"

                    # precip_annual should be positive (mol H2O m^-2 yr^-1)
                    # Ranges from ~5000 (desert, ~100 mm) to ~170000+ (tropical, ~3000 mm)
                    @test precip_annual_val >= FT(0) "precip_annual should be non-negative at $site_name"
                    @test precip_annual_val > FT(0) "precip_annual should be positive at $site_name, got $precip_annual_val"

                    # vpd_gs should be positive (Pa)
                    # Typical growing season VPD: 500-2500 Pa
                    @test vpd_gs_val >= FT(0) "vpd_gs should be non-negative at $site_name"
                    @test vpd_gs_val > FT(0) "vpd_gs should be positive at $site_name, got $vpd_gs_val"

                    # lai_init should be non-negative (m^2 m^-2)
                    # Ranges from 0 (bare) to ~8 (dense forest)
                    @test lai_init_val >= FT(0) "lai_init should be non-negative at $site_name"
                    @test lai_init_val < FT(15) "lai_init should be < 15 at $site_name"

                    # f0 should be in range [0, 1]
                    @test f0_val >= FT(0) "f0 should be >= 0 at $site_name"
                    @test f0_val <= FT(1) "f0 should be <= 1 at $site_name"
                    @test f0_val > FT(0) "f0 should be positive at $site_name, got $f0_val"
                end
            end
        end
    end
end
