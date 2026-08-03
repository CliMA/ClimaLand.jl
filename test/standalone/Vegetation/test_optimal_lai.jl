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
            @test params.tau_long_term isa FT

            # Check expected values from default_parameters.toml (calibrated
            # against MODIS LAI in #1794; headline config promoted in #1815)
            @test params.k ≈ FT(0.5)
            @test params.z ≈ FT(15.0)
            @test params.sigma ≈ FT(0.939)
            @test params.alpha ≈ FT(0.0701)  # ~14 days of memory
            @test params.tau_long_term ≈ FT(6.3072e7)  # 2 years

            # C3/C4 competition is on by default.
            @test params.online_c3c4 ≈ FT(1)

            @test eltype(params) == FT
        end

        @testset "ZhouOptimalLAIModel construction for FT = $FT" begin
            params = Canopy.OptimalLAIParameters{FT}(toml_dict)
            model = Canopy.ZhouOptimalLAIModel{FT}(
                params;
                SAI = FT(0.0),
                RAI = FT(1.0),
                rooting_depth = FT(1.0),
                height = FT(10.0),
            )

            @test model.parameters === params
            @test eltype(model) == FT
            @test model.SAI == FT(0.0)
            @test model.RAI == FT(1.0)

            # Test auxiliary variables: the instantaneous potential GPP and χ, the
            # steady-state LAI target, and the growing-season inputs (GSL, vpd_gs,
            # f0), which are derived from the trailing climate totals in Y each step.
            aux_vars = Canopy.auxiliary_vars(model)
            @test :area_index in aux_vars
            @test :OptVars in aux_vars
            @test :L_opt in aux_vars
            @test :GSL in aux_vars
            @test :vpd_gs in aux_vars
            @test :f0 in aux_vars
            # the daily/annual accumulators are not cache variables; A0_daily,
            # A0_annual, and precip_annual are prognostic in Y
            @test :A0_daily_acc ∉ aux_vars
            @test :A0_annual ∉ aux_vars
            @test :precip_annual ∉ aux_vars

            # The optimal-LAI model carries the running-mean/running-sum climate
            # inputs (A0, precip, PET, VPD·A0, growing-days, per-pathway A0) plus the
            # prognostic LAI as time-integrated prognostic variables in Y. These are
            # what the climate-tracking f0 / vpd_gs / GSL / C3-C4 inputs derive from.
            optlai_prog = (
                :A0_daily,
                :A0_annual,
                :precip_annual,
                :PET_annual,
                :VPDA0_annual,
                :growing_days,
                :A0c3_annual,
                :A0c4_annual,
                :LAI,
            )
            @test Canopy.prognostic_vars(model) == optlai_prog
            @test Canopy.prognostic_types(model) ==
                  ntuple(_ -> FT, length(optlai_prog))
            @test Canopy.prognostic_domain_names(model) ==
                  ntuple(_ -> :surface, length(optlai_prog))
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

        @testset "c3_fraction_from_competition for FT = $FT" begin
            # Regression test for the dynamic C3/C4 competition (two bugs fixed in #1815).
            # The tree-cover proxy must use REALIZED GPP (a0c3·Mc·fapar), not potential
            # (fapar=1): potential GPP saturates the proxy and wrongly suppresses C4 in
            # sparse grasslands. So a sparser canopy (lower realized fapar) → less tree
            # suppression → MORE C4 → LOWER C3 fraction; using fapar=1 would invert this.
            Mc = FT(0.012)  # kg C per mol
            f = Canopy.c3_fraction_from_competition
            c3_sparse = f(FT(100), FT(130), Mc, FT(0.5))
            c3_dense = f(FT(100), FT(130), Mc, FT(1.0))
            @test c3_sparse < c3_dense
            # strong C3 GPP advantage → almost all C3
            @test f(FT(120), FT(40), Mc, FT(0.8)) > FT(0.9)
            # strong C4 advantage in a sparse canopy → almost all C4
            @test f(FT(40), FT(120), Mc, FT(0.3)) < FT(0.1)
            # fraction always in [0, 1]
            for args in (
                (FT(100), FT(130), Mc, FT(0.5)),
                (FT(120), FT(40), Mc, FT(0.8)),
                (FT(40), FT(120), Mc, FT(0.3)),
            )
                v = f(args...)
                @test FT(0) <= v <= FT(1)
            end
        end

        @testset "optimal_lai_initial_conditions for single-point domains for FT = $FT" begin
            # Test that optimal_lai_initial_conditions returns reasonable values
            # for single-point domains at various locations (Fluxnet sites)

            test_sites = [
                # (name, longitude, latitude)
                ("US-MOz (Ozark)", FT(-92.2000), FT(38.7441)),   # Missouri, USA - Deciduous forest
                ("US-Ha1 (Harvard)", FT(-72.1715), FT(42.5378)), # Massachusetts, USA - Mixed forest
                ("Amazon", FT(-60.0), FT(-3.0)),                 # Amazon rainforest
                # Semi-arid Africa. Longitude is kept off 0.0: `Point(longlat)` with
                # long == 0 (or lat == 0) currently builds a degenerate domain.
                ("Sahel", FT(2.5), FT(15.0)),
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
                    # GSL should be positive and reasonable. A year-round growing
                    # season gives the maximum, 12 * 365.25/12 = 365.25 days.
                    @test GSL_val >= FT(0)
                    @test GSL_val <= FT(366)
                    @test GSL_val > FT(0)

                    # A0_annual should be positive (mol CO2 m^-2 yr^-1)
                    # Typical values range from ~50 (arid) to ~500 (tropical rainforest)
                    @test A0_annual_val >= FT(0)
                    @test A0_annual_val > FT(0)
                    @test A0_annual_val < FT(1000)

                    # precip_annual should be positive (mol H2O m^-2 yr^-1)
                    # Ranges from ~5000 (desert, ~100 mm) to ~170000+ (tropical, ~3000 mm)
                    @test precip_annual_val >= FT(0)
                    @test precip_annual_val > FT(0)

                    # vpd_gs should be positive (Pa)
                    # Typical growing season VPD: 500-2500 Pa
                    @test vpd_gs_val >= FT(0)
                    @test vpd_gs_val > FT(0)

                    # lai_init should be non-negative (m^2 m^-2)
                    # Ranges from 0 (bare) to ~8 (dense forest)
                    @test lai_init_val >= FT(0)
                    @test lai_init_val < FT(15)

                    # f0 should be in range [0, 1]
                    @test f0_val >= FT(0)
                    @test f0_val <= FT(1)
                    @test f0_val > FT(0)
                end
            end
        end

        @testset "set_canopy_component_initial_conditions! for FT = $FT" begin
            # The prognostic state is set from the netCDF climatology at `set_ic!`
            # time, i.e. from a file rather than from the model.
            domain =
                Point(; z_sfc = FT(0.0), longlat = (FT(-92.2), FT(38.7441)))
            surface_space = domain.space.surface
            model = Canopy.ZhouOptimalLAIModel{FT}(
                domain,
                toml_dict;
                SAI = FT(0.0),
                RAI = FT(1.0),
                rooting_depth = FT(1.0),
                height = FT(10.0),
            )
            biomass_state = NamedTuple(
                var => ClimaCore.Fields.zeros(surface_space) for
                var in Canopy.prognostic_vars(model)
            )
            Y = ClimaCore.Fields.FieldVector(;
                canopy = (; biomass = biomass_state),
            )
            ClimaLand.Simulations.set_canopy_component_initial_conditions!(
                Y,
                nothing,
                model,
                nothing,
            )
            # the same climatology the IC reads, for the expected values
            ic = Canopy.optimal_lai_initial_conditions(surface_space)
            LAI = Array(parent(Y.canopy.biomass.LAI))[1]
            A0_annual = Array(parent(Y.canopy.biomass.A0_annual))[1]
            A0_daily = Array(parent(Y.canopy.biomass.A0_daily))[1]
            precip_annual = Array(parent(Y.canopy.biomass.precip_annual))[1]
            # LAI starts at the MODIS observation in the same file
            @test LAI ≈ Array(parent(ic.lai_init))[1]
            @test FT(0) < LAI < FT(15)
            # the annual totals start at their (steady-state) climatology, and the
            # one-day total at the corresponding daily share
            @test A0_annual ≈ Array(parent(ic.A0_annual))[1]
            @test A0_daily ≈ A0_annual / FT(365)
            @test precip_annual ≈ Array(parent(ic.precip_annual))[1]
        end
    end
end
