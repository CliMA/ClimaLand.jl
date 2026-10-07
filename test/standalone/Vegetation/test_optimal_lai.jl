using Test
using ClimaLand
import ClimaComms
ClimaComms.@import_required_backends
using ClimaLand.Canopy
using ClimaLand.Domains: Point
import ClimaLand.Parameters as LP
import ClimaParams
using ClimaCore
import NCDatasets

@testset "Optimal LAI Model Tests" begin
    for FT in (Float32, Float64)
        toml_dict = LP.create_toml_dict(FT)

        @testset "OptimalLAIParameters construction for FT = $FT" begin
            # Test parameter construction from TOML
            params = Canopy.OptimalLAIParameters{FT}(toml_dict)

            @test params.k isa FT
            @test params.z_tree isa FT
            @test params.z_grass isa FT
            @test params.sigma isa FT
            @test params.alpha isa FT
            @test params.tau_long_term isa FT

            # Check expected values from default_parameters.toml, calibrated against
            # MODIS LAI
            @test params.k ≈ FT(0.5)
            @test params.z_tree ≈ FT(9.92)
            @test params.z_grass ≈ FT(154)
            @test params.sigma ≈ FT(1.09)
            @test params.alpha ≈ FT(0.202)  # ~15 days of memory
            @test params.f0_max ≈ FT(0.65)
            @test params.tau_long_term ≈ FT(6.3072e7)  # 2 years

            # Lavergne et al. (2022) C3/C4 logistic and tree-cover coefficients
            @test params.c3c4_k ≈ FT(6.63)
            @test params.c3c4_q ≈ FT(0.16)
            @test params.tc_a ≈ FT(15.60)
            @test params.tc_b ≈ FT(1.41)
            @test params.tc_c ≈ FT(-7.72)
            @test params.tc_gpp_ref ≈ FT(2.8)

            # The model's own P-model unit cost ratios (pyrealm defaults)
            @test params.β_c3 ≈ FT(146)
            @test params.β_c4 ≈ FT(146 / 9) rtol = 1e-4
            pmodel = Canopy.PModel{FT}(
                ClimaLand.Domains.Point(;
                    z_sfc = FT(0),
                    longlat = FT.((-60, -3)),
                ),
                toml_dict,
            )
            lai_pmodel =
                Canopy.optimal_lai_pmodel_parameters(pmodel.parameters, params)
            @test lai_pmodel.β_c3 == params.β_c3
            @test lai_pmodel.β_c4 == params.β_c4
            @test lai_pmodel.cstar == pmodel.parameters.cstar

            # logistic of the climate tree share, and the canopy trees retain
            @test params.tree_b0 ≈ FT(-1.07)
            @test params.tree_b_lai ≈ FT(0.0254)
            @test params.tree_b_dry ≈ FT(-0.372)
            @test params.tree_b_temp ≈ FT(0.0917)
            @test params.tree_retention ≈ FT(0.5)
            # season-scaled leaf cost
            @test params.maintenance_share ≈ FT(0.43)
            @test params.maintenance_q10 ≈ FT(2)

            @test eltype(params) == FT
        end

        @testset "climate_tree_share for FT = $FT" begin
            params = Canopy.OptimalLAIParameters{FT}(toml_dict)
            # a humid tropical forest, a desert, a boreal forest
            forest = Canopy.climate_tree_share(FT(5.8), FT(0), FT(23), params)
            desert = Canopy.climate_tree_share(FT(0), FT(12), FT(28), params)
            boreal = Canopy.climate_tree_share(FT(3.5), FT(0), FT(9), params)
            @test forest isa FT
            @test forest > FT(0.7)
            @test desert < FT(0.1)
            @test FT(0) < desert < boreal < forest < FT(1)
            # a longer dry season favours grasses
            @test Canopy.climate_tree_share(FT(3), FT(6), FT(20), params) <
                  Canopy.climate_tree_share(FT(3), FT(2), FT(20), params)
            @test Canopy.climate_tree_share(FT(0), FT(0), FT(0), params) ≈
                  1 / (1 + exp(-params.tree_b0))
        end

        @testset "optimal_chi for FT = $FT" begin
            pmodel = Canopy.PModel{FT}(
                ClimaLand.Domains.Point(;
                    z_sfc = FT(0),
                    longlat = FT.((-60, -3)),
                ),
                toml_dict,
            )
            params = Canopy.OptimalLAIParameters{FT}(toml_dict)
            lai_pmodel =
                Canopy.optimal_lai_pmodel_parameters(pmodel.parameters, params)
            χ(T, vpd, β = lai_pmodel.β_c3) = Canopy.optimal_chi(
                FT(T),
                FT(101325),
                FT(4.2e-4),
                FT(vpd),
                β,
                pmodel.constants,
            )
            @test χ(298, 1000) isa FT
            @test FT(0) < χ(298, 1000) < FT(1)
            # stomata close as the air dries
            @test χ(298, 2000) < χ(298, 1000) < χ(298, 500)
            # a larger cost ratio β keeps the stomata more open; C4 plants are
            # more conservative
            @test χ(298, 1000) > χ(298, 1000, pmodel.parameters.β_c3)
            @test χ(298, 1000, lai_pmodel.β_c4) < χ(298, 1000)
        end

        @testset "moist_season_vpd for FT = $FT" begin
            day = FT(86400)
            vpd(moist_days) = Canopy.moist_season_vpd(
                FT(500) * moist_days * day,
                FT(moist_days),
                FT(2000) * 300 * day,
                FT(300),
                day,
            )
            # the VPD of the moist season, once it lasts a month
            @test vpd(200) ≈ 500
            @test vpd(30) ≈ 500
            # without a moist season, that of the whole growing season
            @test vpd(0) ≈ 2000
            @test vpd(0) > vpd(15) > vpd(30)
        end

        @testset "growing_running_sum_tendency for FT = $FT" begin
            day = FT(86400)
            reduction = ClimaLand.RunningSum(365 * day, 2 * 365 * day)
            # integrate a constant rate from zero with explicit Euler, as the model does
            X, age, Δt, f = FT(0), FT(0), FT(900), 1 / day
            for _ in 1:(30 * 96)
                X +=
                    Δt *
                    Canopy.growing_running_sum_tendency(f, X, age, reduction)
                age += Δt
            end
            # after a month the total is already close to its 365-day steady state
            @test X ≈ 365 rtol = 0.05
            # once older than τ_long, it is a RunningSum
            @test Canopy.growing_running_sum_tendency(
                f,
                FT(100),
                FT(1e9),
                reduction,
            ) ≈ ClimaLand.apply_time_reduction(f, FT(100), reduction)
        end

        @testset "maintenance_rate and leaf_cost_scale for FT = $FT" begin
            T_freeze = FT(273.15)
            @test Canopy.maintenance_rate(T_freeze + 25, T_freeze, FT(2)) ≈ 1
            @test Canopy.maintenance_rate(T_freeze + 15, T_freeze, FT(2)) ≈
                  FT(0.5)
            @test Canopy.maintenance_rate(T_freeze - 5, T_freeze, FT(2)) == 0
            # a year-round 25 °C season keeps the full cost; shorter and colder
            # seasons are cheaper, down to the construction share
            m = FT(0.43)
            @test Canopy.leaf_cost_scale(FT(365), m) ≈ 1
            @test Canopy.leaf_cost_scale(FT(0), m) ≈ 1 - m
            @test Canopy.leaf_cost_scale(FT(50), m) <
                  Canopy.leaf_cost_scale(FT(200), m)
        end

        @testset "leaf_cost for FT = $FT" begin
            z_tree, z_grass = FT(12), FT(100)
            @test Canopy.leaf_cost(FT(1), z_tree, z_grass) ≈ z_tree
            @test Canopy.leaf_cost(FT(0), z_tree, z_grass) ≈ z_grass
            # geometric mean of the two costs at an even share
            @test Canopy.leaf_cost(FT(0.5), z_tree, z_grass) ≈
                  sqrt(z_tree * z_grass)
            @test Canopy.leaf_cost(FT(0.8), z_tree, z_grass) <
                  Canopy.leaf_cost(FT(0.2), z_tree, z_grass)
        end

        @testset "ZhouOptimalLAIModel construction for FT = $FT" begin
            params = Canopy.OptimalLAIParameters{FT}(toml_dict)
            model = Canopy.ZhouOptimalLAIModel{FT}(
                params;
                SAI = FT(0.0),
                RAI = FT(1.0),
                rooting_depth = FT(1.0),
                height = FT(10.0),
                tree_share = FT(0.7),
            )

            @test model.parameters === params
            @test eltype(model) == FT
            @test model.SAI == FT(0.0)
            @test model.RAI == FT(1.0)
            @test model.tree_share == FT(0.7)

            # Cache variables: the potential GPP and χ, the steady-state LAI target,
            # and the growing-season inputs derived from the totals in Y.
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

            # The trailing climate totals and LAI are time-integrated prognostic
            # variables in Y.
            optlai_prog = (
                :A0_daily,
                :A0_annual,
                :precip_annual,
                :PET_annual,
                :VPDgs_annual,
                :growing_days,
                :A0c3_annual,
                :A0c4_annual,
                :GPPc3_annual,
                :LAI,
                :snow_store,
                :water_30d,
                :PET_30d,
                :VPD_moist_annual,
                :moist_days,
                :degree_days,
                :warm_days,
                :maintenance_days,
                :age,
            )
            @test Canopy.prognostic_vars(model) == optlai_prog
            @test Canopy.prognostic_types(model) ==
                  ntuple(_ -> FT, length(optlai_prog))
            @test Canopy.prognostic_domain_names(model) ==
                  ntuple(_ -> :surface, length(optlai_prog))

            # A prognostic tree share uses the same climate totals
            prognostic_tree = Canopy.ZhouOptimalLAIModel{FT}(
                params;
                SAI = FT(0.0),
                RAI = FT(1.0),
                rooting_depth = FT(1.0),
                height = FT(10.0),
                tree_share = Canopy.PrognosticTreeShare(),
            )
            @test Canopy.prognostic_vars(prognostic_tree) == optlai_prog
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
            # The tree-cover term uses the realized C3 GPP (a0c3·fapar), so a
            # sparser canopy means less tree shading, more C4 and a lower C3 fraction.
            Mc = FT(0.012)  # kg C per mol
            params = Canopy.OptimalLAIParameters{FT}(toml_dict)
            f =
                (a3, a4, fapar) -> Canopy.c3_fraction_from_competition(
                    a3,
                    a4,
                    a3 * fapar,
                    FT(365),
                    Mc,
                    params,
                )
            c3_sparse = f(FT(100), FT(130), FT(0.5))
            c3_dense = f(FT(100), FT(130), FT(1.0))
            @test c3_sparse < c3_dense
            # strong C3 GPP advantage → almost all C3
            @test f(FT(120), FT(40), FT(0.8)) > FT(0.9)
            # strong C4 advantage in a sparse canopy → almost all C4
            @test f(FT(40), FT(120), FT(0.3)) < FT(0.1)
            # fraction always in [0, 1]
            for args in (
                (FT(100), FT(130), FT(0.5)),
                (FT(120), FT(40), FT(0.8)),
                (FT(40), FT(120), FT(0.3)),
            )
                v = f(args...)
                @test FT(0) <= v <= FT(1)
            end
        end

        @testset "canopy_composition_from_competition for FT = $FT" begin
            Mc = FT(0.012)  # kg C per mol
            params = Canopy.OptimalLAIParameters{FT}(toml_dict)
            g =
                (a3, a4, fapar) -> Canopy.canopy_composition_from_competition(
                    a3,
                    a4,
                    a3 * fapar,
                    FT(365),
                    Mc,
                    params,
                )
            for args in (
                (FT(100), FT(130), FT(0.5)),
                (FT(120), FT(40), FT(0.8)),
                (FT(40), FT(120), FT(0.3)),
                (FT(400), FT(500), FT(1.0)),
            )
                c = g(args...)
                @test all(FT(0) .<= (c.tree, c.c3_grass, c.c4_grass) .<= FT(1))
                @test c.tree + c.c3_grass + c.c4_grass ≈ FT(1)
                # trees are all C3
                a3, a4, fapar = args
                @test Canopy.c3_fraction_from_competition(
                    a3,
                    a4,
                    a3 * fapar,
                    FT(365),
                    Mc,
                    params,
                ) ≈ c.tree + c.c3_grass
            end
            # below the tree-cover threshold GPP (≈0.6 kg C m^-2 yr^-1) there
            # are no trees, so the open-canopy split is the whole canopy
            sparse = g(FT(40), FT(120), FT(0.3))
            @test sparse.tree == FT(0)
            @test sparse.c4_grass > sparse.c3_grass
            # above canopy closure (tc_gpp_ref = 2.8 kg C m^-2 yr^-1) everything is
            # trees: C4 grasses are shaded out whatever their advantage
            closed = g(FT(400), FT(500), FT(1.0))
            @test closed.tree == FT(1)
            @test closed.c4_grass == FT(0)
            # a sparser canopy lowers the tree share
            @test g(FT(150), FT(150), FT(0.4)).tree <
                  g(FT(150), FT(150), FT(0.9)).tree
            # the same annual GPP over a shorter season, as in boreal forests,
            # supports more trees
            gppc3 = FT(70)  # ≈ 0.84 kg C m^-2 yr^-1
            @test Canopy.tree_share_from_gpp(gppc3, FT(365), Mc, params) <
                  FT(0.1) <
                  Canopy.tree_share_from_gpp(gppc3, FT(150), Mc, params)
        end

        @testset "c4_advantage_for_c3_fraction for FT = $FT" begin
            # The inverse reproduces a C3 fraction the competition can reach.
            Mc = FT(0.012)
            params = Canopy.OptimalLAIParameters{FT}(toml_dict)
            for (a3, fapar, fc3) in (
                (FT(100), FT(0.5), FT(0.3)),
                (FT(150), FT(0.9), FT(0.8)),
                (FT(40), FT(0.3), FT(0.5)),
            )
                gppc3 = a3 * fapar
                tree = Canopy.tree_share_from_gpp(gppc3, FT(365), Mc, params)
                adv = Canopy.c4_advantage_for_c3_fraction(fc3, tree, params)
                @test Canopy.c3_fraction_from_competition(
                    a3,
                    a3 * (1 + adv),
                    gppc3,
                    FT(365),
                    Mc,
                    params,
                ) ≈ fc3 rtol = 1e-3
            end
            # A pure-C3 map needs A0c4 = 0; a small C4 grass share remains.
            adv = Canopy.c4_advantage_for_c3_fraction(FT(1), FT(0), params)
            @test adv == -1
            @test FT(0.9) <
                  Canopy.c3_fraction_from_competition(
                      FT(100),
                      FT(0),
                      FT(0),
                      FT(365),
                      Mc,
                      params,
                  ) <
                  FT(1)
            # More C4 than the open canopy allows saturates at the open canopy.
            gppc3 = FT(150) * FT(0.9)
            tree = Canopy.tree_share_from_gpp(gppc3, FT(365), Mc, params)
            @test FT(0) < tree < FT(1)
            adv = Canopy.c4_advantage_for_c3_fraction(FT(0), tree, params)
            @test isfinite(adv)
            @test Canopy.c3_fraction_from_competition(
                FT(150),
                FT(150) * (1 + adv),
                gppc3,
                FT(365),
                Mc,
                params,
            ) ≈ tree atol = 1e-3
        end

        @testset "f0_from_aridity / aridity_from_f0 for FT = $FT" begin
            f0_max = FT(0.65)
            # f0 peaks at f0_max at the energy-water transition and falls off on
            # both sides, so it is not invertible without choosing a branch.
            @test Canopy.f0_from_aridity(FT(1.9), FT(1), f0_max) ≈ f0_max
            @test Canopy.f0_from_aridity(FT(19), FT(1), f0_max) < f0_max
            @test Canopy.f0_from_aridity(FT(0.19), FT(1), f0_max) < f0_max

            # aridity_from_f0 is the arid-branch inverse, which is what the initial
            # conditions rely on to seed PET so the online f0 starts at the map value.
            for f0 in (FT(0.1), FT(0.3), FT(0.6), f0_max)
                AI = Canopy.aridity_from_f0(f0, f0_max)
                @test AI >= FT(1.9)
                @test Canopy.f0_from_aridity(AI, FT(1), f0_max) ≈ f0 rtol = 1e-5
            end
            # an f0 above the peak has no preimage; the inverse saturates at the peak
            @test Canopy.aridity_from_f0(FT(0.9), f0_max) ≈ FT(1.9)
        end

        @testset "potential_evaporation for FT = $FT" begin
            earth_param_set = LP.LandParameters(toml_dict)
            thermo_params = LP.thermodynamic_parameters(earth_param_set)
            σ = LP.Stefan(earth_param_set)
            M_w = LP.molar_mass_water(earth_param_set)
            λv = LP.LH_v0(earth_param_set)
            ϵ = FT(0.98)
            P = FT(101325)
            pet(; SW = 220, LW = 330, T = 293.15) =
                Canopy.potential_evaporation(
                    FT(SW),
                    FT(LW),
                    FT(T),
                    P,
                    ϵ,
                    σ,
                    M_w,
                    thermo_params,
                )
            mm_per_day(x) = x * M_w * FT(86400)
            Rn_over_λ(; SW = 220, LW = 330, T = 293.15) =
                ((1 - FT(0.17)) * FT(SW) + ϵ * (FT(LW) - σ * FT(T)^4)) /
                (λv * M_w)

            @test pet() > FT(0)
            @test isfinite(pet())
            # Daily-mean forcing of a temperate summer gives a few mm/day.
            @test FT(1) < mm_per_day(pet()) < FT(8)

            # Priestley-Taylor: 1.26 Δ/(Δ+γ) of the net radiation, so above the
            # equilibrium evaporation Δ/(Δ+γ)·Rn but below 1.26·Rn ...
            @test pet() < FT(1.26) * Rn_over_λ()
            @test pet() > FT(0.5) * Rn_over_λ()
            # ... and a larger fraction of Rn in the warm, where Δ/(Δ+γ) is larger.
            @test pet(T = 278.15) / Rn_over_λ(T = 278.15) <
                  pet(T = 303.15) / Rn_over_λ(T = 303.15)

            # Negative (night-time) net radiation does not count.
            @test pet(SW = 0, LW = 250) == FT(0)
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
                    ClimaLand.Simulations.optimal_lai_initial_conditions(
                        surface_space,
                    )

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
            # the same climatology the IC reads, for the expected values
            ic = ClimaLand.Simulations.optimal_lai_initial_conditions(
                surface_space,
            )
            scalar(field) = Array(parent(field))[1]
            max_lai_field = Canopy.modis_max_lai(surface_space)
            # A static map (standing in for the photosynthesis model's) with half
            # of the open canopy C4, which the competition can reach
            Mc = FT(0.0120107)
            k = model.parameters.k
            tree = Canopy.tree_share_from_gpp(
                scalar(ic.A0_annual) * (1 - exp(-k * scalar(max_lai_field))),
                scalar(ic.GSL),
                Mc,
                model.parameters,
            )
            fractional_c3 = tree + (1 - tree) / 2
            ClimaLand.Simulations.set_canopy_component_initial_conditions!(
                Y,
                nothing,
                model,
                nothing,
                ClimaLand.Artifacts.optimal_lai_initial_conditions_path(),
                max_lai_field,
                fractional_c3,
                Mc,
                nothing,
            )
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

            # The climate-responsive accumulators are seeded so each online input
            # reproduces the map value it replaces at t = 0.
            PET_annual = scalar(Y.canopy.biomass.PET_annual)
            f0_max = model.parameters.f0_max
            @test Canopy.f0_from_aridity(PET_annual, precip_annual, f0_max) ≈
                  scalar(ic.f0) rtol = 1e-5
            @test scalar(Y.canopy.biomass.VPDgs_annual) /
                  (scalar(Y.canopy.biomass.growing_days) * FT(86400)) ≈
                  scalar(ic.vpd_gs) rtol = 1e-5
            @test scalar(Y.canopy.biomass.growing_days) ≈ scalar(ic.GSL)
            @test scalar(Y.canopy.biomass.A0c3_annual) ≈ A0_annual
            # the realized C3 GPP is seeded with the fAPAR of the MODIS max LAI,
            # not of the (winter) lai_init snapshot
            max_lai = scalar(max_lai_field)
            GPPc3_annual = scalar(Y.canopy.biomass.GPPc3_annual)
            @test GPPc3_annual ≈ A0_annual * (1 - exp(-k * max_lai))
            @test GPPc3_annual > A0_annual * (1 - exp(-k * LAI))
            # A0c4 is seeded so the competition starts at the static C3 map
            @test Canopy.c3_fraction_from_competition(
                scalar(Y.canopy.biomass.A0c3_annual),
                scalar(Y.canopy.biomass.A0c4_annual),
                GPPc3_annual,
                scalar(Y.canopy.biomass.growing_days),
                Mc,
                model.parameters,
            ) ≈ fractional_c3 rtol = 1e-3
        end

        @testset "set_canopy_component_initial_conditions! with missing data for FT = $FT" begin
            # The file has no f0 south of 60°S, and no data at all over ocean
            for (longlat, all_missing) in
                ((FT(10), FT(-80)) => false, (FT(-150), FT(-10)) => true)
                domain = Point(; z_sfc = FT(0.0), longlat)
                surface_space = domain.space.surface
                model = Canopy.ZhouOptimalLAIModel{FT}(
                    domain,
                    toml_dict;
                    SAI = FT(0.0),
                    RAI = FT(1.0),
                    rooting_depth = FT(1.0),
                    height = FT(1.0),
                )
                Y = ClimaCore.Fields.FieldVector(;
                    canopy = (;
                        biomass = NamedTuple(
                            var => ClimaCore.Fields.zeros(surface_space) for
                            var in Canopy.prognostic_vars(model)
                        )
                    ),
                )
                ic = ClimaLand.Simulations.optimal_lai_initial_conditions(
                    surface_space,
                )
                scalar(field) = Array(parent(field))[1]
                @test isnan(scalar(ic.f0))
                @test isnan(scalar(ic.precip_annual)) == all_missing
                ClimaLand.Simulations.set_canopy_component_initial_conditions!(
                    Y,
                    nothing,
                    model,
                    nothing,
                    ClimaLand.Artifacts.optimal_lai_initial_conditions_path(),
                    Canopy.modis_max_lai(surface_space),
                    FT(1),
                    FT(0.0120107),
                    nothing,
                )
                for var in Canopy.prognostic_vars(model)
                    @test isfinite(scalar(getproperty(Y.canopy.biomass, var)))
                end
                # f0 starts at the peak of f0(AI)
                @test scalar(Y.canopy.biomass.PET_annual) ≈
                      FT(1.9) * scalar(Y.canopy.biomass.precip_annual)
                if all_missing
                    @test scalar(Y.canopy.biomass.LAI) == 0
                    @test scalar(Y.canopy.biomass.A0_annual) == 0
                end
            end
        end

        @testset "set_canopy_component_initial_conditions! from the spin-up state for FT = $FT" begin
            state_path = ClimaLand.Artifacts.optimal_lai_state_path()
            # a forest in Missouri, where the state has data, and the ocean, where
            # the climatology is kept
            for (longlat, from_state) in
                ((FT(-92.2), FT(38.7441)) => true, (FT(-150), FT(-10)) => false)
                domain = Point(; z_sfc = FT(0.0), longlat)
                surface_space = domain.space.surface
                model = Canopy.ZhouOptimalLAIModel{FT}(
                    domain,
                    toml_dict;
                    SAI = FT(0.0),
                    RAI = FT(1.0),
                    rooting_depth = FT(1.0),
                    height = FT(1.0),
                    tree_share = Canopy.PrognosticTreeShare(),
                )
                new_state() = ClimaCore.Fields.FieldVector(;
                    canopy = (;
                        biomass = NamedTuple(
                            var => ClimaCore.Fields.zeros(surface_space) for
                            var in Canopy.prognostic_vars(model)
                        )
                    ),
                )
                scalar(field) = Array(parent(field))[1]
                set_ic!(Y, path) =
                    ClimaLand.Simulations.set_canopy_component_initial_conditions!(
                        Y,
                        nothing,
                        model,
                        nothing,
                        ClimaLand.Artifacts.optimal_lai_initial_conditions_path(),
                        Canopy.modis_max_lai(surface_space),
                        FT(1),
                        FT(0.0120107),
                        path,
                    )
                Y = new_state()
                set_ic!(Y, state_path)
                Y_climatology = new_state()
                set_ic!(Y_climatology, nothing)
                # the nearest cell of the file
                ds = NCDatasets.NCDataset(state_path)
                i = argmin(abs.(ds["lon"][:] .- longlat[1]))
                j = argmin(abs.(ds["lat"][:] .- longlat[2]))
                for name in ClimaLand.Simulations.OPTIMAL_LAI_STATE_NAMES
                    value = scalar(getproperty(Y.canopy.biomass, name))
                    state = ds[String(name)][i, j]
                    @test isfinite(value)
                    climatology =
                        scalar(getproperty(Y_climatology.canopy.biomass, name))
                    @test (value == climatology) == !from_state
                    @test isnan(state) == !from_state
                    from_state && @test value ≈ state
                end
                close(ds)
                if from_state
                    # PET is that of a humid climate, not the arid-branch seed
                    @test scalar(Y.canopy.biomass.PET_annual) <
                          scalar(Y_climatology.canopy.biomass.PET_annual)
                end
                # the 30-day totals start from the annual ones, the yearly sums of
                # the moist and warm seasons from zero
                @test scalar(Y.canopy.biomass.water_30d) ≈
                      scalar(Y.canopy.biomass.precip_annual) * 30 / 365
                @test scalar(Y.canopy.biomass.PET_30d) ≈
                      scalar(Y.canopy.biomass.PET_annual) * 30 / 365
                for name in (
                    :snow_store,
                    :VPD_moist_annual,
                    :moist_days,
                    :degree_days,
                    :warm_days,
                    :maintenance_days,
                    :age,
                )
                    @test scalar(getproperty(Y.canopy.biomass, name)) == 0
                end
            end
        end
    end
end
