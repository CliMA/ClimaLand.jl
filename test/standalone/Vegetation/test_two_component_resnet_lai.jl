# test/standalone/Vegetation/test_two_component_resnet_lai.jl

using Test
using ClimaLand
using ClimaLand.Canopy

@testset "Two-Component ResNet LAI Model Tests" begin
    for FT in (Float32, Float64)
        @testset "TwoComponentResNetLAIParameters construction for FT = $FT" begin
            params = Canopy.TwoComponentResNetLAIParameters{FT}()

            @test params.k isa FT
            @test params.z isa FT
            @test params.sigma isa FT
            @test params.alpha_wood isa FT
            @test params.alpha_herb isa FT
            @test params.w1_a0 isa FT
            @test params.w1_vpd isa FT
            @test params.w1_sw isa FT
            @test params.b1 isa FT
            @test params.w2_h isa FT
            @test params.b2 isa FT
            @test params.w3_out isa FT
            @test params.b3_out isa FT
            @test params.autumn_damp isa FT
            @test params.enf_cold_gain isa FT
            @test params.pulse_thresh isa FT
            @test params.gamma_drought isa FT
            @test params.hydro_weight isa FT
            @test params.beta_dry isa FT
            @test params.green_min isa FT
            @test params.theta_freeze isa FT
            @test params.kappa_cold isa FT
            @test params.theta_ai_scale isa FT
            @test params.tau_thaw isa FT
            @test params.thaw_crit isa FT
            @test params.kappa_thaw isa FT
            @test params.enf_retention_scale isa FT

            @test params.k ≈ FT(0.5)
            @test params.z ≈ FT(12.227)
            @test params.alpha_wood ≈ FT(0.010)
            @test params.alpha_herb ≈ FT(0.59225)
            @test params.gamma_drought ≈ FT(0.5000)
            @test params.hydro_weight ≈ FT(0.1531)
            @test params.beta_dry ≈ FT(3.5000)
            @test params.green_min ≈ FT(0.6500)
            @test params.theta_freeze ≈ FT(273.15)
            @test params.kappa_cold ≈ FT(0.80)
            @test eltype(params) == FT
        end

        @testset "ResNet Dynamic Carbon Allocation for FT = $FT" begin
            params = Canopy.TwoComponentResNetLAIParameters{FT}()

            # Test nominal conditions
            alloc_nom = Canopy.compute_resnet_carbon_allocation(
                FT(0.5),
                FT(0.3),
                FT(0.6),
                params,
            )
            @test alloc_nom isa FT
            @test FT(0.60) <= alloc_nom <= FT(1.50)

            # Test sensitivity to GPP (A0): higher GPP -> higher allocation
            alloc_low_gpp = Canopy.compute_resnet_carbon_allocation(
                FT(0.1),
                FT(0.3),
                FT(0.6),
                params,
            )
            alloc_high_gpp = Canopy.compute_resnet_carbon_allocation(
                FT(0.9),
                FT(0.3),
                FT(0.6),
                params,
            )
            @test alloc_high_gpp > alloc_low_gpp

            # Test sensitivity to soil moisture (S_sw): wetter soil -> higher allocation
            alloc_dry_soil = Canopy.compute_resnet_carbon_allocation(
                FT(0.5),
                FT(0.3),
                FT(0.1),
                params,
            )
            alloc_wet_soil = Canopy.compute_resnet_carbon_allocation(
                FT(0.5),
                FT(0.3),
                FT(0.9),
                params,
            )
            @test alloc_wet_soil > alloc_dry_soil

            # Test sensitivity to atmospheric VPD: high VPD (drying) -> lower allocation
            alloc_low_vpd = Canopy.compute_resnet_carbon_allocation(
                FT(0.5),
                FT(0.1),
                FT(0.6),
                params,
            )
            alloc_high_vpd = Canopy.compute_resnet_carbon_allocation(
                FT(0.5),
                FT(0.9),
                FT(0.6),
                params,
            )
            @test alloc_low_vpd > alloc_high_vpd

            # Test hard bounds under extreme inputs
            alloc_extreme_low = Canopy.compute_resnet_carbon_allocation(
                -FT(100.0),
                FT(100.0),
                -FT(100.0),
                params,
            )
            @test alloc_extreme_low >= FT(0.60)

            alloc_extreme_high = Canopy.compute_resnet_carbon_allocation(
                FT(100.0),
                -FT(100.0),
                FT(100.0),
                params,
            )
            @test alloc_extreme_high <= FT(1.50)
        end

        @testset "Aridity and Seasonality Indices for FT = $FT" begin
            # Aridity index: wet/humid (1500 mm ~ 83,250 mol/m^2/yr) vs arid (300 mm ~ 16,650 mol/m^2/yr)
            ai_wet =
                Canopy.compute_aridity_index(FT(83250.0), FT(200.0), FT(800.0))
            ai_dry =
                Canopy.compute_aridity_index(FT(16650.0), FT(200.0), FT(2000.0))

            @test ai_wet isa FT
            @test ai_dry isa FT
            @test ai_wet > ai_dry
            @test FT(0.05) <= ai_dry <= FT(5.0)
            @test FT(0.05) <= ai_wet <= FT(5.0)

            # Seasonality index: tropical (GSL ~ 365) vs seasonal (GSL ~ 150)
            si_trop = Canopy.compute_seasonality_index(FT(365.25))
            si_temp = Canopy.compute_seasonality_index(FT(150.0))

            @test si_trop isa FT
            @test si_temp isa FT
            @test si_trop ≈ FT(0.0) atol = FT(1e-4)
            @test si_temp > si_trop
            @test FT(0.0) <= si_temp <= FT(1.0)
        end

        @testset "Two-Component Partitioning for FT = $FT" begin
            # Humid / low seasonality vs Arid / high seasonality
            f_wood_humid, f_herb_humid =
                Canopy.compute_two_component_partitioning(FT(2.5), FT(0.1))
            f_wood_arid, f_herb_arid =
                Canopy.compute_two_component_partitioning(FT(0.4), FT(0.9))

            @test f_wood_humid isa FT
            @test f_herb_humid isa FT
            @test f_wood_humid + f_herb_humid ≈ FT(1.0)
            @test f_wood_arid + f_herb_arid ≈ FT(1.0)
            @test f_wood_humid > f_wood_arid
            @test f_herb_arid > f_herb_humid
        end

        @testset "Topsoil Infiltration Pulse Hysteresis for FT = $FT" begin
            s_top = FT(0.0)
            thresh = FT(3.0)
            dt = FT(3600.0) # 1 hour

            # No precipitation
            s1, act1 = Canopy.compute_infiltration_pulse_hysteresis(
                s_top,
                FT(0.0),
                thresh,
                dt,
            )
            @test s1 isa FT
            @test act1 == FT(0.0)

            # Heavy precipitation event (10 mm/hr = 10 / 3600 mm/s)
            s2, act2 = Canopy.compute_infiltration_pulse_hysteresis(
                s_top,
                FT(10.0 / 3600.0),
                thresh,
                dt,
            )
            @test s2 > thresh
            @test act2 > FT(0.0)
            @test act2 <= FT(2.0)
        end

        @testset "Conifer Cold Buffer for FT = $FT" begin
            gain = FT(0.40)
            # Warm summer (+20 C = 293.15 K)
            buf_summer = Canopy.compute_conifer_cold_buffer(FT(293.15), gain)
            @test buf_summer ≈ FT(1.0) atol = FT(1e-2)

            # Cold winter (-20 C = 253.15 K)
            buf_winter = Canopy.compute_conifer_cold_buffer(FT(253.15), gain)
            @test buf_winter < buf_summer
            @test buf_winter ≈ FT(1.0) - gain atol = FT(0.05)
        end

        @testset "Two-Component Dynamics and Time Stepping for FT = $FT" begin
            l_wood = FT(2.0)
            l_herb = FT(0.5)
            l_steady_wood = FT(3.0)
            l_steady_herb = FT(1.5)
            alpha_w = FT(0.010)
            alpha_h = FT(0.592)
            dt = FT(1.0) # 1 day

            new_wood, new_herb = Canopy.compute_two_component_lai(
                l_wood,
                l_herb,
                l_steady_wood,
                l_steady_herb,
                alpha_w,
                alpha_h,
                dt,
            )

            @test new_wood isa FT
            @test new_herb isa FT
            # Both pools should increase towards targets
            @test new_wood > l_wood
            @test new_herb > l_herb
            # Herbaceous pool should adjust much faster than woody pool
            fractional_wood_change =
                (new_wood - l_wood) / (l_steady_wood - l_wood)
            fractional_herb_change =
                (new_herb - l_herb) / (l_steady_herb - l_herb)
            @test fractional_herb_change > fractional_wood_change
        end

        @testset "Cycle 86 Grand Hydro-Thermal Pure Functions for FT = $FT" begin
            # 1. Dryland Phenology Index
            # Wet rainforest: P = 2000 mm -> idx = 0
            idx_wet = Canopy.compute_dryland_phenology_index(FT(0.9), FT(0.0), FT(2000.0))
            # Dry savanna: P = 400 mm -> idx > 0
            idx_sav = Canopy.compute_dryland_phenology_index(FT(0.0), FT(0.0), FT(400.0))
            @test idx_wet isa FT
            @test idx_sav isa FT
            @test idx_wet ≈ FT(0.0)
            @test idx_sav ≈ FT(1.0 - 400.0 / 850.0) atol = FT(1e-4)
            @test FT(0.0) <= idx_sav <= FT(1.0)

            # 2. Thermal Freeze Gate
            # Warm summer: T = 298.15 K (+25 C) -> f_freeze ~ 1.0
            f_freeze_warm = Canopy.compute_thermal_freeze_gate(FT(298.15), FT(273.15), FT(0.80))
            # Deep winter: T = 253.15 K (-20 C) -> f_freeze ~ 0.0
            f_freeze_cold = Canopy.compute_thermal_freeze_gate(FT(253.15), FT(273.15), FT(0.80))
            @test f_freeze_warm isa FT
            @test f_freeze_cold isa FT
            @test f_freeze_warm > FT(0.95)
            @test f_freeze_cold < FT(0.05)

            # 3. Drought Dormancy Gate
            # Wet soil: soil_stress = 1.0 -> f_drought = 1.0
            f_drought_wet = Canopy.compute_drought_dormancy_gate(FT(1.0), idx_sav, FT(3.5))
            # Severe drought: soil_stress = 0.0 -> f_drought = 1 - idx_sav
            f_drought_dry = Canopy.compute_drought_dormancy_gate(FT(0.0), idx_sav, FT(3.5))
            @test f_drought_wet isa FT
            @test f_drought_dry isa FT
            @test f_drought_wet ≈ FT(1.0)
            @test f_drought_dry < f_drought_wet

            # 4. Dual-Trigger Multiplicative Coupling
            f_dormant = Canopy.compute_dual_trigger_dormancy(f_freeze_warm, f_drought_dry)
            @test f_dormant isa FT
            @test f_dormant ≈ f_freeze_warm * f_drought_dry

            # 5. Hydro-Optical Greenness Gating
            f_green_nom = Canopy.compute_hydro_optical_greenness(
                FT(0.0), FT(0.0), idx_sav, FT(1.0), FT(1.0), FT(1.0), FT(1.0), FT(0.65), FT(0.1531)
            )
            f_green_cured = Canopy.compute_hydro_optical_greenness(
                FT(0.0), FT(0.0), idx_sav, FT(1.0), FT(1.0), FT(1.0), FT(0.0), FT(0.65), FT(0.1531)
            )
            @test f_green_nom isa FT
            @test f_green_cured isa FT
            @test f_green_nom ≈ FT(1.0)
            @test f_green_cured < f_green_nom
            @test f_green_cured >= FT(0.65) * (FT(1.0) - FT(0.1531))

            # 6. Drought-Accelerated Woody Turnover
            alpha_w_base = FT(0.010)
            alpha_w_wet = Canopy.compute_drought_accelerated_woody_turnover(alpha_w_base, FT(0.5000), idx_sav, FT(1.0))
            alpha_w_drought = Canopy.compute_drought_accelerated_woody_turnover(alpha_w_base, FT(0.5000), idx_sav, FT(0.0))
            @test alpha_w_wet isa FT
            @test alpha_w_drought isa FT
            @test alpha_w_wet ≈ alpha_w_base
            @test alpha_w_drought > alpha_w_wet

            # 7. Structural Canopy Floor
            floor_rf = Canopy.compute_structural_canopy_floor(
                FT(1.0), FT(0.0), FT(0.0), FT(3.0), FT(1.0), FT(1.8), FT(0.15), FT(0.39384), FT(298.15)
            )
            floor_decid_dormant = Canopy.compute_structural_canopy_floor(
                FT(0.0), FT(0.0), FT(1.0), FT(1.0), FT(0.0), FT(1.8), FT(0.15), FT(0.39384), FT(298.15)
            )
            @test floor_rf isa FT
            @test floor_decid_dormant isa FT
            @test floor_rf >= FT(2.05)
            @test floor_decid_dormant ≈ FT(0.05)
        end

        @testset "TwoComponentResNetLAIModel struct construction for FT = $FT" begin
            params = Canopy.TwoComponentResNetLAIParameters{FT}()
            inputs = (;
                GSL = FT(240.0),
                A0_annual = FT(258.0),
                precip_annual = FT(1000.0),
                vpd_gs = FT(1000.0),
                lai_init = FT(2.0),
                f0 = FT(0.65),
            )
            model = Canopy.TwoComponentResNetLAIModel{FT}(
                params,
                inputs;
                SAI = FT(0.0),
                RAI = FT(1.0),
                rooting_depth = FT(1.0),
                height = FT(10.0),
            )

            @test model.parameters === params
            @test model.optimal_lai_inputs === inputs
            @test eltype(model) == FT
            @test model.SAI == FT(0.0)
            @test model.RAI == FT(1.0)
            @test model.rooting_depth == FT(1.0)
            @test model.height == FT(10.0)

            # Test auxiliary variables
            aux_vars = Canopy.auxiliary_vars(model)
            @test :area_index in aux_vars
            @test :OptVars in aux_vars
            @test :L_opt_wood in aux_vars
            @test :L_opt_herb in aux_vars

            # Test prognostic variables
            prog_vars = Canopy.prognostic_vars(model)
            @test :A0_daily in prog_vars
            @test :A0_annual in prog_vars
            @test :precip_annual in prog_vars
            @test :LAI_wood in prog_vars
            @test :LAI_herb in prog_vars
            @test Canopy.prognostic_types(model) == (FT, FT, FT, FT, FT)
        end
    end
end
