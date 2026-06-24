using Test
import ClimaParams as CP
using ClimaLand.Snow
import ClimaLand
import ClimaLand.Parameters as LP
import Random
using Thermodynamics, SurfaceFluxes

Random.seed!(1234)

for FT in (Float32, Float64)
    toml_dict = LP.create_toml_dict(FT)
    @testset "Snow Parameterizations, FT = $FT" begin
        param_set = LP.LandParameters(toml_dict)
        thermo_params = LP.thermodynamic_parameters(param_set)

        # Density of liquid water (kg/m``^3``)
        _ρ_l = FT(LP.ρ_cloud_liq(param_set))
        # Density of ice water (kg/m``^3``)
        _ρ_i = FT(LP.ρ_cloud_ice(param_set))
        _cp_l = FT(LP.cp_l(param_set))
        _cp_i = FT(LP.cp_i(param_set))
        # Reference temperature (K)
        _T_ref = FT(LP.T_0(param_set))
        # Freezing temperature (K)
        _T_freeze = FT(LP.T_freeze(param_set))
        # Latent heat of fusion at ``T_0`` (J/kg)
        _LH_f0 = FT(LP.LH_f0(param_set))
        # Thermal conductivity of dry air
        κ_air = FT(LP.K_therm(param_set))

        ρ_min = FT(200)
        densitymodel = Snow.MinimumDensityModel(ρ_min)
        z_0m = FT(0.0024)
        α_snow = Snow.ConstantAlbedoModel(FT(0.8))
        scf = Snow.WuWuSnowCoverFractionModel(
            FT(0.106),
            FT(1.81),
            FT(0.08),
            FT(1.77),
            FT(1),
            FT(1),
        )

        # These values should match ClimaParams
        ϵ_snow = FT(0.97)
        z_0b = FT(8e-2)
        θ_r = FT(0.08)
        Ksat = FT(1e-3)
        ΔS = FT(0.1)
        Δt = Float64(180.0)
        parameters = SnowParameters(
            toml_dict,
            Δt,
            density = densitymodel,
            α_snow = α_snow,
            surf_temp =Snow.EquilibriumGradientTemperatureModel{FT}(),
        )
        @test parameters.density.ρ_min == ρ_min
        @test typeof(parameters.density.ρ_min) == FT
        @test parameters.ΔS == ΔS
        @test typeof(parameters.ΔS) == FT
        @test parameters.z_0m == z_0m
        @test typeof(parameters.z_0m) == FT
        @test parameters.z_0b == z_0b
        @test typeof(parameters.z_0b) == FT
        @test parameters.α_snow == α_snow
        @test typeof(parameters.α_snow) == Snow.ConstantAlbedoModel{FT}
        @test parameters.scf == scf
        @test typeof(parameters.scf) == Snow.WuWuSnowCoverFractionModel{FT}
        @test parameters.ϵ_snow == ϵ_snow
        @test typeof(parameters.ϵ_snow) == FT
        @test parameters.Ksat == Ksat
        @test typeof(parameters.Ksat) == FT
        @test parameters.θ_r == θ_r
        @test typeof(parameters.θ_r) == FT
        @test typeof(parameters.κ_snow) == Snow.SturmSnowConductivityModel{FT}

        roughness_model = SurfaceFluxes.ConstantRoughnessParams{FT}(z_0m, z_0b)

        T = FT.([275.0, 272, _T_freeze])
        @test specific_heat_capacity(FT(1.0), param_set) == _cp_l
        @test specific_heat_capacity(FT(0.0), param_set) == _cp_i

        sturm = parameters.κ_snow
        ρ_snow = ρ_min
        x = min(ρ_snow / _ρ_i, sturm.max / _ρ_i)
        @test snow_thermal_conductivity(parameters.κ_snow, ρ_snow, param_set) ==
              sturm.b2 + sturm.m2 * x + sturm.q2 * x^2
        @test typeof(sturm) == Snow.SturmSnowConductivityModel{FT}
        # densities below the threshold, above the threshold, and above the
        # maximum density (where the non-dimensional density is capped)
        for ρ_test in FT.([100, 350, 800])
            x = min(ρ_test / _ρ_i, sturm.max / _ρ_i)
            expected =
                x < sturm.threshold ? sturm.b1 + sturm.m1 * x :
                sturm.b2 + sturm.m2 * x + sturm.q2 * x^2
            @test snow_thermal_conductivity(sturm, ρ_test, param_set) ==
                  expected
        end

        @test all(
            maximum_liquid_mass_fraction.(ρ_min, T, θ_r, Ref(param_set)) .-
            [FT(0), θ_r * _ρ_l / ρ_snow, θ_r * _ρ_l / ρ_snow] .== 0,
        )

        SWE = cat(FT.(rand(10)), FT(0), dims = 1)
        z = SWE * _ρ_l ./ ρ_snow
        @test runoff_timescale.(z, Ksat, FT(Δt)) ≈ max.(Δt, z ./ Ksat)
        ρ_calc = Snow.snow_bulk_density.(SWE, z, param_set)
        @test all(ρ_calc[1:(end - 1)] .≈ ρ_snow)
        @test ρ_calc[end] == _ρ_l
        @test Snow.snow_bulk_density(eps(FT(0)), 2 * eps(FT(0)), param_set) ==
              _ρ_l

        U = energy_from_q_l_and_swe(FT(1), FT(0.5), ΔS, param_set)
        T = snow_bulk_temperature(U, FT(1), FT(0.5), parameters.ΔS, param_set)
        @test T ≈ _T_freeze

        U =
            energy_from_T_and_swe.(
                FT(1),
                FT.([272, 274]),
                parameters.ΔS,
                param_set,
            )
        T =
            snow_bulk_temperature.(
                U,
                FT(1),
                FT.([0.0, 1.0]),
                parameters.ΔS,
                param_set,
            )
        @test all(T .≈ FT.([272, 274]))

        #Tests for surface_temp use unphysical args for straightforward formula checks:
        κ_surf_test = FT(pi * _cp_i)
        ρ_surf_test = FT(86400)
        d0 = Snow.diurnal_damping_depth(κ_surf_test, ρ_surf_test, param_set)
        @test d0 ≈ FT(1)
        d1 = Snow.surface_temp_scaling_length(
            κ_surf_test,
            ρ_surf_test,
            FT(log(2)),
            param_set,
        )
        @test d1 ≈ FT(0.5)
        resid_flux_1 = Snow.surface_residual_flux(
            FT(250),
            κ_surf_test,
            ρ_surf_test,
            FT(log(2)),
            param_set,
        )
        @test resid_flux_1 == FT(0)
        resid_flux_2 = Snow.surface_residual_flux(
            _T_freeze + FT(2),
            κ_surf_test,
            ρ_surf_test,
            FT(log(2)),
            param_set,
        )
        @test resid_flux_2 ≈ FT(-4 * κ_surf_test)

        κ_surf_test = FT(0.08)
        ρ_surf_test = FT(500)
        T_sfc_test = FT(272)
        T_bulk_test = FT(270)
        gustiness = SurfaceFluxes.ConstantGustinessSpec(FT(1))
        SW_net = (parameters.α_snow.α - FT(1)) .* FT(20)
        LW_d = FT(20)
        q_l = FT(0)
        h_sfc = FT(0)
        displ = FT(0)
        P_atmos = FT(101325)
        T_atmos = FT(245)
        q_atmos = FT(0.003)
        u_atmos = FT(3)
        atmos_h = FT(1)
        #No snow - should T_bulk
        result = Snow.solve_for_surface_temp_at_a_point(
            T_sfc_test,
            T_bulk_test,
            FT(0),
            κ_surf_test,
            ρ_surf_test,
            ϵ_snow,
            SW_net,
            LW_d,
            q_l,
            h_sfc,
            displ,
            P_atmos,
            T_atmos,
            q_atmos, 
            u_atmos,
            roughness_model,
            atmos_h,
            gustiness, #gustiness
            param_set,
            Snow.EquilibriumGradientTemperatureModel{FT}())
        @test result ==  T_bulk_test
        #nonzero depth
         result = Snow.solve_for_surface_temp_at_a_point(
            T_sfc_test,
            T_bulk_test,
            FT(1),
            κ_surf_test,
            ρ_surf_test,
            ϵ_snow,
            SW_net,
            LW_d,
            q_l,
            h_sfc,
            displ,
            P_atmos,
            T_atmos,
            q_atmos, 
            u_atmos,
            roughness_model,
            atmos_h,
            gustiness, #gustiness
            param_set,
             Snow.EquilibriumGradientTemperatureModel{FT}())
        surface_flux_params = LP.surface_fluxes_parameters(param_set)
        
        T_sfc= result
        q_sfc = Snow.snow_surface_specific_humidity(
            T_sfc,
            q_l,
            T_atmos,
            P_atmos,
            q_atmos,
            atmos_h - h_sfc,
            surface_flux_params,
            thermo_params
        )
        _σ = LP.Stefan(param_set)
        turb_fluxes = ClimaLand.compute_turbulent_fluxes_at_a_point(P_atmos,
                                                                    T_atmos,
                                                                    q_atmos,
                                                                    u_atmos,
                                                                    atmos_h,
                                                                    T_sfc,
                                                                    q_sfc,
                                                                    roughness_model,
                                                                    nothing,
                                                                    nothing,
                                                                    h_sfc,
                                                                    displ,
                                                                    (args...) -> FT(1),
                                                                    (args...) -> FT(1),
                                                                    gustiness,
                                                                    param_set)
        LW_n = -ϵ_snow * (LW_d - _σ * T_sfc^4)
        d = Snow.surface_temp_scaling_length(κ_surf_test, ρ_surf_test, FT(1), param_set)
        sfc_flux = (SW_net + LW_n + turb_fluxes[1] + turb_fluxes[2])
        residual = (sfc_flux + κ_surf_test * (T_sfc - T_bulk_test)/d)
        @test abs(residual)/abs(sfc_flux) < 0.1
    end
    
    @testset "Alternative parameterizations, FT = $FT" begin
        param_set = LP.LandParameters(toml_dict)
        m = Snow.WuWuSnowCoverFractionModel(
            FT(0.08),
            FT(1.77),
            FT(1),
            FT(2);
            z0 = FT(0.05),
        )
        depth = Array(FT(0):FT(0.01):FT(0.2))
        ρ_snow = FT.(range(0.0, 1000.0, 21))
        scf = zeros(FT, length(depth))
        α_snow = zeros(FT, length(depth))
        p = (;
            snow = (;
                z_snow = depth,
                ρ_snow = ρ_snow,
                snow_cover_fraction = scf,
            ),
            drivers = (; cosθs = FT(0.5)),
        )

        update_snow_cover_fraction!(p, m, nothing, 0.0, param_set, nothing)
        @test m.β_scf == max(FT(1.77 - 0.08 * (2 - 1.5)), FT(1.0))
        @test all(
            @. scf ≈
               min(m.β_scf * depth / FT(0.05) / (depth / FT(0.05) + 1), 1)
        )

        m = Snow.ZenithAngleAlbedoModel(FT(0.6), FT(0.1), FT(2); β = FT(0.7))
        update_snow_albedo!(α_snow, m, nothing, p, 0.0, param_set)
        @test all(
            @. α_snow ≈
               min(1 - FT(0.7) * (ρ_snow / FT(1000) - FT(0.2)), 1) *
               (FT(0.6) + (FT(0.1) * exp(-FT(2) * FT(0.5))))
        )

        # test physical limits
        x0 = FT.([0, 1])
        β = FT.([0, 1])
        for x in x0
            for b in β
                m = Snow.ZenithAngleAlbedoModel(
                    FT(0.6),
                    FT(0.0),
                    FT(0);
                    β = b,
                    x0 = x,
                )
                update_snow_albedo!(α_snow, m, nothing, p, 0.0, param_set)
                @test minimum(α_snow) >= FT(0)
                @test maximum(α_snow) <= FT(1)
            end
        end
    end
end
