using Test
import ClimaParams as CP
using ClimaLand.Snow
import ClimaLand
import ClimaLand.Parameters as LP
import Random

Random.seed!(1234)

for FT in (Float32, Float64)
    toml_dict = LP.create_toml_dict(FT)
    @testset "Snow Parameterizations, FT = $FT" begin
        param_set = LP.LandParameters(toml_dict)

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
        κ_ice = FT(2.21)
        ΔS = FT(0.1)
        Δt = Float64(180.0)
        parameters = SnowParameters(
            toml_dict,
            Δt,
            density = densitymodel,
            α_snow = α_snow,
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
        @test parameters.κ_ice == κ_ice
        @test typeof(parameters.κ_ice) == FT


        T = FT.([275.0, 272, _T_freeze])
        @test snow_surface_temperature.(T) ≈ T
        @test specific_heat_capacity(FT(1.0), parameters) == _cp_l
        @test specific_heat_capacity(FT(0.0), parameters) == _cp_i
        ρ_snow = ρ_min
        @test snow_thermal_conductivity(ρ_snow, parameters) ==
              κ_air +
              (FT(0.07) * (ρ_snow / _ρ_i) + FT(0.93) * (ρ_snow / _ρ_i)^2) *
              (κ_ice - κ_air)
        @test all(
            maximum_liquid_mass_fraction.(ρ_min, T, Ref(parameters)) .-
            [FT(0), θ_r * _ρ_l / ρ_snow, θ_r * _ρ_l / ρ_snow] .== 0,
        )

        SWE = cat(FT.(rand(10)), FT(0), dims = 1)
        z = SWE * _ρ_l ./ ρ_snow
        @test runoff_timescale.(z, Ksat, FT(Δt)) ≈ max.(Δt, z ./ Ksat)
        ρ_calc = snow_bulk_density.(SWE, z, parameters)
        @test all(ρ_calc[1:(end - 1)] .≈ ρ_snow)
        @test ρ_calc[end] == _ρ_l
        @test snow_bulk_density(eps(FT(0)), 2 * eps(FT(0)), parameters) == _ρ_l

        U = energy_from_q_l_and_swe(FT(1), FT(0.5), parameters)
        T = snow_bulk_temperature(U, FT(1), FT(0.5), parameters)
        @test T ≈ _T_freeze

        U = energy_from_T_and_swe.(FT(1), FT.([272, 274]), parameters)
        T = snow_bulk_temperature.(U, FT(1), FT.([0.0, 1.0]), parameters)
        @test all(T .≈ FT.([272, 274]))
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
            snow = (; z_snow = depth, ρ_snow = ρ_snow),
            drivers = (; cosθs = FT(0.5)),
        )

        update_snow_cover_fraction!(scf, m, nothing, p, 0.0, param_set)
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
