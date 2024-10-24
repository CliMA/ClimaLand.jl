using Test
import ClimaParams as CP
using ClimaLand.Snow
import ClimaLand
import ClimaLand.Parameters as LP
import Random

Random.seed!(1234)

for FT in (Float32, Float64)
    @testset "Snow Parameterizations, FT = $FT" begin
        param_set = LP.LandParameters(FT)

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

        ρ_snow = FT(200)
        z_0m = FT(0.0024)
        α_snow = FT(0.8)
        # These values should match ClimaParams
        ϵ_snow = FT(0.99)
        z_0b = FT(0.00024)
        θ_r = FT(0.08)
        Ksat = FT(1e-3)
        κ_ice = FT(2.21)
        ρcD_g = FT(1700 * 2.09e3 * 0.1)
        Δt = Float64(180.0)
        parameters = SnowParameters{FT}(Δt; earth_param_set = param_set)
        @test parameters.ρ_snow == ρ_snow
        @test typeof(parameters.ρ_snow) == FT
        @test parameters.ρcD_g == ρcD_g
        @test typeof(parameters.ρcD_g) == FT
        @test parameters.z_0m == z_0m
        @test typeof(parameters.z_0m) == FT
        @test parameters.z_0b == z_0b
        @test typeof(parameters.z_0b) == FT
        @test parameters.α_snow == α_snow
        @test typeof(parameters.α_snow) == FT
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
        @test all(
            maximum_liquid_mass_fraction.(T, ρ_snow, Ref(parameters)) .==
            [FT(0), θ_r * _ρ_l / ρ_snow, θ_r * _ρ_l / ρ_snow],
        )

        SWE = cat(FT.(rand(100) .+ 0.2), FT(0), dims = 1)
        z = snow_depth.(SWE, ρ_snow, _ρ_l)
        @test all(z .== SWE * _ρ_l ./ ρ_snow)
        @test specific_heat_capacity(FT(1.0), parameters) == _cp_l
        @test specific_heat_capacity(FT(0.0), parameters) == _cp_i
        @test snow_thermal_conductivity(ρ_snow, parameters) ==
              κ_air +
              (FT(7.75e-5) * ρ_snow + FT(1.105e-6) * ρ_snow^2) * (κ_ice - κ_air)
        @test runoff_timescale.(z, Ksat, FT(Δt)) ≈ max.(Δt, z ./ Ksat)


        U = cat(FT.(Array(LinRange(-1e8, 1e7, 100))), FT(0), dims = 1)
        q_l = snow_liquid_mass_fraction.(U, SWE, parameters)
        T_bulk = snow_bulk_temperature.(U, SWE, q_l, parameters)

        @test all(q_l[T_bulk .< _T_freeze] .< eps(FT))
        @test all(q_l[T_bulk .> _T_freeze] .≈ FT(1))
        @test all(q_l[T_bulk .== _T_freeze] .> FT(0.0)) &&
              all(q_l[T_bulk .== _T_freeze] .< FT(1.1))
        Upred =
            (
                _ρ_l .* SWE .* specific_heat_capacity.(q_l, parameters) .+
                ρcD_g
            ) .* (T_bulk .- _T_ref) .- _ρ_l .* SWE .* (1.0f0 .- q_l) .* _LH_f0

        q_lpred = snow_liquid_mass_fraction.(Upred, SWE, parameters)
        T_pred = snow_bulk_temperature.(Upred, SWE, q_lpred, parameters)
        @test all(q_lpred .≈ q_l)
        @test all(T_pred .≈ T_bulk)
        @test ClimaLand.Snow.volumetric_internal_energy_liq(FT, parameters) ==
              _ρ_l * _cp_l * (_T_freeze .- _T_ref)
        temp = FT(273.0)
        U = ClimaLand.Snow.energy_from_T_and_swe(FT(1.0), temp, parameters)
        q_l = snow_liquid_mass_fraction(U, FT(1.0), parameters)
        @test T_bulk =
            snow_bulk_temperature(U, FT(1.0), q_l, parameters) == temp

        q_l = FT(0.5)
        U = ClimaLand.Snow.energy_from_q_l_and_swe(FT(1.0), q_l, parameters)
        @test snow_liquid_mass_fraction(U, FT(1.0), parameters) == q_l

    end
end
