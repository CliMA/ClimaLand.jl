using Test
import CLIMAParameters as CP
using ClimaLSM.Snow
import ClimaLSM
import ClimaLSM.Parameters as LSMP
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

for FT in (Float32, Float64)
    @testset "Snow Parameterizations, FT = $FT" begin
        param_set = create_lsm_parameters(FT)

        # Density of liquid water (kg/m``^3``)
        _ρ_l = FT(LSMP.ρ_cloud_liq(param_set))
        # Density of ice water (kg/m``^3``)
        _ρ_i = FT(LSMP.ρ_cloud_ice(param_set))
        _cp_l = FT(LSMP.cp_l(param_set))
        _cp_i = FT(LSMP.cp_i(param_set))
        # Reference temperature (K)
        _T_ref = FT(LSMP.T_0(param_set))
        # Freezing temperature (K)
        _T_freeze = FT(LSMP.T_freeze(param_set))
        # Latent heat of fusion at ``T_0`` (J/kg)
        _LH_f0 = FT(LSMP.LH_f0(param_set))
        # Thermal conductivity of dry air
        κ_air = FT(LSMP.K_therm(param_set))

        ρ_snow = FT(200)
        z_0m = FT(0.0024)
        z_0b = FT(0.00024)
        α_snow = FT(0.8)
        ϵ_snow = FT(0.99)
        θ_r = FT(0.08)
        Ksat = FT(1e-3)
        κ_ice = FT(2.21)
        Δt = FT(180.0)
        parameters = SnowParameters{FT}(Δt; earth_param_set = param_set)
        @test parameters.ρ_snow == ρ_snow
        @test typeof(parameters.ρ_snow) == FT
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
        @test runoff_timescale.(z, Ksat, Δt) ≈ max.(Δt, z ./ Ksat)


        U = cat(FT.(Array(LinRange(-1e8, 1e7, 100))), FT(0), dims = 1)
        q_l = snow_liquid_mass_fraction.(U, SWE, Ref(parameters))
        T_bulk = snow_bulk_temperature.(U, SWE, q_l, Ref(parameters))

        @test all(q_l[T_bulk .< _T_freeze] .< eps(FT))
        @test all(q_l[T_bulk .> _T_freeze] .≈ FT(1))
        @test all(q_l[T_bulk .== _T_freeze] .> FT(0.0)) &&
              all(q_l[T_bulk .== _T_freeze] .< FT(1.1))
        Upred =
            _ρ_l .* max.(SWE, eps(FT)) .* (
                specific_heat_capacity.(q_l, Ref(parameters)) .*
                (T_bulk .- _T_ref) .- (1.0f0 .- q_l) .* _LH_f0
            )
        q_lpred = snow_liquid_mass_fraction.(Upred, SWE, Ref(parameters))
        T_pred = snow_bulk_temperature.(Upred, SWE, q_lpred, Ref(parameters))
        @test all(q_lpred .≈ q_l)
        @test all(T_pred .≈ T_bulk)
    end
end
