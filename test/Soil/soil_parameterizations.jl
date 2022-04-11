using Test
using UnPack
using CLIMAParameters.Planet: ρ_cloud_liq, ρ_cloud_ice, cp_l, cp_i, T_0, LH_f0
using CLIMAParameters.Atmos.Microphysics: K_therm
using CLIMAParameters: AbstractEarthParameterSet
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM.Soil

@testset "Integrated Energy and Hydrology Parameterizations" begin
    FT = Float64
    struct EarthParameterSet <: AbstractEarthParameterSet end
    param_set = EarthParameterSet()

    # Density of liquid water (kg/m``^3``)
    _ρ_l = FT(ρ_cloud_liq(param_set))
    # Density of ice water (kg/m``^3``)
    _ρ_i = FT(ρ_cloud_ice(param_set))
    # Volum. isobaric heat capacity liquid water (J/m3/K)
    _ρcp_l = FT(cp_l(param_set) * _ρ_l)
    # Volumetric isobaric heat capacity ice (J/m3/K)
    _ρcp_i = FT(cp_i(param_set) * _ρ_i)
    # Reference temperature (K)
    _T_ref = FT(T_0(param_set))
    # Latent heat of fusion at ``T_0`` (J/kg)
    _LH_f0 = FT(LH_f0(param_set))
    # Thermal conductivity of dry air
    κ_air = FT(K_therm(param_set))

    ν = 0.2
    S_s = 1e-3
    θ_r = 0.1
    vg_α = 2.0
    vg_n = 1.4
    K_sat = 1e-5
    ν_ss_om = 0.1
    ν_ss_gravel = 0.1
    ν_ss_quartz = 0.1
    ρc_ds = 1e6
    κ_solid = 0.1
    ρp = 1.0
    κ_sat_unfrozen = 0.0
    κ_sat_frozen = 0.0
    κ_dry = Soil.κ_dry(ρp, ν, κ_solid, κ_air)
    parameters = EnergyHydrologyParameters(;
        κ_dry = κ_dry,
        κ_sat_frozen = κ_sat_frozen,
        κ_sat_unfrozen = κ_sat_unfrozen,
        ρc_ds = ρc_ds,
        ν = ν,
        ν_ss_om = ν_ss_om,
        ν_ss_quartz = ν_ss_quartz,
        ν_ss_gravel = ν_ss_gravel,
        vg_α = vg_α,
        vg_n = vg_n,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
        earth_param_set = param_set,
    )
    # Test that the preset params are set properly
    @unpack α, β, Ω, γ, vg_m, γT_ref = parameters
    @test α == 0.24
    @test β == 18.3
    @test γ == 2.64e-2
    @test Ω == 7.0
    @test γT_ref == 288.0
    @test vg_m == 1.0 - 1.0 / vg_n
    @test typeof.([α, β, γ, Ω, vg_m, γT_ref]) == [FT, FT, FT, FT, FT, FT]

    @test temperature_from_ρe_int(5.4e7, 0.05, 2.1415e6, parameters) ==
          FT(_T_ref + (5.4e7 + 0.05 * _ρ_i * _LH_f0) / 2.1415e6)

    @test volumetric_heat_capacity(0.25, 0.05, parameters) ==
          FT(1e6 + 0.25 * _ρcp_l + 0.05 * _ρcp_i)

    @test volumetric_internal_energy(0.05, 2.1415e6, 300.0, parameters) ==
          FT(2.1415e6 * (300.0 - _T_ref) - 0.05 * _ρ_i * _LH_f0)

    @test κ_sat(0.25, 0.05, 0.57, 2.29) ==
          FT(0.57^(0.25 / (0.05 + 0.25)) * 2.29^(0.05 / (0.05 + 0.25)))

    @test κ_sat(0.0, 0.0, 0.57, 2.29) == FT(0.5 * (0.57 + 2.29))

    @test relative_saturation(0.25, 0.05, 0.4) == FT((0.25 + 0.05) / 0.4)

    # ice fraction = 0
    @test kersten_number(0.0, 0.75, parameters) == FT(
        0.75^((FT(1) + 0.1 - 0.24 * 0.1 - 0.1) / FT(2)) *
        (
            (FT(1) + exp(-18.3 * 0.75))^(-FT(3)) -
            ((FT(1) - 0.75) / FT(2))^FT(3)
        )^(FT(1) - 0.1),
    )

    # ice fraction ~= 0
    @test kersten_number(0.05, 0.75, parameters) == FT(0.75^(FT(1) + 0.1))

    @test thermal_conductivity(1.5, 0.7287, 0.7187) ==
          FT(0.7287 * 0.7187 + (FT(1) - 0.7287) * 1.5)

    @test volumetric_internal_energy_liq(300.0, parameters) ==
          FT(_ρcp_l * (300.0 - _T_ref))

    @test Soil.κ_solid(FT(0.5), FT(0.25), FT(2.0), FT(3.0), FT(2.0)) ≈
          FT(2)^FT(0.5) * FT(2)^FT(0.25) * FT(3.0)^FT(0.25)

    @test Soil.κ_sat_frozen(FT(0.5), FT(0.1), FT(0.4)) ==
          FT(0.5)^FT(0.9) * FT(0.4)^FT(0.1)

    @test Soil.κ_sat_unfrozen(FT(0.5), FT(0.1), FT(0.4)) ==
          FT(0.5)^FT(0.9) * FT(0.4)^FT(0.1)

    @test Soil.κ_dry(ρp, ν, κ_solid, κ_air) ==
          ((FT(0.053) * FT(0.1) - κ_air) * FT(0.8) + κ_air * FT(1.0)) /
          (FT(1.0) - (FT(1.0) - FT(0.053)) * FT(0.8))


    # Impedance factor
    @test impedance_factor(FT(1.0), parameters.Ω) ≈ 1e-7

    # Viscosity Factor
    T = FT.([278.0, 288.0, 298.0])
    @test viscosity_factor.(T, parameters.γ, parameters.γT_ref) ≈
          exp.(parameters.γ .* (T .- parameters.γT_ref))

end

@testset "Richards equation and van Genuchten closures" begin
    FT = Float32
    θ_r = FT(0.2)
    ν = FT(0.4)
    S_s = FT(1e-2)
    vg_α = FT(3.6)
    vg_n = FT(1.56)
    K_sat = FT(2.9e-7)
    vg_m = FT(1.0 - 1.0 / vg_n)
    richards_parameters = RichardsParameters(;
        ν = ν,
        vg_α = vg_α,
        vg_n = vg_n,
        vg_m = vg_m,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
    )
    # Effective _saturation
    θ = FT.([0.3, 0.4, 0.5])
    S = effective_saturation.(ν, θ, θ_r)
    @test S ≈ [0.5, 1.0, 1.5]
    @test eltype(S) == FT

    # Matric Potential and inverse
    va = -((S[1]^(-FT(1) / vg_m) - FT(1)) * vg_α^(-vg_n))^(FT(1) / vg_n)
    ψ = matric_potential.(vg_α, vg_n, vg_m, S[1:2])
    @test ψ ≈ [va, 0.0]
    @test eltype(ψ) == FT

    #Pressure head
    p = pressure_head.(vg_α, vg_n, vg_m, θ_r, θ, ν, S_s)
    @test p ≈ push!(ψ, FT(10.0))
    @test eltype(p) == FT

    # Hydraulic K
    k = hydraulic_conductivity.(K_sat, vg_m, S)
    va =
        (sqrt(S[1]) * (FT(1) - (FT(1) - S[1]^(FT(1) / vg_m))^vg_m)^FT(2)) *
        K_sat
    @test k ≈ [va, K_sat, K_sat]
    @test eltype(k) == FT

    # Volumetric Liquid Fraction
    vlf = volumetric_liquid_fraction.(FT.([0.25, 0.5, 0.75]), FT(0.5))
    @test vlf ≈ FT.([0.25, 0.5, 0.5])
end
