using Test
import CLIMAParameters as CP
using ClimaLSM.Soil
import ClimaLSM
import ClimaLSM.Parameters as LSMP
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

@testset "integrated Energy and Hydrology Parameterizations" begin
    FT = Float64
    param_set = create_lsm_parameters(FT)

    # Density of liquid water (kg/m``^3``)
    _ρ_l = FT(LSMP.ρ_cloud_liq(param_set))
    # Density of ice water (kg/m``^3``)
    _ρ_i = FT(LSMP.ρ_cloud_ice(param_set))
    # Volum. isobaric heat capacity liquid water (J/m3/K)
    _ρcp_l = FT(LSMP.cp_l(param_set) * _ρ_l)
    # Volumetric isobaric heat capacity ice (J/m3/K)
    _ρcp_i = FT(LSMP.cp_i(param_set) * _ρ_i)
    # Reference temperature (K)
    _T_ref = FT(LSMP.T_0(param_set))
    # Latent heat of fusion at ``T_0`` (J/kg)
    _LH_f0 = FT(LSMP.LH_f0(param_set))
    # Thermal conductivity of dry air
    κ_air = FT(LSMP.K_therm(param_set))

    ν = 0.2
    S_s = 1e-3
    θ_r = 0.1
    vg_α = 2.0
    vg_n = 1.4
    hcm = vanGenuchten(; α = vg_α, n = vg_n)
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
    parameters = EnergyHydrologyParameters{FT}(;
        κ_dry = κ_dry,
        κ_sat_frozen = κ_sat_frozen,
        κ_sat_unfrozen = κ_sat_unfrozen,
        ρc_ds = ρc_ds,
        ν = ν,
        ν_ss_om = ν_ss_om,
        ν_ss_quartz = ν_ss_quartz,
        ν_ss_gravel = ν_ss_gravel,
        hydrology_cm = hcm,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
        earth_param_set = param_set,
    )
    # Test that the preset params are set properly
    (; α, β, Ω, hydrology_cm, γ, γT_ref) = parameters
    @test α == 0.24
    @test β == 18.3
    @test γ == 2.64e-2
    @test Ω == 7.0
    @test γT_ref == 288.0
    # Test derived parameters
    @test hydrology_cm.m == 1.0 - 1.0 / vg_n
    @test hydrology_cm.S_c ==
          (1 + ((vg_n - 1) / vg_n)^(1 - 2 * vg_n))^(-hydrology_cm.m)
    @test typeof.([α, β, γ, Ω, γT_ref]) == [FT, FT, FT, FT, FT]

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

    x = [0.0, 0.2, 0.4]
    @test soil_tortuosity.(x, 0.0, 0.3) ≈
          [0.3^1.5, 0.1^2.5 / 0.3, eps(Float64)^2.5 / 0.3]
    x = [0.35, 0.25, 0.0, 0.1]
    y = [0.0, 0.0, 0.1, 0.1 * 0.15 / 0.25]
    @test dry_soil_layer_thickness.(x, 0.25, 0.1) ≈ y
    @test soil_resistance(θ_r, θ_r - eps(FT), FT(0.0), parameters) ≈
          dry_soil_layer_thickness(
        effective_saturation(ν, θ_r - eps(FT), θ_r),
        hcm.S_c,
        parameters.d_ds,
    ) / FT(LSMP.D_vapor(param_set)) / soil_tortuosity(θ_r + eps(FT), 0.0, ν)
    @test soil_resistance(ν, ν + eps(FT), FT(0.0), parameters) ≈
          dry_soil_layer_thickness(
        effective_saturation(ν, ν + eps(FT), θ_r),
        hcm.S_c,
        parameters.d_ds,
    ) / FT(LSMP.D_vapor(param_set)) / soil_tortuosity(ν, 0.0, ν)
end

@testset "Brooks and Corey closure" begin
    FT = Float32
    hcm = BrooksCorey(; ψb = FT(-0.09), c = FT(0.228))
    # Test derived parameter
    @test hcm.S_c == (1 + 1 / FT(0.228))^(-FT(0.228))

    S = FT.([0.5, 1.0, 1.5])
    θ = FT.([0.3, 0.4, 0.5])
    θ_r = FT(0.2)
    ν = FT(0.4)
    K_sat = FT(2.9e-7)
    S_s = FT(1e-2)

    # Matric Potential and inverse
    va = S[1]^(-1 / hcm.c) * hcm.ψb
    ψ = matric_potential.(hcm, S[1:2])
    @test inverse_matric_potential.(hcm, ψ) ≈ S[1:2]
    @test ψ ≈ [va, hcm.ψb]
    @test eltype(ψ) == FT

    #Pressure head
    p = pressure_head.(hcm, θ_r, θ, ν, S_s)
    @test p ≈ push!(ψ, FT(0.1 / 1e-2 + hcm.ψb))
    @test eltype(p) == FT
    # Derivative
    dψ = dψdϑ.(hcm, θ, ν, θ_r, S_s)
    va = @. -hcm.ψb / (hcm.c * (ν - θ_r)) * S^(-(1 + 1 / hcm.c))
    @test dψ ≈ [va[1], 1 / S_s, 1 / S_s]

    # Hydraulic K
    k = @. hydraulic_conductivity.(hcm, K_sat, S)
    va = S[1]^(2 / hcm.c + 3) * K_sat
    @test k ≈ [va, K_sat, K_sat]
    @test eltype(k) == FT


end

@testset "Richards equation and van Genuchten closures" begin
    FT = Float32
    θ_r = FT(0.2)
    ν = FT(0.4)
    S_s = FT(1e-2)
    vg_α = FT(3.6)
    vg_n = FT(1.56)
    hcm = vanGenuchten(; α = vg_α, n = vg_n)
    K_sat = FT(2.9e-7)
    vg_m = FT(1.0 - 1.0 / vg_n)
    richards_parameters = RichardsParameters(;
        ν = ν,
        hydrology_cm = hcm,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
    )
    # Effective _saturation
    θ = FT.([0.3, 0.4, 0.5])
    S = Soil.effective_saturation.(ν, θ, θ_r)
    @test S ≈ [0.5, 1.0, 1.5]
    @test eltype(S) == FT

    # Matric Potential and inverse
    va = -((S[1]^(-FT(1) / vg_m) - FT(1)) * vg_α^(-vg_n))^(FT(1) / vg_n)
    ψ = matric_potential.(hcm, S[1:2])
    @test inverse_matric_potential.(hcm, ψ) ≈ S[1:2]
    @test ψ ≈ [va, 0.0]
    @test eltype(ψ) == FT
    # Derivative
    dψ = dψdϑ.(hcm, θ, ν, θ_r, S_s)
    va = FT(
        1.0 / (vg_α * vg_m * vg_n) / (ν - θ_r) *
        (S[1]^(-1 / vg_m) - 1)^(1 / vg_n - 1) *
        S[1]^(-1 / vg_m - 1),
    )
    @test dψ ≈ [va, 1 / S_s, 1 / S_s]
    #Pressure head
    p = pressure_head.(hcm, θ_r, θ, ν, S_s)
    @test p ≈ push!(ψ, FT(10.0))
    @test eltype(p) == FT

    # Hydraulic K
    k = hydraulic_conductivity.(hcm, K_sat, S)
    va =
        (sqrt(S[1]) * (FT(1) - (FT(1) - S[1]^(FT(1) / vg_m))^vg_m)^FT(2)) *
        K_sat
    @test k ≈ [va, K_sat, K_sat]
    @test eltype(k) == FT

    # Volumetric Liquid Fraction
    vlf = volumetric_liquid_fraction.(FT.([0.25, 0.5, 0.75]), FT(0.5), FT(0.0))
    @test vlf ≈ FT.([0.25, 0.5, 0.5])
end

@testset "Freezing and Thawing" begin
    FT = Float64
    param_set = create_lsm_parameters(FT)

    # Density of liquid water (kg/m``^3``)
    _ρ_l = FT(LSMP.ρ_cloud_liq(param_set))
    # Density of ice water (kg/m``^3``)
    _ρ_i = FT(LSMP.ρ_cloud_ice(param_set))
    # Volum. isobaric heat capacity liquid water (J/m3/K)
    _ρcp_l = FT(LSMP.cp_l(param_set) * _ρ_l)
    # Volumetric isobaric heat capacity ice (J/m3/K)
    _ρcp_i = FT(LSMP.cp_i(param_set) * _ρ_i)
    # Reference temperature (K)
    _T_ref = FT(LSMP.T_0(param_set))
    # Latent heat of fusion at ``T_0`` (J/kg)
    _LH_f0 = FT(LSMP.LH_f0(param_set))
    # Gravitational acceleration
    _grav = FT(LSMP.grav(param_set))
    # T freeze
    _T_freeze = FT(LSMP.T_freeze(param_set))

    # Thermal conductivity of dry air
    κ_air = FT(LSMP.K_therm(param_set))

    ν = 0.2
    S_s = 1e-3
    θ_r = 0.1
    vg_α = 2.0
    vg_n = 1.4
    hcm = vanGenuchten(; α = vg_α, n = vg_n)

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
    parameters = EnergyHydrologyParameters{FT}(;
        κ_dry = κ_dry,
        κ_sat_frozen = κ_sat_frozen,
        κ_sat_unfrozen = κ_sat_unfrozen,
        ρc_ds = ρc_ds,
        ν = ν,
        ν_ss_om = ν_ss_om,
        ν_ss_quartz = ν_ss_quartz,
        ν_ss_gravel = ν_ss_gravel,
        hydrology_cm = hcm,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
        earth_param_set = param_set,
    )

    Δz = FT(1.0)
    @test thermal_time(ρc_ds, Δz, κ_dry) == ρc_ds * Δz^2 / κ_dry
    θ_l = FT.([0.11, 0.15, ν])
    θ_i = FT(0.0)
    T = FT(273)
    θtot = @.(_ρ_i / _ρ_l * θ_i + θ_l)
    ψ0 = @. matric_potential(hcm, Soil.effective_saturation(ν, θtot, θ_r))
    ψT = @.(_LH_f0 / _T_freeze / _grav * (T - _T_freeze))
    θ_star = @. inverse_matric_potential(hcm, ψ0 + ψT) * (ν - θ_r) + θ_r
    @test (
        phase_change_source.(θ_l, θ_i, T, FT(1.0), Ref(parameters)) ≈
        θ_l .- θ_star
    )
    @test sum(
        phase_change_source.(θ_l, θ_i, T, FT(1.0), Ref(parameters)) .> 0.0,
    ) == 3
    # try θ_l = 0.1

    θ_l = FT.([0.11, 0.15, ν])
    θ_i = FT(0.0)
    T = FT(274)
    θtot = @.(_ρ_i / _ρ_l * θ_i + θ_l)
    ψ0 = @. matric_potential(hcm, Soil.effective_saturation(ν, θtot, θ_r))
    ψT = FT(0.0)
    θ_star = @. inverse_matric_potential(hcm, ψ0 + ψT) * (ν - θ_r) + θ_r
    @test (
        phase_change_source.(θ_l, θ_i, T, FT(1.0), Ref(parameters)) ≈
        zeros(FT, 3)
    )
    @test (θ_star ≈ θ_l)



    θ_l = FT(0.01)
    θ_i = FT.([0.05, 0.08])
    T = FT(274)
    θtot = @.(_ρ_i / _ρ_l * θ_i + θ_l)
    ψ0 = @. matric_potential(hcm, Soil.effective_saturation(ν, θtot, θ_r))
    ψT = FT(0.0)
    θ_star = @. inverse_matric_potential(hcm, ψ0 + ψT) * (ν - θ_r) + θ_r
    @test (
        phase_change_source.(θ_l, θ_i, T, FT(1.0), Ref(parameters)) ≈
        θ_l .- θ_star
    )
    @test sum(
        phase_change_source.(θ_l, θ_i, T, FT(1.0), Ref(parameters)) .< 0.0,
    ) == 2


    θ_l = FT(0.11)
    θ_i = FT.([0.05, 0.08])
    T = FT(260)
    θtot = @.(_ρ_i / _ρ_l * θ_i + θ_l)
    ψ0 = @. matric_potential(hcm, Soil.effective_saturation(ν, θtot, θ_r))
    ψT = @.(_LH_f0 / _T_freeze / _grav * (T - _T_freeze))
    θ_star = @. inverse_matric_potential(hcm, ψ0 + ψT) * (ν - θ_r) + θ_r
    @test (
        phase_change_source.(θ_l, θ_i, T, FT(1.0), Ref(parameters)) ≈
        θ_l .- θ_star
    )
    @test sum(
        phase_change_source.(θ_l, θ_i, T, FT(1.0), Ref(parameters)) .> 0.0,
    ) == 2

end
