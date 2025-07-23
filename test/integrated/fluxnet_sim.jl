using Dates
using Test
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput

using ClimaLand
using FluxnetSimulations
const FT = Float64


@testset "US-MOz domain info + parameters" begin
    site_ID = FluxnetSimulations.replace_hyphen("US-MOz")
    (; dz_tuple,
        nelements,
        zmin
        zmax) = FluxnetSimulations.get_domain_info(Val(site_ID))

    @test dz_tuple == (FT(1.5), FT(0.025))
    @test nelements == 20
    @test zmin == FT(-10)
    @test xmax == FT(0)

    (; time_offset, lat, long, atmos_h, soil_ν, soil_K_sat,
        soil_S_s, soil_vg_n, soil_vg_α, θ_r, ν_ss_quartz, ν_ss_om,
        ν_ss_gravel, z_0m_soil, z_0b_soil, soil_ϵ, soil_α_PAR, soil_α_NIR,
        Ω, ld, G_Function, α_PAR_leaf, λ_γ_PAR, τ_PAR_leaf,
        α_NIR_leaf, τ_NIR_leaf, ϵ_canopy, ac_canopy, g1, Drel,
        g0, Vcmax25, SAI, f_root_to_shoot, K_sat_plant, ψ63,
        Weibull_param, a, conductivity_model, retention_model,
        plant_ν, plant_S_s, rooting_depth, n_stem, n_leaf,
        h_leaf, h_stem, h_canopy, z0_m, z0_b) = FluxnetSimulations.get_parameters(site_ID)

    @test time_offset == 5
    @test lat == FT(42.5378)
    @test long == FT(-72.1715)
    @test atmos_h == FT(30)
    @test soil_ν == FT(0.5)
    @test soil_K_sat == FT(4e-7)
    @test soil_S_s == FT(1e-3)
    @test soil_vg_n == FT(2.05)
    @test soil_vg_α == FT(0.04)
    @test θ_r == FT(0.067)
    @test ν_ss_quartz == FT(0.1)
    @test ν_ss_om == FT(0.1)
    @test ν_ss_gravel == FT(0.0)
    @test z_0m_soil == FT(0.01)
    @test z_0b_soil == FT(0.001)
    @test soil_ϵ == FT(0.98)
    @test soil_α_PAR == FT(0.2)
    @test soil_α_NIR == FT(0.2)
    @test Ω == FT(0.69)
    @test ld == FT(0.5)
    @test G_Function == ConstantGFunction(ld)
    @test α_PAR_leaf == FT(0.1)
    @test λ_γ_PAR == FT(5e-7)
    @test τ_PAR_leaf == FT(0.05)
    @test α_NIR_leaf == FT(0.45)
    @test τ_NIR_leaf == FT(0.25)
    @test ϵ_canopy == FT(0.97)
    @test ac_canopy == FT(2.5e3)
    @test g1 == FT(141)
    @test Drel == FT(1.6)
    @test g0 == FT(1e-4)
    @test Vcmax25 == FT(9e-5)
    @test SAI == FT(1.0)
    @test f_root_to_shoot == FT(3.5)
    @test K_sat_plant == 5e-9
    @test ψ63 == FT(-4 / 0.0098)
    @test Weibull_param == FT(4)
    @test a == FT(0.05 * 0.0098)
    @test conductivity_model == PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
    @test retention_model == PlantHydraulics.LinearRetentionCurve{FT}(a)
    @test plant_ν == FT(2.46e-4)
    @test plant_S_s == FT(1e-2 * 0.0098)
    @test rooting_depth == FT(0.5)
    @test n_stem == Int64(1)
    @test n_leaf == Int64(1)
    @test h_leaf == FT(12)
    @test h_stem == FT(14)
    @test z0_m == FT(0.13) * h_canopy
    @test z0_b == FT(0.1) * z0_m
end