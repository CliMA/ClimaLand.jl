using Dates
using Test
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput

using ClimaLand
using ClimaLand.Canopy
using ClimaLand.PlantHydraulics
using ClimaLand.Domains
using ClimaLand.Soil

const FT = Float64

using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations

"""
    create_site_column(FT, lat, long)

Creates and returns a Column domain for user-given latitude and longitude coordinates.
"""
function create_site_column(FT, lat, long, dz_tuple, nelements, zmin, zmax)

    zlim = (zmin, zmax)
    longlat = (long, lat)

    return Domains.Column(; zlim, nelements, longlat, dz_tuple)
end

@testset "US-Ha1 domain info + parameters" begin
    site_ID = FluxnetSimulations.replace_hyphen("US-Ha1")

    # domain information
    (; dz_tuple, nelements, zmin, zmax) =
        FluxnetSimulations.get_domain_info(FT, Val(site_ID))

    @test dz_tuple == (FT(1.5), FT(0.025))
    @test nelements == 20
    @test zmin == FT(-10)
    @test zmax == FT(0)

    # geographical info
    (; time_offset, lat, long, atmos_h) =
        FluxnetSimulations.get_location(FT, Val(site_ID))

    @test time_offset == -5
    @test lat == FT(42.5378)
    @test long == FT(-72.1715)
    @test atmos_h == FT(30)

    # parameters
    (;
        soil_ν,
        soil_K_sat,
        soil_S_s,
        soil_hydrology_cm,
        θ_r,
        ν_ss_quartz,
        ν_ss_om,
        ν_ss_gravel,
        z_0m_soil,
        z_0b_soil,
        soil_ϵ,
        soil_albedo,
        Ω,
        χl,
        G_Function,
        α_PAR_leaf,
        λ_γ_PAR,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        ϵ_canopy,
        ac_canopy,
        g1,
        Drel,
        g0,
        Vcmax25,
        SAI,
        f_root_to_shoot,
        K_sat_plant,
        ψ63,
        Weibull_param,
        a,
        conductivity_model,
        retention_model,
        plant_ν,
        plant_S_s,
        rooting_depth,
        n_stem,
        n_leaf,
        h_leaf,
        h_stem,
        h_canopy,
        z0_m,
        z0_b,
    ) = FluxnetSimulations.get_parameters(FT, Val(site_ID), Vcmax25 = FT(1e-4))

    @test soil_ν == FT(0.5)
    @test soil_K_sat == FT(4e-7)
    @test soil_S_s == FT(1e-3)
    @test θ_r == FT(0.067)
    @test ν_ss_quartz == FT(0.1)
    @test ν_ss_om == FT(0.1)
    @test ν_ss_gravel == FT(0.0)
    @test z_0m_soil == FT(0.01)
    @test z_0b_soil == FT(0.001)
    @test soil_ϵ == FT(0.98)
    @test Ω == FT(0.69)
    @test χl == FT(0.5)
    @test G_Function == ConstantGFunction(χl)
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
    @test Vcmax25 == FT(1e-4)
    @test SAI == FT(1.0)
    @test f_root_to_shoot == FT(3.5)
    @test K_sat_plant == 5e-9
    @test ψ63 == FT(-4 / 0.0098)
    @test Weibull_param == FT(4)
    @test a == FT(0.05 * 0.0098)
    @test conductivity_model ==
          PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
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

@testset "US-MOz domain info + parameters" begin
    site_ID = FluxnetSimulations.replace_hyphen("US-MOz")

    # domain information
    (; dz_tuple, nelements, zmin, zmax) =
        FluxnetSimulations.get_domain_info(FT, Val(site_ID))

    @test dz_tuple == (FT(1.5), FT(0.1))
    @test nelements == 20
    @test zmin == FT(-10)
    @test zmax == FT(0)

    # geographical info
    (; time_offset, lat, long, atmos_h) =
        FluxnetSimulations.get_location(FT, Val(site_ID), time_offset = 5)
    @test time_offset == 5
    @test lat == FT(38.7441)
    @test long == FT(-92.2000)
    @test atmos_h == FT(32)

    # parameters
    (;
        soil_ν,
        soil_K_sat,
        soil_S_s,
        soil_hydrology_cm,
        θ_r,
        ν_ss_quartz,
        ν_ss_om,
        ν_ss_gravel,
        z_0m_soil,
        z_0b_soil,
        soil_ϵ,
        soil_albedo,
        Ω,
        χl,
        G_Function,
        α_PAR_leaf,
        λ_γ_PAR,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        ϵ_canopy,
        ac_canopy,
        g1,
        Drel,
        g0,
        Vcmax25,
        SAI,
        f_root_to_shoot,
        K_sat_plant,
        ψ63,
        Weibull_param,
        a,
        conductivity_model,
        retention_model,
        plant_ν,
        plant_S_s,
        rooting_depth,
        n_stem,
        n_leaf,
        h_leaf,
        h_stem,
        h_canopy,
        z0_m,
        z0_b,
    ) = FluxnetSimulations.get_parameters(FT, Val(site_ID))

    # selected parameters from each "model group" for testing
    @test soil_ν == FT(0.55)
    @test soil_ϵ == FT(0.98)
    @test Ω == FT(0.69)
    @test ac_canopy == FT(5e2)
    @test g0 == FT(1e-4)
    @test Vcmax25 == FT(6e-5)
    @test SAI == FT(1.0)
    @test h_stem == FT(9)
    @test z0_m == FT(0.13) * h_canopy

end

@testset "US-NR1 domain info + parameters" begin
    site_ID = FluxnetSimulations.replace_hyphen("US-NR1")

    # domain information
    (; dz_tuple, nelements, zmin, zmax) =
        FluxnetSimulations.get_domain_info(FT, Val(site_ID))

    @test dz_tuple == (FT(1.25), FT(0.05))
    @test nelements == 20
    @test zmin == FT(-10)
    @test zmax == FT(0)

    (; time_offset, lat, long, atmos_h) =
        FluxnetSimulations.get_location(FT, Val(site_ID))

    @test time_offset == -7
    @test lat == FT(40.0329)
    @test long == FT(-105.5464)
    @test atmos_h == FT(21.5)

    # parameters
    (;
        soil_ν,
        soil_K_sat,
        soil_S_s,
        soil_hydrology_cm,
        θ_r,
        ν_ss_quartz,
        ν_ss_om,
        ν_ss_gravel,
        z_0m_soil,
        z_0b_soil,
        soil_ϵ,
        soil_albedo,
        Ω,
        χl,
        G_Function,
        α_PAR_leaf,
        λ_γ_PAR,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        ϵ_canopy,
        ac_canopy,
        g1,
        Drel,
        g0,
        Vcmax25,
        SAI,
        f_root_to_shoot,
        K_sat_plant,
        ψ63,
        Weibull_param,
        a,
        conductivity_model,
        retention_model,
        plant_ν,
        plant_S_s,
        rooting_depth,
        n_stem,
        n_leaf,
        h_leaf,
        h_stem,
        h_canopy,
        z0_m,
        z0_b,
    ) = FluxnetSimulations.get_parameters(FT, Val(site_ID), Ω = FT(1))

    # selected parameters from each "model group" for testing
    @test soil_ν == FT(0.45)
    @test soil_ϵ == FT(0.98)
    @test Ω == FT(1)
    @test ac_canopy == FT(3e3)
    @test g0 == FT(1e-4)
    @test Vcmax25 == FT(9e-5)
    @test SAI == FT(1.0)
    @test h_stem == FT(7.5)
    @test z0_m == FT(0.13) * h_canopy
end


@testset "US-Var domain info + parameters" begin
    site_ID = FluxnetSimulations.replace_hyphen("US-Var")

    # domain information
    (; dz_tuple, nelements, zmin, zmax) =
        FluxnetSimulations.get_domain_info(FT, Val(site_ID))

    @test dz_tuple == FT.((0.05, 0.02))
    @test nelements == 24
    @test zmin == FT(-0.5)
    @test zmax == FT(0)


    (; time_offset, lat, long, atmos_h) =
        FluxnetSimulations.get_location(FT, Val(site_ID))
    @test time_offset == -8
    @test lat == FT(38.4133)
    @test long == FT(-120.9508)
    @test atmos_h == FT(2)

    # parameters
    (;
        soil_ν,
        soil_K_sat,
        soil_S_s,
        soil_hydrology_cm,
        θ_r,
        ν_ss_quartz,
        ν_ss_om,
        ν_ss_gravel,
        z_0m_soil,
        z_0b_soil,
        soil_ϵ,
        soil_albedo,
        Ω,
        χl,
        G_Function,
        α_PAR_leaf,
        λ_γ_PAR,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        ϵ_canopy,
        ac_canopy,
        g1,
        Drel,
        g0,
        Vcmax25,
        SAI,
        f_root_to_shoot,
        K_sat_plant,
        ψ63,
        Weibull_param,
        a,
        conductivity_model,
        retention_model,
        plant_ν,
        plant_S_s,
        rooting_depth,
        n_stem,
        n_leaf,
        h_leaf,
        h_stem,
        h_canopy,
        z0_m,
        z0_b,
    ) = FluxnetSimulations.get_parameters(FT, Val(site_ID), g0 = FT(5e-4))

    # selected parameters from each "model group" for testing
    @test soil_ν == FT(0.5)
    @test soil_ϵ == FT(0.98)
    @test Ω == FT(0.75)
    @test ac_canopy == FT(745)
    @test g0 == FT(5e-4)
    @test Vcmax25 == FT(2.5e-5)
    @test SAI == FT(0)
    @test h_stem == FT(0)
    @test z0_m == FT(0.13) * h_canopy
end

@testset "generic site domain info + parameters" begin
    site_ID = "AU-Emr"

    (; time_offset, lat, long, atmos_h) =
        FluxnetSimulations.get_location(site_ID)

    @test lat == FT(-23.8587)
    @test long == FT(148.4746)
    @test time_offset == -10
    @test atmos_h == FT(6.7)

    # domain information
    (; dz_tuple, nelements, zmin, zmax) = FluxnetSimulations.get_domain_info(FT)

    @test dz_tuple == (FT(1.5), FT(0.1))
    @test nelements == 20
    @test zmin == FT(-10)
    @test zmax == FT(0)

    domain = create_site_column(FT, lat, long, dz_tuple, nelements, zmin, zmax)

    (;
        soil_ν,
        soil_K_sat,
        soil_S_s,
        soil_hydrology_cm,
        θ_r,
        ν_ss_quartz,
        ν_ss_om,
        ν_ss_gravel,
        z_0m_soil,
        z_0b_soil,
        soil_ϵ,
        soil_albedo,
        Ω,
        χl,
        G_Function,
        α_PAR_leaf,
        λ_γ_PAR,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        ϵ_canopy,
        ac_canopy,
        g1,
        Drel,
        g0,
        Vcmax25,
        SAI,
        f_root_to_shoot,
        K_sat_plant,
        ψ63,
        Weibull_param,
        a,
        conductivity_model,
        retention_model,
        plant_ν,
        plant_S_s,
        rooting_depth,
        n_stem,
        n_leaf,
        h_leaf,
        h_stem,
        h_canopy,
        z0_m,
        z0_b,
    ) = FluxnetSimulations.get_parameters(FT, site_ID, domain; θ_r = FT(0.067))

    @test θ_r == FT(0.067)
    @test Ω == 0.75
    @test g0 == FT(100.0)
    @test χl == FT(-0.3)
    @test h_leaf == FT(0.5)
    @test Vcmax25 == FT(2.4e-5)
end
nothing
