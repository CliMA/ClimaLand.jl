using Dates
using Test
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput

using ClimaLand
using ClimaLand.Canopy

const FT = Float64

using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations

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
    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(site_ID))

    @test time_offset == -5
    @test lat == FT(42.5378)
    @test long == FT(-72.1715)

    (; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID))
    @test atmos_h == FT(30)

    # parameters
    (;
        soil_ν,
        soil_K_sat,
        soil_S_s,
        soil_vg_n,
        soil_vg_α,
        θ_r,
        ν_ss_quartz,
        ν_ss_om,
        ν_ss_gravel,
        z_0m_soil,
        z_0b_soil,
        soil_ϵ,
        soil_α_PAR,
        soil_α_NIR,
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
        h_canopy,
    ) = FluxnetSimulations.get_parameters(FT, Val(site_ID), Vcmax25 = FT(1e-4))

    @test soil_ν == FT(0.5)
    @test soil_K_sat == FT(4e-7)
    @test soil_S_s == FT(1e-3)
    @test soil_vg_n == FT(2.05)
    @test soil_vg_α == FT(4.0)
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
    @test K_sat_plant == 7e-8
    @test ψ63 == FT(-4 / 0.0098)
    @test Weibull_param == FT(4)
    @test a == FT(0.05 * 0.0098)
    @test conductivity_model ==
          Canopy.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
    @test retention_model == Canopy.LinearRetentionCurve{FT}(a)
    @test plant_ν == FT(2.46e-4)
    @test plant_S_s == FT(1e-2 * 0.0098)
    @test rooting_depth == FT(0.5)
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
    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(site_ID))
    @test time_offset == -6
    @test lat == FT(38.7441)
    @test long == FT(-92.2000)

    (; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID))
    @test atmos_h == FT(32)

    # parameters
    (;
        soil_ν,
        soil_K_sat,
        soil_S_s,
        soil_vg_n,
        soil_vg_α,
        θ_r,
        ν_ss_quartz,
        ν_ss_om,
        ν_ss_gravel,
        z_0m_soil,
        z_0b_soil,
        soil_ϵ,
        soil_α_PAR,
        soil_α_NIR,
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
        h_canopy,
    ) = FluxnetSimulations.get_parameters(FT, Val(site_ID))

    # selected parameters from each "model group" for testing
    @test soil_ν == FT(0.55)
    @test soil_ϵ == FT(0.98)
    @test Ω == FT(0.69)
    @test ac_canopy == FT(5e2)
    @test g0 == FT(1e-4)
    @test Vcmax25 == FT(6e-5)
    @test SAI == FT(1.0)

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

    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(site_ID))

    @test time_offset == -7
    @test lat == FT(40.0329)
    @test long == FT(-105.5464)

    (; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID))
    @test atmos_h == FT(21.5)

    # parameters
    (;
        soil_ν,
        soil_K_sat,
        soil_S_s,
        soil_vg_n,
        soil_vg_α,
        θ_r,
        ν_ss_quartz,
        ν_ss_om,
        ν_ss_gravel,
        z_0m_soil,
        z_0b_soil,
        soil_ϵ,
        soil_α_PAR,
        soil_α_NIR,
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
        h_canopy,
    ) = FluxnetSimulations.get_parameters(FT, Val(site_ID), Ω = FT(1))

    # selected parameters from each "model group" for testing
    @test soil_ν == FT(0.45)
    @test soil_ϵ == FT(0.98)
    @test Ω == FT(1)
    @test ac_canopy == FT(3e3)
    @test g0 == FT(1e-4)
    @test Vcmax25 == FT(9e-5)
    @test SAI == FT(1.0)
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

    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(site_ID))
    @test time_offset == -8
    @test lat == FT(38.4133)
    @test long == FT(-120.9508)

    (; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID))
    @test atmos_h == FT(2)

    # parameters
    (;
        soil_ν,
        soil_K_sat,
        soil_S_s,
        soil_vg_n,
        soil_vg_α,
        θ_r,
        ν_ss_quartz,
        ν_ss_om,
        ν_ss_gravel,
        z_0m_soil,
        z_0b_soil,
        soil_ϵ,
        soil_α_PAR,
        soil_α_NIR,
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
        h_canopy,
    ) = FluxnetSimulations.get_parameters(FT, Val(site_ID), g0 = FT(5e-4))

    # selected parameters from each "model group" for testing
    @test soil_ν == FT(0.5)
    @test soil_ϵ == FT(0.98)
    @test Ω == FT(0.75)
    @test ac_canopy == FT(745)
    @test g0 == FT(5e-4)
    @test Vcmax25 == FT(2.5e-5)
    @test SAI == FT(0)
end

@testset "Sensor depths" begin
    depths = FluxnetSimulations.get_sensor_depths(FT, Val(:US_Var))
    @test depths.tsoil == FT.((0.02, 0.04, 0.08, 0.16, 0.32))
    @test depths.swc == FT.((0, 0.10, 0.20))
    # Sites without documented depths fall back to unknown
    @test FluxnetSimulations.get_sensor_depths(FT, Val(:US_MOz)) ==
          (; tsoil = nothing, swc = nothing)
end

@testset "FLUXNET2015 metadata" begin
    metadata_path = joinpath(mktempdir(), "metadata_DD_clean.csv")
    write(
        metadata_path,
        """
        site_id,latitude,longitude,utc_offset,annual_temp,annual_precip,canopy_height,atmospheric_sensor_heights,swc_depths,ts_depths
        BE-Vie,50.3049,5.9981,1.0,7.8,1062.0,30.0,40.0;52.0,NaN,NaN
        AU-ASM,-22.283,133.249,9.5,,,6.5,11.6,NaN,NaN
        XX-Nah,10.0,20.0,-3.0,,,,,NaN,NaN
        XX-Nol,,20.0,-3.0,,,,11.0,NaN,NaN
        """,
    )
    kw = (; fluxnet2015_metadata_path = metadata_path)

    info = FluxnetSimulations.get_site_info("BE-Vie"; kw...)
    @test info.time_offset === 1
    @test info.atmospheric_sensor_height == [40.0, 52.0]

    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(:BE_Vie); kw...)
    @test (time_offset, lat, long) == (1, FT(50.3049), FT(5.9981))
    @test FluxnetSimulations.get_fluxtower_height(FT, Val(:BE_Vie); kw...) ==
          (; atmos_h = FT(52))
    @test FluxnetSimulations.get_location(FT, Val(:AU_ASM); kw...).time_offset ==
          9.5
    @test FluxnetSimulations.get_fluxtower_height(FT, Val(:AU_ASM); kw...) ==
          (; atmos_h = FT(11.6))
    @test FluxnetSimulations.get_canopy_height("BE-Vie"; kw...) == 30.0
    @test_throws ErrorException FluxnetSimulations.get_fluxtower_height(
        FT,
        Val(:XX_Nah);
        kw...,
    )
    @test_throws ErrorException FluxnetSimulations.get_location(
        FT,
        Val(:XX_Nol);
        kw...,
    )
    @test_throws ErrorException FluxnetSimulations.get_site_info(
        "XX-Abc";
        kw...,
    )

    # A site with a hardcoded configuration does not read the metadata
    @test FluxnetSimulations.get_location(FT, Val(:US_MOz)).time_offset == -6
end

@testset "get_data_dates with required_columns" begin
    (start_date, stop_date) = FluxnetSimulations.get_data_dates(
        "US-MOz",
        -6;
        duration = Day(1),
        required_columns = FluxnetSimulations.FLUXNET_FORCING_COLUMNS,
    )
    @test stop_date - start_date == Day(1)
    @test start_date >= first(FluxnetSimulations.get_data_dates("US-MOz", -6))
end
