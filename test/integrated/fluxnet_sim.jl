using Dates
using Test
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput

using ClimaLand
using ClimaLand.Canopy
import ClimaLand.Parameters as LP
import ClimaCore
import ClimaComms
ClimaComms.@import_required_backends

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

@testset "US-MOz forcing, data dates and comparison data" begin
    site_ID = "US-MOz"
    site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)
    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(site_ID_val))
    (; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))
    toml_dict = LP.create_toml_dict(FT)

    # Dates of the available data
    start_date, stop_date =
        FluxnetSimulations.get_data_dates(site_ID, time_offset)
    @test start_date < stop_date
    offset_start, offset_stop = FluxnetSimulations.get_data_dates(
        site_ID,
        time_offset;
        duration = Day(10),
        start_offset = Day(1),
    )
    @test offset_start == start_date + Day(1)
    @test offset_stop == offset_start + Day(10)
    @test_throws AssertionError FluxnetSimulations.get_data_dates(
        site_ID,
        time_offset;
        start_offset = -Day(1),
    )
    @test_throws AssertionError FluxnetSimulations.get_data_dates(
        site_ID,
        time_offset;
        duration = Year(100),
    )
    @test FluxnetSimulations.get_data_dt(site_ID) == 1800.0

    # Forcing constructed from the site data
    (; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
        site_ID,
        lat,
        long,
        time_offset,
        atmos_h,
        start_date,
        toml_dict,
        FT,
    )
    @test atmos isa ClimaLand.PrescribedAtmosphere
    @test atmos.h == atmos_h
    @test radiation isa ClimaLand.PrescribedRadiativeFluxes
    domain = ClimaLand.Domains.Column(;
        zlim = FT.((-1.0, 0.0)),
        nelements = 5,
        longlat = (long, lat),
    )
    surface_space = domain.space.surface
    function evaluate(input, t)
        field = ClimaCore.Fields.zeros(surface_space)
        evaluate!(field, input, t)
        return first(parent(field))
    end
    t = 3 * 3600.0
    @test 220 < evaluate(atmos.T, t) < 330
    @test 5e4 < evaluate(atmos.P, t) < 1.1e5
    @test 0 <= evaluate(atmos.q, t) < 0.05
    @test evaluate(atmos.u, t) >= 0
    @test evaluate(atmos.liquid_precip, t) <= 0
    @test evaluate(atmos.snow_precip, t) <= 0
    @test 3e-4 < evaluate(atmos.c_co2, t) < 6e-4
    @test evaluate(radiation.SW_d, t) >= 0
    @test evaluate(radiation.LW_d, t) > 0

    # Without precipitation partitioning, all precipitation is liquid
    (; atmos) = FluxnetSimulations.prescribed_forcing_fluxnet(
        site_ID,
        lat,
        long,
        time_offset,
        atmos_h,
        start_date,
        toml_dict,
        FT;
        split_precip = false,
    )
    @test evaluate(atmos.snow_precip, t) == 0

    # Observations for comparison with model output
    comparison = FluxnetSimulations.get_comparison_data(site_ID, time_offset)
    n = length(comparison.UTC_datetime)
    @test n > 0
    for name in (:gpp, :lhf, :shf, :swu, :lwu, :swc, :tsoil)
        @test haskey(comparison, name)
        @test length(getproperty(comparison, name)) == n
        @test !any(==(-9999), getproperty(comparison, name))
    end
    @test all(comparison.swc .<= 1)
    @test all(comparison.tsoil .> 200)
end

@testset "US-MOz initial conditions from observations" begin
    site_ID = "US-MOz"
    site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)
    (; time_offset, lat, long) =
        FluxnetSimulations.get_location(FT, Val(site_ID_val))
    (; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))
    toml_dict = LP.create_toml_dict(FT)
    start_date, _ = FluxnetSimulations.get_data_dates(site_ID, time_offset)
    forcing = FluxnetSimulations.prescribed_forcing_fluxnet(
        site_ID,
        lat,
        long,
        time_offset,
        atmos_h,
        start_date,
        toml_dict,
        FT,
    )
    domain = ClimaLand.Domains.Column(;
        zlim = FT.((-2.0, 0.0)),
        nelements = 10,
        longlat = (long, lat),
    )
    LAI = TimeVaryingInput((t) -> FT(1.0))
    land = ClimaLand.LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        FT(900);
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
    )
    set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
        site_ID,
        start_date,
        time_offset,
        land,
    )
    Y, p, _ = ClimaLand.initialize(land)
    set_ic!(Y, p, 0.0, land)

    (; ν, θ_r) = land.soil.parameters
    ϑ_l = parent(Y.soil.ϑ_l)
    @test all(ϑ_l .> parent(θ_r))
    @test all(ϑ_l .<= parent(ν))
    @test all(parent(Y.soil.θ_i) .== 0)
    T_soil = ClimaLand.Soil.temperature_from_ρe_int.(
        Y.soil.ρe_int,
        Y.soil.θ_i,
        ClimaLand.Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            land.soil.parameters.ρc_ds,
            land.soil.parameters.earth_param_set,
        ),
        land.soil.parameters.earth_param_set,
    )
    @test all(250 .< parent(T_soil) .< 320)
    @test all(parent(Y.snow.S) .== 0)
    @test all(parent(Y.soilco2.SOC) .== 5)
    @test all(parent(Y.canopy.hydraulics.ϑ_l) .> 0)
    @test all(250 .< parent(Y.canopy.energy.T) .< 320)
end
