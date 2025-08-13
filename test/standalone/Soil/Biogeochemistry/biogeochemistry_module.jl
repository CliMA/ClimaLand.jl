using Test
import ClimaComms
ClimaComms.@import_required_backends
using Dates
import ClimaParams as CP
using ClimaCore
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil.Biogeochemistry
import ClimaLand
import ClimaLand.Parameters as LP

for FT in (Float32, Float64)
    @testset "Soil co2 biogeochemistry sources, FT = $FT" begin
        # Prognostic variables
        P_sfc = (t) -> 101e3
        T_soil = (z, t) -> eltype(z)(t)
        θ_l = (z, t) -> eltype(z)(0.3)
        θ_i = (z, t) -> eltype(z)(0)
        Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(
            TimeVaryingInput((t) -> 5),
        )
        D_ref = FT(0.0)
        parameters = SoilCO2ModelParameters(FT; D_ref)

        nelems = 50 # number of layers in the vertical
        zmin = FT(-1) # 0 to 1 m depth
        zmax = FT(0.0)
        soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
        top_bc = SoilCO2StateBC((p, t) -> 3000.0)
        bot_bc = SoilCO2StateBC((p, t) -> 100.0)
        sources = (MicrobeProduction{FT}(),)
        boundary_conditions = (; top = top_bc, bottom = bot_bc)

        # Make a PrescribedAtmosphere - we only care about atmos_p though
        precipitation_function = (t) -> 1.0
        snow_precip = (t) -> 1.0
        atmos_T = (t) -> 1.0
        atmos_u = (t) -> 1.0
        atmos_q = (t) -> 1.0
        atmos_p = (t) -> 100000.0
        UTC_DATETIME = Dates.now()
        atmos_h = FT(30)
        atmos_co2 = (t) -> 1.0
        earth_param_set = LP.LandParameters(FT)

        atmos = ClimaLand.PrescribedAtmosphere(
            TimeVaryingInput(precipitation_function),
            TimeVaryingInput(snow_precip),
            TimeVaryingInput(atmos_T),
            TimeVaryingInput(atmos_u),
            TimeVaryingInput(atmos_q),
            TimeVaryingInput(atmos_p),
            UTC_DATETIME,
            atmos_h,
            earth_param_set;
            c_co2 = TimeVaryingInput(atmos_co2),
        )
        ν = FT(0.6)
        θ_r = FT(0.0)
        α = FT(0.1)
        n = FT(2)
        hcm = ClimaLand.Soil.vanGenuchten{FT}(; α = α, n = n)
        prescribed_met = PrescribedMet{FT}(T_soil, θ_l, ν, θ_r, hcm)
        soil_drivers = SoilDrivers(prescribed_met, Csom, atmos)

        model = SoilCO2Model{FT}(
            soil_domain,
            soil_drivers;
            parameters,
            boundary_conditions,
            sources,
        )

        Y, p, coords = initialize(model)
        dY = similar(Y)
        exp_tendency! = make_exp_tendency(model)
        t = Float64(1)
        Y.soilco2.C .= FT(4)
        set_initial_cache! = make_set_initial_cache(model)
        set_initial_cache!(p, Y, t)
        exp_tendency!(dY, Y, p, t)
        @test dY.soilco2.C ≈ p.soilco2.Sm
        @test p.soilco2.top_bc_wvec ==
              ClimaCore.Geometry.WVector.(p.soilco2.top_bc)
        @test p.soilco2.bottom_bc_wvec ==
              ClimaCore.Geometry.WVector.(p.soilco2.bottom_bc)
    end


    @testset "Soil co2 biogeochemistry diffusion, FT = $FT" begin
        # Prognostic variables
        P_sfc = (t) -> 101e3
        T_soil = (z, t) -> eltype(z)(303)
        θ_l = (z, t) -> eltype(z)(0.3)
        θ_i = (z, t) -> eltype(z)(0.0)
        Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(
            TimeVaryingInput((t) -> 5),
        )

        parameters = SoilCO2ModelParameters(FT)
        C = FT(4)
        nelems = 50 # number of layers in the vertical
        zmin = FT(-1) # 0 to 1 m depth
        zmax = FT(0.0)
        soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
        top_bc = SoilCO2StateBC((p, t) -> C)
        bot_bc = SoilCO2StateBC((p, t) -> C)
        sources = ()
        boundary_conditions = (; top = top_bc, bottom = bot_bc)

        # Make a PrescribedAtmosphere - we only care about atmos_p though
        precipitation_function = (t) -> 1.0
        snow_precip = (t) -> 1.0
        atmos_T = (t) -> 1.0
        atmos_u = (t) -> 1.0
        atmos_q = (t) -> 1.0
        atmos_p = (t) -> 100000.0
        UTC_DATETIME = Dates.now()
        atmos_h = FT(30)
        atmos_co2 = (t) -> 1.0
        earth_param_set = LP.LandParameters(FT)

        atmos = ClimaLand.PrescribedAtmosphere(
            TimeVaryingInput(precipitation_function),
            TimeVaryingInput(snow_precip),
            TimeVaryingInput(atmos_T),
            TimeVaryingInput(atmos_u),
            TimeVaryingInput(atmos_q),
            TimeVaryingInput(atmos_p),
            UTC_DATETIME,
            atmos_h,
            earth_param_set;
            c_co2 = TimeVaryingInput(atmos_co2),
        )
        ν = FT(0.6)
        θ_r = FT(0.0)
        α = FT(0.1)
        n = FT(2)
        hcm = ClimaLand.Soil.vanGenuchten{FT}(; α = α, n = n)
        soil_drivers = SoilDrivers(
            PrescribedMet{FT}(T_soil, θ_l, ν, θ_r, hcm),
            Csom,
            atmos, # need to create some functions
        )

        model = SoilCO2Model{FT}(
            soil_domain,
            soil_drivers;
            parameters,
            boundary_conditions,
            sources,
        )

        Y, p, coords = initialize(model)
        dY = similar(Y)
        exp_tendency! = make_exp_tendency(model)
        t = Float64(1)
        Y.soilco2.C .= FT(C)
        set_initial_cache! = make_set_initial_cache(model)
        set_initial_cache!(p, Y, t)
        exp_tendency!(dY, Y, p, t)
        @test sum(dY.soilco2.C) ≈ FT(0.0)
        @test p.soilco2.top_bc_wvec ==
              ClimaCore.Geometry.WVector.(p.soilco2.top_bc)
        @test p.soilco2.bottom_bc_wvec ==
              ClimaCore.Geometry.WVector.(p.soilco2.bottom_bc)
    end
end
