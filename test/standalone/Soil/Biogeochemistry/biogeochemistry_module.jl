using Test
using Dates
import CLIMAParameters as CP
using ClimaCore
using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil.Biogeochemistry
import ClimaLSM
import ClimaLSM.Parameters as LSMP

include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

for FT in (Float32, Float64)
    @testset "Soil co2 biogeochemistry sources, FT = $FT" begin
        earth_param_set = create_lsm_parameters(FT)
        # Prognostic variables
        P_sfc = (t) -> 101e3
        T_soil = (z, t) -> eltype(z)(t)
        θ_l = (z, t) -> eltype(z)(0.3)
        θ_i = (z, t) -> eltype(z)(0)
        Csom = (z, t) -> eltype(z)(5.0) # 3 [kg C m-3] soil organic C content at depth z
        D_ref = FT(0.0)
        parameters = SoilCO2ModelParameters{FT}(;
            D_ref = D_ref,
            earth_param_set = earth_param_set,
        )

        nelems = 50 # number of layers in the vertical
        zmin = FT(-1) # 0 to 1 m depth
        zmax = FT(0.0)
        soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
        top_bc = SoilCO2StateBC((p, t) -> 3000.0)
        bot_bc = SoilCO2StateBC((p, t) -> 100.0)
        sources = (MicrobeProduction{FT}(),)
        boundary_conditions =
            (; top = (CO2 = top_bc,), bottom = (CO2 = bot_bc,))

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

        atmos = ClimaLSM.PrescribedAtmosphere(
            precipitation_function,
            snow_precip,
            atmos_T,
            atmos_u,
            atmos_q,
            atmos_p,
            UTC_DATETIME,
            atmos_h;
            c_co2 = atmos_co2,
        )

        soil_drivers = SoilDrivers(
            PrescribedMet{FT}(T_soil, θ_l),
            PrescribedSOC{FT}(Csom),
            atmos,
        )

        model = SoilCO2Model{FT}(;
            parameters = parameters,
            domain = soil_domain,
            sources = sources,
            boundary_conditions = boundary_conditions,
            drivers = soil_drivers,
        )

        Y, p, coords = initialize(model)
        dY = similar(Y)
        exp_tendency! = make_exp_tendency(model)
        t = Float64(1)
        Y.soilco2.C .= FT(4)
        set_initial_aux_state! = make_set_initial_aux_state(model)
        set_initial_aux_state!(p, Y, t)
        exp_tendency!(dY, Y, p, t)
        @test dY.soilco2.C ≈ p.soilco2.Sm
    end


    @testset "Soil co2 biogeochemistry diffusion, FT = $FT" begin
        earth_param_set = create_lsm_parameters(FT)
        # Prognostic variables
        P_sfc = (t) -> 101e3
        T_soil = (z, t) -> eltype(z)(303)
        θ_l = (z, t) -> eltype(z)(0.3)
        θ_i = (z, t) -> eltype(z)(0.0)
        Csom = (z, t) -> eltype(z)(5.0) # 3 [kg C m-3] soil organic C content at depth z

        parameters =
            SoilCO2ModelParameters{FT}(; earth_param_set = earth_param_set)
        C = FT(4)
        nelems = 50 # number of layers in the vertical
        zmin = FT(-1) # 0 to 1 m depth
        zmax = FT(0.0)
        soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
        top_bc = SoilCO2StateBC((p, t) -> C)
        bot_bc = SoilCO2StateBC((p, t) -> C)
        sources = ()
        boundary_conditions =
            (; top = (CO2 = top_bc,), bottom = (CO2 = bot_bc,))

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

        atmos = ClimaLSM.PrescribedAtmosphere(
            precipitation_function,
            snow_precip,
            atmos_T,
            atmos_u,
            atmos_q,
            atmos_p,
            UTC_DATETIME,
            atmos_h;
            c_co2 = atmos_co2,
        )

        soil_drivers = SoilDrivers(
            PrescribedMet{FT}(T_soil, θ_l),
            PrescribedSOC{FT}(Csom),
            atmos, # need to create some functions
        )

        model = SoilCO2Model{FT}(;
            parameters = parameters,
            domain = soil_domain,
            sources = sources,
            boundary_conditions = boundary_conditions,
            drivers = soil_drivers,
        )

        Y, p, coords = initialize(model)
        dY = similar(Y)
        exp_tendency! = make_exp_tendency(model)
        t = Float64(1)
        Y.soilco2.C .= FT(C)
        set_initial_aux_state! = make_set_initial_aux_state(model)
        set_initial_aux_state!(p, Y, t)
        exp_tendency!(dY, Y, p, t)
        @test sum(dY.soilco2.C) ≈ FT(0.0)
    end
end
