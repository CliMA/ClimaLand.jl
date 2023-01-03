using Test
using UnPack
import CLIMAParameters as CP

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using DiffEqCallbacks
using OrdinaryDiffEq: ODEProblem, solve, Euler
using ClimaCore
using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil.Biogeochemistry
import ClimaLSM
import ClimaLSM.Parameters as LSMP

include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

@testset "Soil co2 biogeochemistry sources" begin
    FT = Float32
    earth_param_set = create_lsm_parameters(FT)
    # Prognostic variables
    T_soil = (z, t) -> eltype(t)(303)
    θ_l = (z, t) -> eltype(t)(0.3)
    θ_i = (z, t) -> eltype(t)(0.0)
    θ_ant_roots = (z, t) -> eltype(t)(0.3)
    θ_ant_microbe = (z, t) -> eltype(t)(0.3)
    T_ant_soil = (z, t) -> eltype(t)(303)
    Cr = (z, t) -> eltype(t)(10.0) # 2 [kg C m-3] Cr is root biomass carbon, see figure S5
    Csom = (z, t) -> eltype(t)(5.0) # 3 [kg C m-3] soil organic C content at depth z
    Cmic = (z, t) -> eltype(t)(1.0) # 4 [kg C m-3] Microbial C pool, ~ 1 % of Csom at DETECT site
    D_ref = FT(0.0)
    parameters = SoilCO2ModelParameters{FT}(;D_ref = D_ref, earth_param_set = earth_param_set)

    nelems = 50 # number of layers in the vertical
    zmin = FT(-1) # 0 to 1 m depth
    zmax = FT(0.0)
    soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
    top_bc = SoilCO2StateBC((p, t) -> eltype(t)(3000.0))
    bot_bc = SoilCO2StateBC((p, t) -> eltype(t)(100.0))
    sources = (RootProduction{FT}(), MicrobeProduction{FT}())
    boundary_conditions = (; CO2 = (top = top_bc, bottom = bot_bc))

    soil_drivers = PrescribedSoil(
        T_soil,
        θ_l,
        T_ant_soil,
        θ_ant_microbe,
        θ_ant_roots,
        Cr,
        Csom,
        Cmic,
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
    ode! = make_ode_function(model)
    t = FT(1)
    ode!(dY, Y, p, t)
    @test dY.soilco2.C ≈ p.soilco2.Sm .+ p.soilco2.Sr
end


@testset "Soil co2 biogeochemistry sources" begin
    FT = Float32
    earth_param_set = create_lsm_parameters(FT)
    # Prognostic variables
    T_soil = (z, t) -> eltype(t)(303)
    θ_l = (z, t) -> eltype(t)(0.3)
    θ_i = (z, t) -> eltype(t)(0.0)
    θ_ant_roots = (z, t) -> eltype(t)(0.3)
    θ_ant_microbe = (z, t) -> eltype(t)(0.3)
    T_ant_soil = (z, t) -> eltype(t)(303)
    Cr = (z, t) -> eltype(t)(10.0) # 2 [kg C m-3] Cr is root biomass carbon, see figure S5
    Csom = (z, t) -> eltype(t)(5.0) # 3 [kg C m-3] soil organic C content at depth z
    Cmic = (z, t) -> eltype(t)(1.0) # 4 [kg C m-3] Microbial C pool, ~ 1 % of Csom at DETECT site

    parameters = SoilCO2ModelParameters{FT}(;earth_param_set = earth_param_set)
    C = FT(4)
    nelems = 50 # number of layers in the vertical
    zmin = FT(-1) # 0 to 1 m depth
    zmax = FT(0.0)
    soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
    top_bc = SoilCO2StateBC((p, t) -> eltype(t)(C))
    bot_bc = SoilCO2StateBC((p, t) -> eltype(t)(C))
    sources = ()
    boundary_conditions = (; CO2 = (top = top_bc, bottom = bot_bc))

    soil_drivers = PrescribedSoil(
        T_soil,
        θ_l,
        T_ant_soil,
        θ_ant_microbe,
        θ_ant_roots,
        Cr,
        Csom,
        Cmic,
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
    ode! = make_ode_function(model)
    t = FT(1)
    Y.soilco2.C .= FT(C)
    ode!(dY, Y, p, t)
    @test sum(dY.soilco2.C) ≈ FT(0.0)
end
