using Test
using UnPack
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using DiffEqCallbacks
using OrdinaryDiffEq: ODEProblem, solve, Euler
using ClimaCore
using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil.Biogeochemistry

@testset "Soil co2 biogeochemistry sources" begin
    FT = Float32

    # Parameters should be supplied in m/kg/s (Pa... etc)
    P_sfc = FT(101e3) # 1 [Pa] pressure just above the soil surface at time t
    Rb = FT(6e-5) # 5 [kg m^3] total root biomass C in soil column
    α1r = FT(11.65) # 6 [-]
    α2r = FT(20.7) # 7 [-]
    α3r = FT(-164.2) # 8 [-]
    Vb = FT(0.0015) # 11 [kg C m-3 s-1] value of Vmax at 10 °C and mean environmental conditions
    α1m = FT(14.05) # 12 [-]
    α2m = FT(11.05) # 13 [-]
    α3m = FT(-87.6) # 14 [-]
    Km = FT(10e-5) # 15 [kg C m-3 s-1] Michaelis-Menten half saturation constant
    CUE = FT(0.8) # 16 [kg C kg-1 C-1] microbial carbon use efficiency
    soluble_fraction = FT(0.004) # 17 [-] fraction of soil organic C that is soluble
    D_liq = FT(3.17) # 18 [-] Diffusivity of soil C substrate in liquid
    Estar = FT(324.6) # 19 [Kelvin] temperature sensitivity parameter
    T_ref = FT(273.15) # 20 [Kelvin] temperature sensitivity-related parameter
    α4 = FT(-4.7) # 21 [-]
    T_ref_soil = FT(283.15) # 22 [Kelvin] ref temperature for other param e.g., Rb
    α5 = FT(4.547) # 23 [-]
    ν = FT(0.556) # 26 [m3 m-3] soil porosity
    θ_a100 = FT(0.1846) # 25 air filled porosity at soil water potential of -100 cm H2O
    D_ref = FT(0)#FT(1.39e-5) # 27 [m2 s-1] diffusion coefficient for CO2 in air at standard T and P
    P_ref = FT(101325) # 28 [Pa] standard pressure
    b = FT(4.547) # 29 [-] parameter related to pore size distribution

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

    parameters = SoilCO2ModelParameters(;
        P_sfc = P_sfc,
        Rb = Rb,
        α1r = α1r,
        α2r = α2r,
        α3r = α3r,
        Vb = Vb,
        α1m = α1m,
        α2m = α2m,
        α3m = α3m,
        Km = Km,
        CUE = CUE,
        soluble_fraction = soluble_fraction,
        D_liq = D_liq,
        Estar = Estar,
        T_ref = T_ref,
        α4 = α4,
        T_ref_soil = T_ref_soil,
        α5 = α5,
        ν = ν,
        θ_a100 = θ_a100,
        D_ref = D_ref,
        P_ref = P_ref,
        b = b,
    )

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

    # Parameters should be supplied in m/kg/s (Pa... etc)
    P_sfc = FT(101e3) # 1 [Pa] pressure just above the soil surface at time t
    Rb = FT(6e-5) # 5 [kg m^3] total root biomass C in soil column
    α1r = FT(11.65) # 6 [-]
    α2r = FT(20.7) # 7 [-]
    α3r = FT(-164.2) # 8 [-]
    Vb = FT(0.0015) # 11 [kg C m-3 s-1] value of Vmax at 10 °C and mean environmental conditions
    α1m = FT(14.05) # 12 [-]
    α2m = FT(11.05) # 13 [-]
    α3m = FT(-87.6) # 14 [-]
    Km = FT(10e-5) # 15 [kg C m-3 s-1] Michaelis-Menten half saturation constant
    CUE = FT(0.8) # 16 [kg C kg-1 C-1] microbial carbon use efficiency
    soluble_fraction = FT(0.004) # 17 [-] fraction of soil organic C that is soluble
    D_liq = FT(3.17) # 18 [-] Diffusivity of soil C substrate in liquid
    Estar = FT(324.6) # 19 [Kelvin] temperature sensitivity parameter
    T_ref = FT(273.15) # 20 [Kelvin] temperature sensitivity-related parameter
    α4 = FT(-4.7) # 21 [-]
    T_ref_soil = FT(283.15) # 22 [Kelvin] ref temperature for other param e.g., Rb
    α5 = FT(4.547) # 23 [-]
    ν = FT(0.556) # 26 [m3 m-3] soil porosity
    θ_a100 = FT(0.1846) # 25 air filled porosity at soil water potential of -100 cm H2O
    D_ref = FT(1.39e-5) # 27 [m2 s-1] diffusion coefficient for CO2 in air at standard T and P
    P_ref = FT(101325) # 28 [Pa] standard pressure
    b = FT(4.547) # 29 [-] parameter related to pore size distribution

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

    parameters = SoilCO2ModelParameters(;
        P_sfc = P_sfc,
        Rb = Rb,
        α1r = α1r,
        α2r = α2r,
        α3r = α3r,
        Vb = Vb,
        α1m = α1m,
        α2m = α2m,
        α3m = α3m,
        Km = Km,
        CUE = CUE,
        soluble_fraction = soluble_fraction,
        D_liq = D_liq,
        Estar = Estar,
        T_ref = T_ref,
        α4 = α4,
        T_ref_soil = T_ref_soil,
        α5 = α5,
        ν = ν,
        θ_a100 = θ_a100,
        D_ref = D_ref,
        P_ref = P_ref,
        b = b,
    )
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
