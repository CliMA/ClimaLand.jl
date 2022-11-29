using Test
using UnPack
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM.Soil.Biogeochemistry

@testset "Soil CO2 production and transport" begin
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
    α4 = FT(4.7) # 21 [-]
    T_ref_soil = FT(283.15) # 22 [Kelvin] ref temperature for other param e.g., Rb
    α5 = FT(4.547) # 23 [-]
    ν = FT(0.556) # 26 [m3 m-3] soil porosity
    θ_a100 = FT(0.1846) # 25 air filled porosity at soil water potential of -100 cm H2O
    D_ref = FT(1.39e-5) # 27 [m2 s-1] diffusion coefficient for CO2 in air at standard T and P
    P_ref = FT(101325) # 28 [Pa] standard pressure
    b = FT(4.547) # 29 [-] parameter related to pore size distribution

    # Prognostic variables
    T_soil = FT(303)
    θ_l = FT(0.3)
    θ_i = FT(0.0)
    θ_w = θ_l + θ_i
    θ_ant_roots = FT(0.3)
    θ_ant_microbe = FT(0.3)
    T_ant_soil = FT(303)
    Cr = FT(10.0) # 2 [kg C m-3] Cr is root biomass carbon, see figure S5
    Csom = FT(5.0) # 3 [kg C m-3] soil organic C content at depth z
    Cmic = FT(1.0) # 4 [kg C m-3] Microbial C pool, ~ 1 % of Csom at DETECT site

    parameters = DETECTModelParameters(;
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

    # Test that parameterizations functions are working properly
    θ_a = volumetric_air_content(θ_w, parameters)
    @test θ_a == parameters.ν - θ_w
    @test typeof(θ_a) == FT

    D = co2_diffusivity(T_soil, P_sfc, θ_w, parameters)
    @test D ==
          (
              parameters.D_ref *
              (T_soil / T_ref)^FT(1.75) *
              (parameters.P_ref / parameters.P_sfc)
          ) *
          (FT(2)parameters.θ_a100^FT(3) + FT(0.04)parameters.θ_a100) *
          (θ_a / parameters.θ_a100)^(FT(2) + FT(3) / parameters.b)
    @test typeof(D) == FT

    fr = root_source_moisture_coeff(θ_l, θ_ant_roots, parameters)
    @test fr == exp(α1r * θ_l + α2r * θ_ant_roots + α3r * θ_l * θ_ant_roots)
    @test typeof(fr) == FT

    E0 = energy_activation(T_ant_soil, parameters)
    @test E0 == Estar + α4 * T_ant_soil
    @test typeof(E0) == FT

    g = source_temperature_coeff(T_soil, T_ant_soil, parameters)
    @test g ==
          exp(E0 * (FT(1) / (T_ref_soil - T_ref - FT(1) / (T_soil - T_ref))))
    @test typeof(g) == FT

    source = root_source(T_soil, T_ant_soil, θ_l, θ_ant_roots, Cr, parameters)
    @test source == Rb * Cr * fr * g
    @test typeof(source) == FT

    fm = microbe_source_moisture_coeff(θ_l, θ_ant_microbe, parameters)
    @test fm == exp(α1m * θ_l + α2m * θ_ant_microbe + α3m * θ_l * θ_ant_microbe)
    @test typeof(fm) == FT

    Vmax = decomposition_potential(
        T_soil,
        T_ant_soil,
        θ_l,
        θ_ant_microbe,
        parameters,
    )
    @test Vmax == Vb * fm * g
    @test typeof(Vmax) == FT

    Csol = soluble_soil_carbon(θ_l, Csom, parameters)
    @test eltype(Csol) == FT
    @test Csol == Csom * soluble_fraction * θ_l^FT(3) * D_liq

    ms = microbe_source(
        T_soil,
        T_ant_soil,
        θ_l,
        θ_ant_microbe,
        Csom,
        Cmic,
        parameters,
    )
    @test ms == Vmax * Csol / (Km + Csol) * Cmic * (FT(1) - CUE)
    @test typeof(ms) == FT
end
