using Test
using UnPack
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM.DETECT
import ClimaLSM

@testset "Soil CO2 production and transport" begin
    FT = Float64

    Pz = FT(101.0) # 1 [kPa] pressure just above the soil surface at time t
    Cᵣ = FT(10.0) # 2 [mg C cm-3] Cᵣ is root biomass carbon, see figure S5
    Csom = FT(5.0) # 3 [mg C cm-3] soil organic C content at depth z
    Cmic = FT(1.0) # 4 [mg C cm-3] Microbial C pool, ~ 1 % of Csom at DETECT site
    Rᵦ = FT(6e-5) # 5 [mg C cm-2] total root biomass C in a 1 m deep by 1 cm2 soil column
    α₁ᵣ = FT(11.65) # 6 [-]
    α₂ᵣ = FT(20.7) # 7 [-]
    α₃ᵣ = FT(-164.2) # 8 [-]
    Sᶜ = FT(711.6) # 9 [mg C cm-2] total soil organic C in a 1 m deep by 1 cm2 soil column
    Mᶜ = FT(12.3) # 10 [mg C cm-2] total microbial biomass C in a 1 m deep by 1 cm2 soil column
    Vᵦ = FT(0.0015) # 11 [mg C cm-3 h-1] value of Vmax at 10 °C and mean environmental conditions
    α₁ₘ = FT(14.05) # 12 [-]
    α₂ₘ = FT(11.05) # 13 [-]
    α₃ₘ = FT(-87.6) # 14 [-]
    Kₘ = FT(10e-5) # 15 [mg C cm-3 h-1] Michaelis-Menten half saturation constant
    CUE = FT(0.8) # 16 [mg C mg-1 C-1] microbial carbon use efficiency
    pf = FT(0.004) # 17 [-] fraction of soil organic C that is soluble
    Dₗᵢ = FT(3.17) # 18 [-] Diffusivity of soil C substrate in liquid
    E₀ₛ = FT(324.6) # 19 [Kelvin] temperature sensitivity parameter
    T₀ = FT(227.5) # 20 [Kelvin] temperature sensitivity-related parameter
    α₄ = FT(4.7) # 21 [-]
    Tᵣₑ = FT(283.0) # 22 [Kelvin] ref temperature for other param e.g., Rᵦ
    α₅ = FT(4.547) # 23 [-]
    BD = FT(1.12) # 24 [g cm-3] soil bulk density
    ϕ₁₀₀ = FT(0.1846) # 25 air filled porosity at soil water potential of -100 cm H2O (~ 10 kPa)
    PD = FT(2.52) # 26 [g cm-3] soil particle density
    Dstp = FT(1.39e-5) # 27 [m2 s-1] diffusion coefficient for CO2 in air at standard T and P
    P₀ = FT(101.325) # 28 [kPa] standard pressure
    b = FT(4.547) # 29 [-] parameter related to pore size distribution

    Tₛ = FT(303)
    θ = FT(0.3)
    θₐᵣ = FT(0.3)
    θₐₘ = FT(0.3)
    Tₛₐ = FT(303)
    

    parameters = DETECTParameters{FT}(Pz, Cᵣ, Csom, Cmic, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, pf, Dₗᵢ, E₀ₛ, T₀, α₄, Tᵣₑ, α₅, BD, ϕ₁₀₀, PD, Dstp, P₀, b)

    # Test that parameterizations functions are working properly
    @test diffusion_coefficient(Dstp, Tₛ, T₀, P₀, Pz) == Dstp * (Tₛ/T₀)^FT(1.75) * (P₀/Pz)
          Dg₀ = diffusion_coefficient(Dstp, Tₛ, T₀, P₀, Pz)
    @test air_filled_soil_porosity(BD, PD, θ) == (1 - BD/PD) - θ
          ϕ = air_filled_soil_porosity(BD, PD, θ)
    @test CO₂_diffusivity(Dg₀, ϕ₁₀₀, ϕ, b) == Dg₀ * (FT(2)ϕ₁₀₀^FT(3) + FT(0.04)ϕ₁₀₀) * (ϕ/ϕ₁₀₀)^(FT(2) + FT(3)/b)
    @test root_θ_adj(θ, θₐᵣ, α₁ᵣ, α₂ᵣ, α₃ᵣ) == exp(α₁ᵣ*θ + α₂ᵣ*θₐᵣ + α₃ᵣ*θ*θₐᵣ)
          fᵣ = root_θ_adj(θ, θₐᵣ, α₁ᵣ, α₂ᵣ, α₃ᵣ)
    @test energy_act(E₀ₛ, α₄, Tₛₐ) == E₀ₛ + α₄*Tₛₐ
          E₀ = energy_act(E₀ₛ, α₄, Tₛₐ)
    @test temp_adj(E₀, Tᵣₑ, T₀, Tₛ) == exp(E₀ * (FT(1)/(Tᵣₑ - T₀) - FT(1)/(Tₛ - T₀)))
          g = temp_adj(E₀, Tᵣₑ, T₀, Tₛ)
    @test root_source(Rᵦ, Cᵣ, fᵣ, g) == Rᵦ * Cᵣ * fᵣ * g
    @test microbe_θ_adj(α₁ₘ, α₂ₘ, α₃ₘ, θ, θₐₘ) == exp(α₁ₘ*θ + α₂ₘ*θₐₘ + α₃ₘ*θ*θₐₘ)
          fₘ = microbe_θ_adj(α₁ₘ, α₂ₘ, α₃ₘ, θ, θₐₘ)
    @test fVmax(Vᵦ, fₘ, g) == Vᵦ * fₘ * g
          Vmax = fVmax(Vᵦ, fₘ, g)
    @test fCsol(Csom, Dₗᵢ, θ, pf) == Csom * pf * θ^FT(3) * Dₗᵢ
          Csol = fCsol(Csom, Dₗᵢ, θ, pf)
    @test microbe_source(Vmax, Csol, Kₘ, Cmic, CUE) == Vmax * Csol/(Kₘ + Csol) * Cmic * (FT(1) - CUE)
end

