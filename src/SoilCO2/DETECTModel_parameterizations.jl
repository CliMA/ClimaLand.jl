# define functions for DETECT auxiliary variables

export diffusion_coefficient,
       air_filled_soil_porosity,
       CO₂_diffusivity,
       root_θ_adj,
       temp_adj,
       energy_act,
       root_source,
       microbe_θ_adj,
       Vmax,
       Csol,
       microbe_source,
       S

# 1. CO2 Diffusivity
"""
	diffusion_coefficient(Dstp::FT, Ts::FT, T₀::FT, P₀::FT, P::FT) where {FT}

Diffusion coefficient for CO₂ in air at standard temperature (T₀, 273 K) and pressure 
(P₀, 101.325 kPa); Tₛ is the soil temperature (K) at depth z and time t, and P(t) is the
air pressure (kPa) just above the soil surface at time t. 
"""
function diffusion_coefficient(Dstp::FT, Tₛ::FT, T₀::FT, P₀::FT, Pz::FT) where {FT}
	Dg₀ = Dstp * (Tₛ/T₀)^FT(1.75) * (P₀/Pz)
	return Dg₀
end


"""
	air_filled_soil_porosity(BD::FT, PD::FT, θ::FT) where {FT}

Air filled soil porosity, which is related to the total soil porosity (ϕt) and 
volumetric soil water content (θ).
"""
function air_filled_soil_porosity(BD::FT, PD::FT, θ::FT) where {FT}
	ϕₜ = 1 - BD/PD # soil porosity
	ϕ = ϕₜ - θ
	return ϕ
end


"""
	CO₂_diffusivity(Dg₀::FT, ϕ₁₀₀::FT, ϕ::FT, b::FT) where {FT}

The diffusivity of CO₂ within the soil (Dgs).
"""
function CO2_diffusivity(Dg₀::FT, ϕ₁₀₀::FT, ϕ::FT, b::FT) where {FT}
	Dgs = Dg₀ * (FT(2)ϕ₁₀₀^FT(3) + FT(0.04)ϕ₁₀₀) * (ϕ/ϕ₁₀₀)^(FT(2) + FT(3)/b)
	return Dgs
end

# 2.1 CO2 source: root
"""
	root_θ_adj(α₁ᵣ::FT, α₂ᵣ::FT, α₃ᵣ::FT, θ::FT, θₐᵣ::FT) where {FT}

Soil moisture scaling function for roots, accounting for antecedent conditions. 
"""
function root_θ_adj(α₁ᵣ::FT, α₂ᵣ::FT, α₃ᵣ::FT, θ::FT, θₐᵣ::FT) where {FT}
	fᵣ = exp(α₁ᵣ*θ + α₂ᵣ*θₐᵣ + α₃ᵣ*θ*θₐᵣ)
	return fᵣ
end


"""
	energy_act(E₀ₛ::FT, α₄::FT, Tₛₐ::FT) where {FT}

Analogous to an energy of activation term that governs the apparent
temperature sensitivity of Sᵣ, accounting for antecedent condition 
of soil temperature. 
"""
function energy_act(E₀ₛ::FT, α₄::FT, Tₛₐ::FT) where {FT}
	E₀ = E₀ₛ + α₄*Tₛₐ
	return E₀
end


"""
	temp_adj(E₀::FT, Tᵣₑ::FT, T₀::FT, Tₛ::FT, T₀::FT) where {FT}

Temperature scaling function, motivated by Lloyd and Taylor (1994).
"""
function temp_adj(E₀::FT, Tᵣₑ::FT, T₀::FT, Tₛ::FT) where {FT}
	g = exp(E₀ * (FT(1)/(Tᵣₑ - T₀) - FT(1)/(Tₛ - T₀)))
	return g
end


"""
	root_source(Rᵦ::FT, Cᵣ::FT, fᵣ::FT, gᵣ::FT) where {FT}

CO₂ production in the soil by roots, in depth and time. 
"""
function root_source(Rᵦ::FT, Cᵣ::FT, fᵣ::FT, g::FT) where {FT}
	Sᵣ = Rᵦ * Cᵣ * fᵣ * g
	return Sᵣ
end

# 2.2 CO2 source: microbe
"""
	microbe_θ_adj(α₁ₘ::FT, α₂ₘ::FT, α₃ₘ::FT, θ::FT, θₐₘ::FT) where {FT}

Soil moisture scaling function for microbe, accounting for antecedent conditions.
"""
function microbe_θ_adj(α₁ₘ::FT, α₂ₘ::FT, α₃ₘ::FT, θ::FT, θₐₘ::FT) where {FT}
	fₘ = exp(α₁ₘ*θ + α₂ₘ*θₐₘ + α₃ₘ*θ*θₐₘ) 
	return fₘ
end


"""
	fVmax(Vᵦ::FT, fₘ::FT, g::FT) where {FT}

Maximum potential decomposition rate. 
"""
function fVmax(Vᵦ::FT, fₘ::FT, g::FT) where {FT}
	Vmax = Vᵦ * fₘ * g # g defined in 2.1
	return Vmax
end


"""
	fCsol(Csom::FT, p::FT, θ::FT, Dₗᵢ::FT) where {FT}

Soluble soil-C pool. 
"""
function fCsol(Csom::FT, pf::FT, θ::FT, Dₗᵢ::FT) where {FT}
	Csol = Csom * pf * θ^FT(3) * Dₗᵢ
	return Csol
end


"""
	microbe_source(Vmax::FT, Csol::FT, Km::FT, Cmic::FT, CUE::FT) where {FT}

CO₂ production in the soil by microbes, in depth and time.
"""
function microbe_source(Vmax::FT, Csol::FT, Km::FT, Cmic::FT, CUE::FT) where {FT}
	Sₘ = Vmax * Csol/(Kₘ + Csol) * Cmic * (FT(1) - CUE)
	return Sₘ
end


"""
	fS(Sᵣ::FT, Sₘ::FT) where {FT}

Total CO₂ production in the soil (roots + microbes).
"""
function fS(Sᵣ::FT, Sₘ::FT) where {FT}
	S = Sᵣ + Sₘ
	return S
end

