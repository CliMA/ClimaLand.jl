using UnPack

export volumetric_air_content,
    co2_diffusivity,
    root_source_moisture_coeff,
    source_temperature_coeff,
    energy_activation,
    root_source,
    microbe_source_moisture_coeff,
    decomposition_potential,
    soluble_soil_carbon,
    microbe_source

"""
    volumetric_air_content(θ_w::FT, params::DETECTModelParameters{FT}) where {FT}
    
Computes the volumetric air content (`θ_a`) in the soil, 
which is related to the total soil porosity (`ν`) and 
volumetric soil water content (`θ_w = θ_l+θ_i`).
"""
function volumetric_air_content(
    θ_w::FT,
    params::DETECTModelParameters{FT},
) where {FT}
    θ_a = params.ν - θ_w
    return θ_a
end


"""
    co2_diffusivity(
                    T_soil::FT,
                    P_sfc::FT,
                    θ_a::FT,
                    params::DETECTModelParameters{FT}
                    ) where {FT}

Computes the diffusivity of CO₂ within the soil (D).

First, D0 is computed using the temperature within the soil (`T_soil` in K) and 
pressure at the surface of the soil (`P_sfc` in Pa), using reference
values of `T_ref` and `P_ref` (273 K and 101325 Pa). Here, `θ_a` is the Volumetric
air content and `θ_a100` is the volumetric air content at a soil water potential of 
100cm, and b is the pore size distribution of the soil.
"""
function co2_diffusivity(
    T_soil::FT,
    P_sfc::FT,
    θ_w::FT,
    params::DETECTModelParameters{FT},
) where {FT}
    @unpack D_ref, T_ref, P_ref, θ_a100, b, ν = params
    θ_a = volumetric_air_content(θ_w, params)
    D0 = D_ref * (T_soil / T_ref)^FT(1.75) * (P_ref / P_sfc)
    D =
        D0 *
        (FT(2)θ_a100^FT(3) + FT(0.04)θ_a100) *
        (θ_a / θ_a100)^(FT(2) + FT(3) / b)
    return D
end


"""
    root_source_moisture_coeff(
                               θ_l::FT,
                               θ_ant_roots::FT,
                               params::DETECTModelParameters{FT}
                               ) where {FT}    

Returns the scaling coefficient/adjustment factor for the root co2 source term, 
accounting for antecedent "ant" conditions of soil moisture `θ_l`.
"""
function root_source_moisture_coeff(
    θ_l::FT,
    θ_ant_roots::FT,
    params::DETECTModelParameters{FT},
) where {FT}
    @unpack α1r, α2r, α3r = params
    coeff = exp(α1r * θ_l + α2r * θ_ant_roots + α3r * θ_l * θ_ant_roots)
    return coeff
end

# For microbe and root source terms

"""
    energy_activation(T_ant_soil::FT, params::DETECTModelParameters{FT}) where {FT}
    
Computes the energy activation temperature based on the antecendent soil temperature `T_ant_soil`.

This corresponds to an energy of activation term that governs the apparent
temperature sensitivity of the root and microbe source term.
"""
function energy_activation(
    T_ant_soil::FT,
    params::DETECTModelParameters{FT},
) where {FT}
    @unpack α4, Estar = params
    E0 = Estar + α4 * T_ant_soil
    return E0
end


"""
    source_temperature_coeff(
        T_soil::FT,
        T_ant_soil::FT,
        params::DETECTModelParameters{FT},
    ) where {FT}
    
Computes the temperature scaling factor function (unitless), which is motivated by Lloyd and Taylor (1994).
"""
function source_temperature_coeff(
    T_soil::FT,
    T_ant_soil::FT,
    params::DETECTModelParameters{FT},
) where {FT}
    @unpack T_ref, T_ref_soil = params
    E0 = energy_activation(T_ant_soil, params)
    g = exp(E0 * (FT(1) / (T_ref_soil - T_ref - FT(1) / (T_soil - T_ref))))
    return g
end


"""
    root_source(
                T_soil::FT,
                T_ant_soil::FT,
                θ_l::FT,
                θ_ant_roots::FT,
                Cr::FT,
                params::DETECTModelParameters{FT}
                ) where {FT}
                
Computes the CO₂ production in the soil by roots, in depth and time (kg C / m^3/s).
"""
function root_source(
    T_soil::FT,
    T_ant_soil::FT,
    θ_l::FT,
    θ_ant_roots::FT,
    Cr::FT,
    params::DETECTModelParameters{FT},
) where {FT}
    g = source_temperature_coeff(T_soil, T_ant_soil, params)
    coeff = root_source_moisture_coeff(θ_l, θ_ant_roots, params)
    source = params.Rb * Cr * coeff * g
    return source
end

# 2.2 CO2 source: microbe
"""
    microbe_source_moisture_coeff(
                                  θ_l::FT,
                                  θ_ant_microbe::FT,
                                  params::DETECTModelParameters{FT}
                                  ) where {FT}    

Returns the scaling coefficient/adjustment factor for the microbe co2 source term, 
accounting for antecedent "ant" conditions of soil moisture `θ_l`.
"""
function microbe_source_moisture_coeff(
    θ_l::FT,
    θ_ant_microbe::FT,
    params::DETECTModelParameters{FT},
) where {FT}
    @unpack α1m, α2m, α3m = params
    coeff = exp(α1m * θ_l + α2m * θ_ant_microbe + α3m * θ_l * θ_ant_microbe)
    return coeff
end


"""
    decomposition_potential(
                            T_soil::FT,
                            T_ant_soil::FT,
                            θ_l::FT,
                            θ_ant_microbe::FT,
                            params::DETECTModelParameters{FT}
                            ) where {FT}

Computes the maximum potential decomposition rate (kg C/m^3/s). 
"""
function decomposition_potential(
    T_soil::FT,
    T_ant_soil::FT,
    θ_l::FT,
    θ_ant_microbe::FT,
    params::DETECTModelParameters{FT},
) where {FT}
    g = source_temperature_coeff(T_soil, T_ant_soil, params)
    fm = microbe_source_moisture_coeff(θ_l, θ_ant_microbe, params)
    Vmax = params.Vb * fm * g # g defined in 2.1
    return Vmax
end


"""
    soluble_soil_carbon(
                        θ_l::FT,
                        Csom::FT,
                        params::DETECTModelParameters{FT}
                        ) where {FT}        
Computes the soluble soil-C pool (kg C/m^3).
"""
function soluble_soil_carbon(
    θ_l::FT,
    Csom::FT,
    params::DETECTModelParameters{FT},
) where {FT}
    @unpack soluble_fraction, D_liq = params
    Csol = Csom * soluble_fraction * θ_l^FT(3) * D_liq
    return Csol
end


"""
    microbe_source(
                   T_soil::FT,
                   T_ant_soil::FT,
                   θ_l::FT,
                   θ_ant_microbe::FT,
                   Csom::FT,
                   Cmic::FT,
                   params::DETECTModelParameters{FT},
                   ) where {FT}

Computes the CO₂ production in the soil by microbes, in depth and time (kg C / m^3/s).
"""
function microbe_source(
    T_soil::FT,
    T_ant_soil::FT,
    θ_l::FT,
    θ_ant_microbe::FT,
    Csom::FT,
    Cmic::FT,
    params::DETECTModelParameters{FT},
) where {FT}
    @unpack CUE, Km = params
    Vmax =
        decomposition_potential(T_soil, T_ant_soil, θ_l, θ_ant_microbe, params)
    Csol = soluble_soil_carbon(θ_l, Csom, params)
    Sm = Vmax * Csol / (Km + Csol) * Cmic * (FT(1) - CUE)
    return Sm
end
