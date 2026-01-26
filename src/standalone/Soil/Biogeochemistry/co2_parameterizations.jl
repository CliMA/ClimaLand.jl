export volumetric_air_content,
    co2_diffusivity,
    microbe_source,
    o2_concentration,
    o2_fraction_from_concentration,
    o2_availability,
    henry_constant,
    beta_gas,
    effective_porosity


"""
    microbe_source(T_soil::FT,
                   θ_l::FT,
                   Csom::FT,
                   O2_avail::FT,
                   params::SoilCO2ModelParameters{FT}
                   ) where {FT}

Computes the CO₂ production in the soil by microbes, in depth and time (kg C / m^3/s), using
the Dual Arrhenius Michaelis Menten model (Davidson et al., 2012).
O2_avail is a dimensionless O₂ availability metric that accounts for tortuosity effects.
"""
function microbe_source(
    T_soil::FT,
    θ_l::FT,
    Csom::FT,
    O2_avail::FT,
    params::SoilCO2ModelParameters{FT},
) where {FT}
    (; α_sx, Ea_sx, kM_sx, kM_o2, D_liq, p_sx, earth_param_set) = params
    R = FT(LP.gas_constant(earth_param_set))
    Vmax = α_sx * exp(-Ea_sx / (R * T_soil)) # Maximum potential rate of respiration
    Sx = p_sx * Csom * D_liq * max(θ_l, FT(0))^3 # All soluble substrate, kgC m⁻³
    MM_sx = Sx / (kM_sx + Sx) # Availability of substrate factor, 0-1
    # Use pre-computed O2 availability (includes tortuosity effects)
    MM_o2 = O2_avail / (kM_o2 + O2_avail) # Oxygen limitation factor, 0-1
    R_sm = Vmax * MM_sx * MM_o2 # Respiration, kg C m⁻³ s⁻¹
    return R_sm
end


"""
    o2_concentration(O2_f::FT,
                     T_soil::FT,
                     P_sfc::FT,
                     params::SoilCO2ModelParameters{FT},
                     ) where {FT}

Computes the O₂ mass concentration in air (kg O2/m³ air) from the volumetric fraction O2_f,
using the ideal gas law.

The O2 mass concentration in the air phase is:
    ρ_O2_air = O2_f * P * M_O2 / (R * T)

where:
- O2_f : volumetric fraction of O2 in air (dimensionless, ~0.21)
- P: pressure (Pa)
- M_O2: molar mass of O2 (kg/mol) - from parameters
- R: universal gas constant (J/(mol·K)) - from ClimaParams
- T: temperature (K)

Note: This returns concentration per m³ of air, not per m³ of soil.
For diffusion in soil, the effective concentration per m³ of soil would be θ_a * ρ_O2_air,
but that multiplication is handled separately in the diffusion equation.
"""
function o2_concentration(
    O2_f::FT,
    T_soil::FT,
    P_sfc::FT,
    params::SoilCO2ModelParameters{FT},
) where {FT}
    R = FT(LP.gas_constant(params.earth_param_set))
    M_O2 = params.M_O2
    # O2 mass concentration in air phase (kg O2/m³ air)
    ρ_O2_air = O2_f * P_sfc * M_O2 / (R * T_soil)
    return ρ_O2_air
end


"""
    o2_fraction_from_concentration(ρ_O2_air::FT,
                                    T_soil::FT,
                                    P_sfc::FT,
                                    params::SoilCO2ModelParameters{FT},
                                    ) where {FT}

Computes the O₂ volumetric fraction (dimensionless) from the O₂ mass concentration in air,
using the ideal gas law. This is the inverse of o2_concentration.

The O2 volumetric fraction is:
    O2_f =ρ_O2_air * R * T / (P * M_O2)

where:
- ρ_O2_air: O2 mass concentration in air (kg O2/m³ air)
- P: pressure (Pa)
- M_O2: molar mass of O2 (kg/mol) - from parameters
- R: universal gas constant (J/(mol·K)) - from ClimaParams
- T: temperature (K)
"""
function o2_fraction_from_concentration(
    ρ_O2_air::FT,
    T_soil::FT,
    P_sfc::FT,
    params::SoilCO2ModelParameters{FT},
) where {FT}
    R = FT(LP.gas_constant(params.earth_param_set))
    M_O2 = params.M_O2
    # O2 volumetric fraction (dimensionless)
    O2_f = ρ_O2_air * R * T_soil / (P_sfc * M_O2)
    return O2_f
end


"""
    o2_availability(O2_f::FT,
                    θ_a::FT,
                    D_oa::FT,
                    ) where {FT}

Computes the dimensionless O₂ availability for microbial kinetics using
the Millington-Quirk tortuosity model.

The O2 availability accounts for diffusion limitations in porous media:
    O2_avail = D_oa * O2_f * θ_a^(4/3)

where:
- O2_f: volumetric fraction of O2 in air (dimensionless, ~0.21)
- θ_a: volumetric air content (m³ air / m³ soil)
- D_oa: oxygen diffusion coefficient in air (dimensionless)
- θ_a^(4/3): Millington-Quirk tortuosity factor

This is used in Michaelis-Menten kinetics for microbial respiration.
"""
function o2_availability(O2_f::FT, θ_a::FT, D_oa::FT) where {FT}
    # Dimensionless O2 availability with tortuosity
    O2_avail = D_oa * O2_f * θ_a^(FT(4 / 3))
    return O2_avail
end


"""
    henry_constant(K_H_298::FT, dln_K_H_dT::FT, T::FT) where {FT}

Compute temperature-dependent Henry's law constant using van 't Hoff equation.
Returns K_H in mol/(m³·Pa).

The temperature dependence follows:
    K_H(T) = K_H(T_ref) * exp[dln_K_H_dT * (1/T - 1/T_ref)]

where T_ref = 298.15 K and dln_K_H_dT is the temperature coefficient.

Reference: Sander (2015), Atmos. Chem. Phys., 15, 4399-4981.
"""
function henry_constant(K_H_298::FT, dln_K_H_dT::FT, T::FT) where {FT}
    T_ref = FT(298.15)
    return K_H_298 * exp(dln_K_H_dT * (FT(1) / T - FT(1) / T_ref))
end


"""
    beta_gas(K_H::FT, R::FT, T::FT) where {FT}

Compute dimensionless Henry's law factor β = K_H * R * T.

This converts liquid water storage to air-equivalent storage capacity.
For CO2 at 20-25°C: β ≈ 0.7-0.9 (significant buffering)
For O2 at 20°C: β ≈ 0.03 (less buffering, but still helps)

Arguments:
- K_H: Henry's law constant (mol/(m³·Pa))
- R: Universal gas constant (J/(mol·K))
- T: Temperature (K)
"""
function beta_gas(K_H::FT, R::FT, T::FT) where {FT}
    return K_H * R * T
end


"""
    effective_porosity(θ_a::FT, θ_l::FT, β::FT) where {FT}

Compute effective porosity accounting for gas and dissolved storage.

    θ_eff = θ_a + β * θ_l

When θ_a → 0 (saturated soil) but θ_l > 0, θ_eff remains finite,
preventing blow-up in concentration calculations.

Arguments:
- θ_a: Volumetric air content (m³/m³)
- θ_l: Volumetric liquid water content (m³/m³)
- β: Dimensionless Henry's law factor
"""
function effective_porosity(θ_a::FT, θ_l::FT, β::FT) where {FT}
    return θ_a + β * θ_l
end


"""
    volumetric_air_content(θ_w::FT,
                           ν::FT,
                           ) where {FT}

Computes the volumetric air content (`θ_a`) in the soil,
which is related to the total soil porosity (`ν`) and
volumetric soil water content (`θ_w = θ_l+θ_i`).

Note: The effective_porosity function provides numerical stability
when θ_a approaches zero by accounting for dissolved gas in liquid water.
"""
function volumetric_air_content(θ_w::FT, ν::FT) where {FT}
    θ_a = max(ν - θ_w, FT(0))  # Physical floor only
    return θ_a
end

"""
    co2_diffusivity(
                    T_soil::FT,
                    θ_w::FT,
                    P_sfc::FT,
                    θ_a100::FT,
                    b::FT,
                    ν::FT,
                    params::SoilCO2ModelParameters{FT},
                    ) where {FT}

Computes the diffusivity of CO₂ within the soil (D).

First, D0 is computed using the temperature within the soil (`T_soil` in K) and
pressure at the surface of the soil (`P_sfc` in Pa), using reference
values of `T_ref` and `P_ref` (273 K and 101325 Pa). Here, `θ_a` is the
volumetric air content and `θ_a100` is the volumetric air content
at a soil water potential of
100cm, and b is the pore size distribution of the soil.

This parameterization is from Ryan et al., GMD 11, 1909-1928, 2018,
https://doi.org/10.5194/gmd-11-1909-2018.
"""
function co2_diffusivity(
    T_soil::FT,
    θ_w::FT,
    P_sfc::FT,
    θ_a100::FT,
    b::FT,
    ν::FT,
    params::SoilCO2ModelParameters{FT},
) where {FT}
    (; D_ref, earth_param_set) = params
    T_ref = FT(LP.T_0(earth_param_set))
    P_ref = FT(LP.P_ref(earth_param_set))
    θ_a = volumetric_air_content(θ_w, ν)
    D0 = D_ref * max((T_soil / T_ref), 0)^FT(1.75) * (P_ref / P_sfc)
    D =
        D0 *
        (FT(2)θ_a100^FT(3) + FT(0.04)θ_a100) *
        (θ_a / θ_a100)^(FT(2) + FT(3) / b)
    return D
end
