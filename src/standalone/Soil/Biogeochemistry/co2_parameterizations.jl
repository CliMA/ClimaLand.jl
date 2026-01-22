export volumetric_air_content,
    co2_diffusivity,
    microbe_source,
    o2_concentration,
    o2_fraction_from_concentration,
    o2_availability


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
    # Clamp Csom and O2_avail to non-negative to prevent negative respiration
    # when prognostic variables slightly undershoot zero
    Csom_clamped = max(Csom, FT(0))
    O2_avail_clamped = max(O2_avail, FT(0))
    Sx = p_sx * Csom_clamped * D_liq * max(θ_l, FT(0))^3 # All soluble substrate, kgC m⁻³
    MM_sx = Sx / (kM_sx + Sx) # Availability of substrate factor, 0-1
    # Use pre-computed O2 availability (includes tortuosity effects)
    MM_o2 = O2_avail_clamped / (kM_o2 + O2_avail_clamped) # Oxygen limitation factor, 0-1
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
    # Clamp inputs to non-negative to prevent NaN from fractional exponent
    # when prognostic variables slightly undershoot zero
    O2_f_clamped = max(O2_f, FT(0))
    θ_a_clamped = max(θ_a, FT(0))
    # Dimensionless O2 availability with tortuosity
    O2_avail = D_oa * O2_f_clamped * θ_a_clamped^(FT(4 / 3))
    return O2_avail
end


"""
    volumetric_air_content(θ_w::FT,
                           ν::FT,
                           ) where {FT}

Computes the volumetric air content (`θ_a`) in the soil,
which is related to the total soil porosity (`ν`) and
volumetric soil water content (`θ_w = θ_l+θ_i`).
"""
function volumetric_air_content(θ_w::FT, ν::FT) where {FT}
    θ_a = max(ν - θ_w, FT(0))
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
    # Clamp the ratio to non-negative to prevent NaN from fractional exponent
    # when θ_a is very small or θ_a100 has numerical issues
    θ_a_ratio = max(θ_a / max(θ_a100, eps(FT)), FT(0))
    D =
        D0 *
        (FT(2)θ_a100^FT(3) + FT(0.04)θ_a100) *
        θ_a_ratio^(FT(2) + FT(3) / b)
    return D
end
