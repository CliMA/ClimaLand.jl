export volumetric_air_content, co2_diffusivity, microbe_source, o2_concentration


"""
    microbe_source(T_soil::FT,
                   θ_l::FT,
                   Csom::FT,
                   O2_a::FT,
                   ν::FT,
                   params::SoilCO2ModelParameters{FT}
                   ) where {FT}

Computes the CO₂ production in the soil by microbes, in depth and time (kg C / m^3/s), using
the Dual Arrhenius Michaelis Menten model (Davidson et al., 2012).
O2_a is the prognostic volumetric fraction of O₂ in the soil air.
"""
function microbe_source(
    T_soil::FT,
    θ_l::FT,
    Csom::FT,
    O2_a::FT,
    ν::FT,
    params::SoilCO2ModelParameters{FT},
) where {FT}
    (; α_sx, Ea_sx, kM_sx, kM_o2, D_liq, p_sx, D_oa, earth_param_set) =
        params
    R = FT(LP.gas_constant(earth_param_set))
    Vmax = α_sx * exp(-Ea_sx / (R * T_soil)) # Maximum potential rate of respiration
    Sx = p_sx * Csom * D_liq * max(θ_l, FT(0))^3 # All soluble substrate, kgC m⁻³
    MM_sx = Sx / (kM_sx + Sx) # Availability of substrate factor, 0-1
    # Compute O2 from prognostic O2_a
    O2 = D_oa * O2_a * (max((ν - θ_l), 0)^(FT(4 / 3))) # Oxygen concentration
    MM_o2 = O2 / (kM_o2 + O2) # Oxygen limitation factor, 0-1
    R_sm = Vmax * MM_sx * MM_o2 # Respiration, kg C m⁻³ s⁻¹
    return R_sm
end


"""
    o2_concentration(O2_a::FT,
                     θ_a::FT,
                     T_soil::FT,
                     P_sfc::FT,
                     ) where {FT}

Computes the O₂ mass concentration in soil (kg/m³) from the volumetric fraction O2_a,
using the ideal gas law.

The O2 mass concentration in the soil air space is:
    ρ_O2 = θ_a * f_O2 * P * M_O2 / (R * T)

where:
- θ_a: volumetric air content (m³ air / m³ soil)
- f_O2 (O2_a): volumetric fraction of O2 in air (dimensionless, ~0.21)
- P: pressure (Pa)
- M_O2: molar mass of O2 (0.032 kg/mol)
- R: universal gas constant (8.314462618 J/(mol·K))
- T: temperature (K)
"""
function o2_concentration(
    O2_a::FT,
    θ_a::FT,
    T_soil::FT,
    P_sfc::FT,
) where {FT}
    R = FT(8.314462618)  # J/(mol·K)
    M_O2 = FT(0.032)      # kg/mol
    # O2 mass concentration in soil (kg/m³)
    ρ_O2 = θ_a * O2_a * P_sfc * M_O2 / (R * T_soil)
    return ρ_O2
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
    D =
        D0 *
        (FT(2)θ_a100^FT(3) + FT(0.04)θ_a100) *
        (θ_a / θ_a100)^(FT(2) + FT(3) / b)
    return D
end
