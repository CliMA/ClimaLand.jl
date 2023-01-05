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
    microbe_source(T_soil::FT,
                   θ_l::FT,
                   Csom::FT,
                   params::SoilCO2ModelParameters{FT}
                   ) where {FT}

Computes the CO₂ production in the soil by microbes, in depth and time (kg C / m^3/s), using
the Dual Arrhenius Michaelis Menten model (Davidson et al., 2012).
"""
function microbe_source(
    T_soil::FT,
    θ_l::FT,
    Csom::FT,
    params::SoilCO2ModelParameters{FT},
) where {FT}
    @unpack α_sx,
    Ea_sx,
    kM_sx,
    kM_o2,
    ν,
    D_liq,
    p_sx,
    D_oa,
    O2_a,
    earth_param_set = params
    R = FT(LSMP.gas_constant(earth_param_set))
    Vmax = α_sx * exp(-Ea_sx / (R * T_soil)) # Maximum potential rate of respiration
    Sx = p_sx * Csom * D_liq * θ_l^3 # All soluble substrate, kgC m⁻³
    MM_sx = Sx / (kM_sx + Sx) # Availability of substrate factor, 0-1 
    O2 = D_oa * O2_a * ((ν - θ_l)^(FT(4 / 3))) # Oxygen concentration
    MM_o2 = O2 / (kM_o2 + O2) # Oxygen limitation factor, 0-1
    R_sm = Vmax * MM_sx * MM_o2 # Respiration, kg C m⁻³ s⁻¹
    return R_sm
end


"""
    volumetric_air_content(θ_w::FT,
                           params::SoilCO2ModelParameters{FT}
                           ) where {FT}
    
Computes the volumetric air content (`θ_a`) in the soil, 
which is related to the total soil porosity (`ν`) and 
volumetric soil water content (`θ_w = θ_l+θ_i`).
"""
function volumetric_air_content(
    θ_w::FT,
    params::SoilCO2ModelParameters{FT},
) where {FT}
    θ_a = params.ν - θ_w
    return θ_a
end

"""
    co2_diffusivity(
                    T_soil::FT,
                    θ_w::FT,
                    params::SoilCO2ModelParameters{FT}
                    ) where {FT}

Computes the diffusivity of CO₂ within the soil (D).

First, D0 is computed using the temperature within the soil (`T_soil` in K) and 
pressure at the surface of the soil (`P_sfc` in Pa), using reference
values of `T_ref` and `P_ref` (273 K and 101325 Pa). Here, `θ_a` is the 
volumetric air content and `θ_a100` is the volumetric air content 
at a soil water potential of 
100cm, and b is the pore size distribution of the soil.
"""
function co2_diffusivity(
    T_soil::FT,
    θ_w::FT,
    params::SoilCO2ModelParameters{FT},
) where {FT}
    @unpack P_sfc, D_ref, θ_a100, b, ν, earth_param_set = params
    T_ref = FT(LSMP.T_0(earth_param_set))
    P_ref = FT(LSMP.P_ref(earth_param_set))
    θ_a = volumetric_air_content(θ_w, params)
    D0 = D_ref * (T_soil / T_ref)^FT(1.75) * (P_ref / P_sfc)
    D =
        D0 *
        (FT(2)θ_a100^FT(3) + FT(0.04)θ_a100) *
        (θ_a / θ_a100)^(FT(2) + FT(3) / b)
    return D
end
