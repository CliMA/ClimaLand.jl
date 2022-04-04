## To do: use clima parameters; remove hardwired parameters
export volumetric_internal_energy,
    temperature_from_ρe_int,
    volumetric_internal_energy_liq,
    volumetric_heat_capacity,
    k_solid,
    ksat_unfrozen,
    ksat_frozen,
    thermal_conductivity,
    relative_saturation,
    kersten_number

"""
    volumetric_heat_capacity(
        θ_l::FT,
        θ_i::FT,
        ρc_ds::FT,
        param_set::AbstractParameterSet
    ) where {FT}
Compute the expression for volumetric heat capacity.
"""
function volumetric_heat_capacity(
    θ_l::FT,
    θ_i::FT,
    ρc_ds::FT,
    param_set::AbstractParameterSet,
) where {FT}
    _ρ_i = FT(ρ_cloud_ice(param_set))
    ρcp_i = FT(cp_i(param_set) * _ρ_i)

    _ρ_l = FT(ρ_cloud_liq(param_set))
    ρcp_l = FT(cp_l(param_set) * _ρ_l)

    ρc_s = ρc_ds + θ_l * ρcp_l + θ_i * ρcp_i
    return ρc_s
end
"""
    temperature_from_ρe_int(ρe_int::FT, θ_i::FT, ρc_s::FT
                            param_set::AbstractParameterSet) where {FT}

A pointwise function for computing the temperature from the volumetric
internal energy, volumetric ice content, and volumetric heat capacity of
the soil.
"""
function temperature_from_ρe_int(ρe_int::FT, θ_i::FT, ρc_s::FT
                                 param_set::AbstractParameterSet) where {FT}
    _ρ_i = FT(ρ_cloud_ice(param_set))
    _T_ref = FT(T_0(param_set))
    _LH_f0 = FT(LH_f0(param_set))
    T = _T_ref + (ρe_int + θ_i * _ρ_i * _LH_f0) / ρc_s
    return T
end

"""
    volumetric_internal_energy(θ_i::FT, ρc_s::FT, T::FT,
                               param_set::AbstractParameterSet) where {FT}

A pointwise function for computing the volumetric internal energy of the soil,
given the volumetric ice content, volumetric heat capacity, and temperature.
"""
function volumetric_internal_energy(θ_i::FT, ρc_s::FT, T::FT
                                    param_set::AbstractParameterSet) where {FT}
    _ρ_i = FT(ρ_cloud_ice(param_set))
    _LH_f0 = FT(LH_f0(param_set))
    _T_ref = FT(T_0(param_set))
    ρe_int = ρc_s * (T - _T_ref) - θ_i * _ρ_i * _LH_f0
    return ρe_int
end

"""
    volumetric_internal_energy_liq(T::FT, param_set::AbstractParameterSet) where {FT}

A pointwise function for computing the volumetric internal energy
of the liquid water in the soil, given the temperature T.
"""
function volumetric_internal_energy_liq(T::FT, param_set::AbstractParameterSet) where {FT}
    _T_ref = FT(T_0(param_set))
    _ρ_l = FT(ρ_cloud_liq(param_set))
    ρcp_l = FT(cp_l(param_set) * _ρ_l)
    ρe_int_l = ρcp_l * (T - _T_ref)
    return ρe_int_l
end



"""
    function k_solid(
        soil_composition::PS,
        κ_quartz::FT,
        κ_minerals::FT,
        κ_om::FT,
    ) where {FT}
Computes the thermal conductivity of the solid material in soil.
The `_ss_` subscript denotes that the volumetric fractions of the soil
components are referred to the soil solid components, not including the pore
space.
"""
function k_solid(
    soil_composition::PS
    κ_quartz::FT,
    κ_minerals::FT,
    κ_om::FT,
    ) where {PS,FT}
    return κ_om^soil_composition.ν_ss_om *
           κ_quartz^soil_composition.ν_ss_quartz *
           κ_minerals^(FT(1) -
           soil_composition.ν_ss_om -
           soil_composition.ν_ss_quartz)
end


"""
    function ksat_frozen(
        κ_solid::FT,
        porosity::FT,
        κ_ice::FT
    ) where {FT}

Computes the thermal conductivity for saturated frozen soil.
"""
function ksat_frozen(κ_solid::FT, porosity::FT, κ_ice::FT) where {FT}
    return κ_solid^(FT(1.0) - porosity) * κ_ice^(porosity)
end

"""
    function ksat_unfrozen(
        κ_solid::FT,
        porosity::FT,
        κ_l::FT
    ) where {FT}

Computes the thermal conductivity for saturated unfrozen soil.
"""
function ksat_unfrozen(κ_solid::FT, porosity::FT, κ_l::FT) where {FT}
    return κ_solid^(FT(1.0) - porosity) * κ_l^porosity
end

"""
    saturated_thermal_conductivity(
        θ_l::FT,
        θ_i::FT,
        κ_sat_unfrozen::FT,
        κ_sat_frozen::FT
    ) where {FT}

Compute the expression for saturated thermal conductivity of soil matrix.
"""
function saturated_thermal_conductivity(
    θ_l::FT,
    θ_i::FT,
    κ_sat_unfrozen::FT,
    κ_sat_frozen::FT,
) where {FT}
    θ_w = θ_l + θ_i
    return FT(κ_sat_unfrozen^(θ_l / θ_w) * κ_sat_frozen^(θ_i / θ_w))
end

"""
    thermal_conductivity(
        κ_dry::FT,
        K_e::FT,
        κ_sat::FT
    ) where {FT}

Compute the expression for thermal conductivity of soil matrix.
"""
function thermal_conductivity(κ_dry::FT, K_e::FT, κ_sat::FT) where {FT}
    κ = K_e * κ_sat + (FT(1) - K_e) * κ_dry
    return κ
end


"""
    relative_saturation(
            θ_l::FT,
            θ_i::FT,
            porosity::FT
    ) where {FT}

Compute the expression for relative saturation. 

This is referred to as θ_sat in Balland and Arp's paper.
"""
function relative_saturation(θ_l::FT, θ_i::FT, porosity::FT) where {FT}
    return (θ_l + θ_i) / porosity
end

"""
    kersten_number(
        θ_i::FT,
        S_r::FT,
        soil_param_functions::PS
    ) where {FT, PS}

Compute the expression for the Kersten number.
"""
function kersten_number(
    θ_i::FT,
    S_r::FT,
    soil_composition::PS,
    soil_heat_params::PS2
) where {FT, PS,PS2}
    a = soil_heat_params.a
    b = soil_heat_params.b
    ν_ss_om = soil_composition.ν_ss_om
    ν_ss_quartz = soil_composition.ν_ss_quartz
    ν_ss_gravel = soil_composition.ν_ss_gravel

    if θ_i < eps(FT)
        K_e =
            S_r^((FT(1) + ν_ss_om - a * ν_ss_quartz - ν_ss_gravel) / FT(2)) *
            (
                (FT(1) + exp(-b * S_r))^(-FT(3)) -
                ((FT(1) - S_r) / FT(2))^FT(3)
            )^(FT(1) - ν_ss_om)
    else
        K_e = S_r^(FT(1) + ν_ss_om)
    end
    return K_e
end


"""
    function k_dry(
        param_set::AbstractParameterSet
        soil_composition::PS,
        soil_heat_params::PS2
    ) where {PS,PS2}

Computes the thermal conductivity of dry soil.
"""
function k_dry(
    param_set::AbstractParameterSet,
    soil_composition::PS,
    soil_heat_params::PS2
) where {PS,PS2}
    κ_dry_parameter = soil_heat_params.κ_dry_parameter
    FT = typeof(κ_dry_parameter)
    porosity = soil_composition.ν
    ρp = soil_composition.ρp
    κ_solid = soil_heat_params.κ_solid
    κ_air = FT(K_therm(param_set))
    ρb_val = ρb_ss(porosity, ρp)
    numerator = (κ_dry_parameter * κ_solid - κ_air) * ρb_val + κ_air * ρp
    denom = ρp - (FT(1.0) - κ_dry_parameter) * ρb_val
    return numerator / denom
end
