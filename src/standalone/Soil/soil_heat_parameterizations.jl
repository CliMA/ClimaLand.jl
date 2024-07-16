export volumetric_internal_energy,
    temperature_from_ρe_int,
    volumetric_internal_energy_liq,
    volumetric_heat_capacity,
    κ_solid,
    κ_sat_unfrozen,
    κ_sat_frozen,
    thermal_conductivity,
    relative_saturation,
    κ_sat,
    kersten_number,
    κ_dry,
    thermal_time,
    phase_change_source

"""
    thermal_time(ρc::FT, Δz::FT, κ::FT) where {FT}

Returns the thermal timescale for temperature differences across
a typical thickness Δz to equilibrate.

Clip to 10x the timestep in order to respect CFL condition,
given that phase change is stepped explicitly.
"""
function thermal_time(ρc::FT, Δz::FT, κ::FT, Δt::FT) where {FT}
    return max(ρc * Δz^2 / κ, 10 * Δt)
end

"""
    phase_change_source(
        θ_l::FT,
        θ_i::FT,
        T::FT,
        τ::FT,
        ν::FT,
        θ_r::FT,
        hydrology_cm::C,
        earth_param_set::EP,
    ) where {FT, EP, C}
Returns the source term (1/s) used for converting liquid water
and ice into each other during phase changes. Note that
there are unitless prefactors multiplying this term in the 
equations.

Note that these equations match what is in Dall'Amico (for θstar,
ψ(T), ψw0). We should double check them in the case where we have
ϑ_l > θ_l, but they should be very close to the form we want regardless.
"""
function phase_change_source(
    θ_l::FT,
    θ_i::FT,
    T::FT,
    τ::FT,
    ν::FT,
    θ_r::FT,
    hydrology_cm::C,
    earth_param_set::EP,
) where {FT, EP, C}
    _ρ_i = FT(LP.ρ_cloud_ice(earth_param_set))
    _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
    _LH_f0 = FT(LP.LH_f0(earth_param_set))
    _T_freeze = FT(LP.T_freeze(earth_param_set))
    _grav = FT(LP.grav(earth_param_set))
    # According to Dall'Amico (text above equation 1), ψw0 corresponds
    # to the matric potential corresponding to the total water content (liquid and ice).
    θtot = min(_ρ_i / _ρ_l * θ_i + θ_l, ν)
    # This is consistent with Equation (22) of Dall'Amico
    ψw0 = matric_potential(hydrology_cm, effective_saturation(ν, θtot, θ_r))

    ψT = _LH_f0 / _T_freeze / _grav * (T - _T_freeze) * heaviside(_T_freeze - T)
    # Equation (23) of Dall'Amico
    θstar = inverse_matric_potential(hydrology_cm, ψw0 + ψT) * (ν - θ_r) + θ_r

    return (θ_l - θstar) / τ
end


"""
    volumetric_heat_capacity(
        θ_l::FT,
        θ_i::FT,
        ρc_ds::FT,
        earth_param_set::EP,
    ) where {FT,EP}

Compute the expression for volumetric heat capacity.
"""
function volumetric_heat_capacity(
    θ_l::FT,
    θ_i::FT,
    ρc_ds::FT,
    earth_param_set::EP,
) where {FT, EP}
    _ρ_i = FT(LP.ρ_cloud_ice(earth_param_set))
    ρcp_i = FT(LP.cp_i(earth_param_set) * _ρ_i)

    _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
    ρcp_l = FT(LP.cp_l(earth_param_set) * _ρ_l)

    ρc_s = ρc_ds + θ_l * ρcp_l + θ_i * ρcp_i
    return ρc_s
end
"""
    temperature_from_ρe_int(ρe_int::FT, θ_i::FT, ρc_s::FT
                            earth_param_set::EP) where {FT, EP}

A pointwise function for computing the temperature from the volumetric
internal energy, volumetric ice content, and volumetric heat capacity of
the soil.
"""
function temperature_from_ρe_int(
    ρe_int::FT,
    θ_i::FT,
    ρc_s::FT,
    earth_param_set::EP,
) where {FT, EP}

    _ρ_i = FT(LP.ρ_cloud_ice(earth_param_set))
    _T_ref = FT(LP.T_0(earth_param_set))
    _LH_f0 = FT(LP.LH_f0(earth_param_set))
    T = _T_ref + (ρe_int + θ_i * _ρ_i * _LH_f0) / ρc_s
    return T
end

"""
    volumetric_internal_energy(θ_i::FT, ρc_s::FT, T::FT,
                                 earth_param_set::EP) where {FT, EP}

A pointwise function for computing the volumetric internal energy of the soil,
given the volumetric ice content, volumetric heat capacity, and temperature.
"""
function volumetric_internal_energy(
    θ_i::FT,
    ρc_s::FT,
    T::FT,
    earth_param_set::EP,
) where {FT, EP}
    _ρ_i = FT(LP.ρ_cloud_ice(earth_param_set))
    _LH_f0 = FT(LP.LH_f0(earth_param_set))
    _T_ref = FT(LP.T_0(earth_param_set))
    ρe_int = ρc_s * (T - _T_ref) - θ_i * _ρ_i * _LH_f0
    return ρe_int
end

"""
    volumetric_internal_energy_liq(T::FT, earth_param_set::EP) where {FT, EP}


A pointwise function for computing the volumetric internal energy
of the liquid water in the soil, given the temperature T.
"""
function volumetric_internal_energy_liq(
    T::FT,
    earth_param_set::EP,
) where {FT, EP}

    _T_ref = FT(LP.T_0(earth_param_set))
    _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
    ρcp_l = FT(LP.cp_l(earth_param_set) * _ρ_l)
    ρe_int_l = ρcp_l * (T - _T_ref)
    return ρe_int_l
end

"""
    κ_sat(
        θ_l::FT,
        θ_i::FT,
        κ_sat_unfrozen::FT,
        κ_sat_frozen::FT
    ) where {FT}

Compute the expression for saturated thermal conductivity of soil matrix.
"""
function κ_sat(
    θ_l::FT,
    θ_i::FT,
    κ_sat_unfrozen::FT,
    κ_sat_frozen::FT,
) where {FT}
    θ_w = θ_l + θ_i
    if θ_w < eps(FT)
        return (κ_sat_unfrozen + κ_sat_frozen) / FT(2)
    else
        return FT(κ_sat_unfrozen^(θ_l / θ_w) * κ_sat_frozen^(θ_i / θ_w))
    end

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
            ν::FT
    ) where {FT}

Compute the expression for relative saturation. 
This is referred to as θ_sat in Balland and Arp's paper.
"""
function relative_saturation(θ_l::FT, θ_i::FT, ν::FT) where {FT}
    return max((θ_l + θ_i) / ν, eps(FT))
end

"""
    kersten_number(
        θ_i::FT,
        S_r::FT,
        α::FT,
        β::FT,
        ν_ss_om::FT,
        ν_ss_quartz::FT,
        ν_ss_gravel::FT,
        ) where {FT}

Compute the expression for the Kersten number, using the Balland
and Arp model.
"""
function kersten_number(
    θ_i::FT,
    S_r::FT,
    α::FT,
    β::FT,
    ν_ss_om::FT,
    ν_ss_quartz::FT,
    ν_ss_gravel::FT,
) where {FT}
    if θ_i < eps(FT)
        K_e =
            S_r^((FT(1) + ν_ss_om - α * ν_ss_quartz - ν_ss_gravel) / FT(2)) *
            (
                (FT(1) + exp(-β * S_r))^(-FT(3)) -
                ((FT(1) - S_r) / FT(2))^FT(3)
            )^(FT(1) - ν_ss_om)
    else
        K_e = S_r^(FT(1) + ν_ss_om)
    end
    return K_e
end



### Functions to be computed ahead of time, and values stored in
### EnergyHydrologyParameters
"""
    κ_solid(ν_ss_om::FT,
            ν_ss_quartz::FT,
            κ_om::FT,
            κ_quartz::FT,
            κ_minerals::FT) where {FT}

Computes the thermal conductivity of the solid material in soil.
The `_ss_` subscript denotes that the volumetric fractions of the soil
components are referred to the soil solid components, not including the pore
space.
"""
function κ_solid(
    ν_ss_om::FT,
    ν_ss_quartz::FT,
    κ_om::FT,
    κ_quartz::FT,
    κ_minerals::FT,
) where {FT}
    return κ_om^ν_ss_om *
           κ_quartz^ν_ss_quartz *
           κ_minerals^(FT(1) - ν_ss_om - ν_ss_quartz)
end


"""
    function κ_sat_frozen(
        κ_solid::FT,
        ν::FT,
        κ_ice::FT
    ) where {FT}
Computes the thermal conductivity for saturated frozen soil.
"""
function κ_sat_frozen(κ_solid::FT, ν::FT, κ_ice::FT) where {FT}
    return κ_solid^(FT(1.0) - ν) * κ_ice^(ν)
end

"""
    function κ_sat_unfrozen(
        κ_solid::FT,
        ν::FT,
        κ_l::FT
    ) where {FT}
Computes the thermal conductivity for saturated unfrozen soil.
"""
function κ_sat_unfrozen(κ_solid::FT, ν::FT, κ_l::FT) where {FT}
    return κ_solid^(FT(1.0) - ν) * κ_l^ν
end


"""
    function κ_dry(ρp::FT,
                   ν::FT,
                   κ_solid::FT,
                   κ_air::FT;
                   a = FT(0.053)) where {FT}

Computes the thermal conductivity of dry soil according
to the model of Balland and Arp.
"""
function κ_dry(ρp::FT, ν::FT, κ_solid::FT, κ_air::FT; a = FT(0.053)) where {FT}
    # compute the dry soil bulk density from particle density
    ρb_dry = (FT(1.0) - ν) * ρp
    numerator = (a * κ_solid - κ_air) * ρb_dry + κ_air * ρp
    denom = ρp - (FT(1.0) - a) * ρb_dry
    return numerator / denom
end
