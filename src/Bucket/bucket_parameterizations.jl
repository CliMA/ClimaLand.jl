"""
     heaviside(x::FT)::FT where {FT}

Computes the heaviside function.
"""
function heaviside(x::FT)::FT where {FT}
    if x > eps(FT)
        return FT(1.0)
    else
        return FT(0.0)
    end
end

"""
    beta_factor(W::FT, σS::FT, W_f::FT) where {FT}
Computes the beta factor which scales the evaporation from the potential
rate.
"""
function beta_factor(W::FT, σS::FT, W_f::FT) where {FT}
    snow_cover_fraction = heaviside(σS)
    return (1 - snow_cover_fraction) * β(W, W_f) + snow_cover_fraction * 1
end

"""
    partition_surface_fluxes(
        σS::FT,
        T_sfc::FT,
        τ::FT,
        snow_cover_fraction::FT,
        E::FT,
        F_sfc::FT,
        _LH_f0::FT,
        _T_freeze::FT,
    ) where{FT}

Partitions the surface fluxes in a flux for melting snow, a flux for sublimating snow,
and a ground heat flux. 

All fluxes are positive if they are in the direction from land towards
the atmosphere.
"""
function partition_surface_fluxes(
    σS::FT,
    T_sfc::FT,
    τ::FT,
    snow_cover_fraction::FT,
    E::FT,
    F_sfc::FT,
    _ρLH_f0::FT,
    _T_freeze::FT,
) where {FT}
    F_available_to_melt = F_sfc + _ρLH_f0 * E # Eqn (20), negative as towards ground
    if σS > -F_available_to_melt * τ / _ρLH_f0 # Eqn (21)
        F_melt = F_available_to_melt
    else
        F_melt = -σS * _ρLH_f0 / τ
    end
    F_melt =
        F_melt * heaviside(T_sfc - _T_freeze) * heaviside(-F_available_to_melt)
    F_into_snow = -_ρLH_f0 * E * snow_cover_fraction + F_melt # F_melt is already multiplied by σ
    G = (F_sfc - F_into_snow) # Eqn 20
    return (; F_melt = F_melt, F_into_snow = F_into_snow, G = G)
end

"""
    surface_albedo(α_sfc::FT, α_snow::FT, σS::FT, σS_c::FT)::FT where {FT <: AbstractFloat}

Returns the bulk surface albedo, linearly interpolating between the albedo
of snow and of the surface, based on the snow water equivalent S relative to
the parameter S_c.

The linear interpolation is taken from Lague et al 2019.
"""
function surface_albedo(
    α_sfc::FT,
    α_snow::FT,
    σS::FT,
    σS_c::FT,
)::FT where {FT <: AbstractFloat}
    safe_σS::FT = max(σS, eps(FT))
    return (1 - σS / (σS + σS_c)) * α_sfc + σS / (σS + σS_c) * α_snow
end

"""
    infiltration_at_point(W::FT, M::FT, P::FT, E::FT, W_f::FT)::FT where {FT <: AbstractFloat}

Returns the infiltration given the current water content of the bucket W, 
the snow melt volume flux M, the precipitation volume flux P, the liquid evaporative volume
flux E, and the bucket capacity W_f. Positive values indicate increasing
soil moisture; the infiltration is the magnitude of the water
flux into the soil.

Extra inflow when the bucket is at capacity runs off.
Note that P and M are positive by definition, while E can be 
positive (evaporation) or negative (condensation).
"""
function infiltration_at_point(
    W::FT,
    M::FT,
    P::FT,
    E::FT,
    W_f::FT,
)::FT where {FT <: AbstractFloat}
    if (W > W_f) & (P + M - E > 0) # Equation (4b) of text
        return FT(0)
    else
        return FT(P + M - E)
    end
end

"""
    β(W::FT, W_f::FT) where {FT}

Returns the coefficient which scales the saturated
specific humidity at the surface based on the bucket water
levels, which is then used
 to obtain the
true specific humidity of the soil surface <= q_sat.
"""
function β(W::FT, W_f::FT) where {FT}
    safe_W = max(0, W)
    if safe_W < FT(0.75) * W_f
        safe_W / (FT(0.75) * W_f)
    else
        FT(1.0)
    end
end

"""
    saturation_specific_humidity(T::FT, σS::FT, ρ_sfc::FT, parameters::PE)::FT where {FT, PE}

Computes the saturation specific humidity for the land surface, over ice
if snow is present (σS>0), and over water for a snow-free surface.
"""
function saturation_specific_humidity(
    T::FT,
    σS::FT,
    ρ_sfc::FT,
    parameters::PE,
)::FT where {FT, PE}
    thermo_params = LSMP.thermodynamic_parameters(parameters.earth_param_set)
    return (1 - heaviside(σS)) * Thermodynamics.q_vap_saturation_generic(
        thermo_params,
        T,
        ρ_sfc,
        Thermodynamics.Liquid(),
    ) +
           heaviside(σS) * Thermodynamics.q_vap_saturation_generic(
        thermo_params,
        T,
        ρ_sfc,
        Thermodynamics.Ice(),
    )
end
