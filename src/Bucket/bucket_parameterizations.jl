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
    surface_albedo(albedo::BulkAlbedo{FT}, coords, S::FT, S_c::FT)::FT where {FT <: AbstractFloat}

Returns the bulk surface albedo, linearly interpolating between the albedo
of snow and of soil, based on the snow water equivalent S relative to
the parameter S_c.

The linear interpolation is taken from Lague et al 2019.
"""
function surface_albedo(
    albedo::BulkAlbedo{FT},
    coords,
    S::FT,
    S_c::FT,
)::FT where {FT <: AbstractFloat}
    (; α_snow, α_soil) = albedo
    α_soil_values = α_soil(coords)
    safe_S::FT = max(S, eps(FT))
    return (FT(1.0) - S / (S + S_c)) * α_soil_values + S / (S + S_c) * α_snow
end

"""
    infiltration_at_point(W::FT, M::FT, P::FT, E::FT, W_f::FT)::FT where {FT <: AbstractFloat}

Returns the soil infiltration given the current water content of the bucket W, 
the snow melt volume flux M, the precipitation volume flux P, the liquid evaporative volume
flux E, and the bucket capacity W_f. Positive values indicate increasing
soil moisture; the soil infiltration is the magnitude of the water
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
    safe_W = max(FT(0.0), W)
    if safe_W < FT(0.75) * W_f
        safe_W / (FT(0.75) * W_f)
    else
        FT(1.0)
    end
end

"""
    q_sat(T::FT, S::FT, ρ_sfc::FT, parameters::PE)::FT where {FT, PE}

Computes the saturation specific humidity for the land surface, over ice
if snow is present (S>0), and over water for a snow-free surface.
"""
function q_sat(T::FT, S::FT, ρ_sfc::FT, parameters::PE)::FT where {FT, PE}
    return (FT(1.0) - heaviside(S)) * Thermodynamics.q_vap_saturation_generic(
        parameters.earth_param_set,
        T,
        ρ_sfc,
        Thermodynamics.Liquid(),
    )
    +heaviside(S) * Thermodynamics.q_vap_saturation_generic(
        parameters.earth_param_set,
        T,
        ρ_sfc,
        Thermodynamics.Ice(),
    )
end
