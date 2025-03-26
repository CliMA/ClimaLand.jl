## Surface characteristics functions
"""
    surface_albedo(model::BucketModel, Y, p)

Returns the bulk surface albedo, which gets updated in `update_aux`
via `next_albedo`.
"""
function ClimaLand.surface_albedo(model::BucketModel{FT}, Y, p) where {FT}
    return p.bucket.α_sfc
end

"""
   ClimaLand.surface_emissivity(model::BucketModel{FT}, Y, p)

Returns the emissivity for the bucket model (1.0).
"""
function ClimaLand.surface_emissivity(model::BucketModel{FT}, Y, p) where {FT}
    return FT(1)
end

"""
    ClimaLand.surface_temperature(model::BucketModel, Y, p)

a helper function which returns the surface temperature for the bucket
model, which is stored in the aux state.
"""
function ClimaLand.surface_temperature(model::BucketModel, Y, p, t)
    return p.bucket.T_sfc
end

"""
    ClimaLand.surface_specific_humidity(atmos, model::BucketModel, Y, p)

a helper function which returns the surface specific humidity for the bucket
model, which is stored in the aux state.
"""
function ClimaLand.surface_specific_humidity(
    atmos,
    model::BucketModel,
    Y,
    p,
    _...,
)
    return p.bucket.q_sfc
end

"""
    ClimaLand.surface_height(model::BucketModel, Y, p)

a helper function which returns the surface height for the bucket
model, which is zero currently.
"""
function ClimaLand.surface_height(model::BucketModel{FT}, Y, p) where {FT}
    return FT(0)
end

"""
    ClimaLand.surface_evaporative_scaling(model::BucketModel, Y, p)

a helper function which computes and returns the surface evaporative scaling
 factor for the bucket model.
"""
function ClimaLand.surface_evaporative_scaling(model::BucketModel, Y, p)
    beta =
        beta_factor.(
            Y.bucket.W,
            Y.bucket.σS,
            model.parameters.f_bucket * model.parameters.W_f,
            model.parameters.f_snow * model.parameters.σS_c,
            model.parameters.p,
        )
    return beta
end

"""
    beta_factor(W::FT, σS::FT, fW_f::FT, fσS_c::FT, p::FT) where {FT}

Computes the beta factor which scales the evaporation/sublimation from the potential
rate. The beta factor is given by:

β = (x/x_c)^p x < x_c
    1         otherwise

where x = W and x_c = f_bucket * W_f for the bucket,
and x = σS and x_c = f_snow *σS_c for snow.

"""
function beta_factor(W::FT, σS::FT, fW_f::FT, fσS_c::FT, p::FT) where {FT}
    snow_cover_fraction = heaviside(σS)
    return (1 - snow_cover_fraction) * β(W, fW_f, p) +
           snow_cover_fraction * β(σS, fσS_c, p)
end

"""
    partition_snow_surface_fluxes(
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
and a ground heat flux. Fluxes are given over snow covered area and multiplied by snow cover
fraciton elsewhere.

All fluxes are positive if they are in the direction from land towards
the atmosphere.
"""
function partition_snow_surface_fluxes(
    σS::FT,
    T_sfc::FT,
    τ::FT,
    snow_cover_fraction::FT,
    E::FT,
    F_sfc::FT,
    _ρLH_f0::FT,
    _T_freeze::FT,
) where {FT}
    σS = max(FT(0), σS) # clip to zero in case small and negative
    S = snow_cover_fraction > FT(0) ? σS / snow_cover_fraction : FT(0)
    F_available_to_melt = F_sfc + _ρLH_f0 * E # Eqn (23), negative as towards ground
    if S > -F_available_to_melt * τ / _ρLH_f0 # Eqn (24)
        F_melt = F_available_to_melt
    else
        F_melt = -S * _ρLH_f0 / τ
    end
    F_melt =
        F_melt * heaviside(T_sfc - _T_freeze) * heaviside(-F_available_to_melt)
    F_into_snow = -_ρLH_f0 * E + F_melt
    G_under_snow = (F_sfc - F_into_snow) # Eqn 22
    return (;
        F_melt = F_melt,
        F_into_snow = F_into_snow,
        G_under_snow = G_under_snow,
    )
end


"""
    infiltration_at_point(W::FT, M::FT, P::FT, E::FT, W_f::FT)::FT where {FT <: AbstractFloat}

Returns the infiltration given the current water content of the bucket W,
the snow melt volume flux M, the precipitation volume flux P, the liquid evaporative volume
flux E, and the bucket capacity W_f.

Extra inflow when the bucket is at capacity runs off.
Note that all fluxes are positive if towards the atmosphere.
"""
function infiltration_at_point(
    W::FT,
    M::FT,
    P::FT,
    E::FT,
    W_f::FT,
)::FT where {FT <: AbstractFloat}
    f = FT(P + M + E)
    if (W > W_f) & (f < 0) # Equation (2) of text
        return FT(0) # we are at capacity and the net flux is negative
    else
        return f # can be positive (water loss to atmos) or negative
    end
end

"""
    β(x::FT, x_c::FT, p::FT) where {FT}

Returns the coefficient which scales the evaporation or
sublimation rate based on the bucket water
or snow levels.

Over ground, x_c is a default of 75% of the capacity, since the ground
evaporation rate remains near the potential rate until water has dropped
sufficiently.

Over snow, x_c taken a default of 10% of the value around which snow starts to become patchy,
since snow sublimates at the potential rate in general. We use the β function
mainly to damp sublimation to zero for vanishing snowpacks. Note that the snow cover fraction
returns zero for 0 < σS < eps(FT) while this function returns a nonzero function.
"""
function β(x::FT, x_c::FT, p::FT) where {FT}
    safe_x = max(0, x)
    if safe_x < x_c
        (safe_x / x_c)^p
    else
        FT(1.0)
    end
end

"""
    saturation_specific_humidity(T::FT, ρ_sfc::FT, thermo_parameters::TPE)::FT where {FT, TPE}

Computes the saturation specific humidity for the land surface, over ice
if the temperature is below freezing, and over water otherwise.
"""
function saturation_specific_humidity(
    T::FT,
    ρ_sfc::FT,
    thermo_params::TPE,
)::FT where {FT, TPE}
    if T > TP.T_freeze(thermo_params)
        Thermodynamics.q_vap_saturation_generic(
            thermo_params,
            T,
            ρ_sfc,
            Thermodynamics.Liquid(),
        )
    else
        Thermodynamics.q_vap_saturation_generic(
            thermo_params,
            T,
            ρ_sfc,
            Thermodynamics.Ice(),
        )
    end
end
