export insolation, solar_flux_and_cos_sza

"""
    insolation(θ::FT, d::FT, param_set::IP.AIP) where {FT <: Real}

Returns the insolation given the zenith angle and earth-sun distance
param_set is an AbstractParameterSet from ClimaParams.jl.
"""
function insolation(θ::FT, d::FT, param_set::IP.AIP) where {FT <: Real}
    S0::FT = IP.tot_solar_irrad(param_set)
    d0::FT = IP.orbit_semimaj(param_set)
    # set max. zenith angle to π/2, insolation should not be negative
    if θ > FT(π) / 2
        θ = FT(π) / 2
    end
    # weighted irradiance (3.12)
    S = S0 * (d0 / d)^2
    # TOA insolation (3.15)
    F = S * cos(θ)
    return F
end

"""
    solar_flux_and_cos_sza(date::DateTime,
                      od::OrbitalData,
                      longitude::FT,
                      latitude::FT,
                      param_set::IP.AIP) where {FT <: Real}

Returns the top-of-atmosphere (TOA) solar flux, i.e. 
the total solar irradiance (TSI) weighted by the earth-sun distance
and cos(solar zenith angle) for input to RRTMGP.jl
param_set is an AbstractParameterSet from ClimaParams.jl.
"""
function solar_flux_and_cos_sza(
    date::DateTime,
    od::OrbitalData,
    longitude::FT,
    latitude::FT,
    param_set::IP.AIP,
) where {FT <: Real}
    S0::FT = IP.tot_solar_irrad(param_set)
    d0::FT = IP.orbit_semimaj(param_set)
    args = (
        Insolation.helper_instantaneous_zenith_angle(date, od, param_set)...,
        longitude,
        latitude,
    )
    # θ = solar zenith angle, ζ = solar azimuth angle, d = earth-sun distance
    θ, ζ, d = instantaneous_zenith_angle(args...)
    # set max. zenith angle to π/2, insolation should not be negative
    if θ > FT(π) / 2
        θ = FT(π) / 2
    end
    μ = cos(θ)
    # TOA solar flux (3.12)
    S = S0 * (d0 / d)^2
    # return inputs needed for RRTMGP.jl
    return S, μ
end
