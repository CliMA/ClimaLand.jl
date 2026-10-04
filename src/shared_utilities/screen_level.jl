"""
Reconstruction of near-surface ("screen-level") air temperature, specific
humidity, and wind speed from the Monin-Obukhov profiles of the turbulent flux
solves of the land surfaces. The land models store the inputs of this
reconstruction with their turbulent fluxes (see
`turbulent_fluxes_from_output`), and the diagnostics `t2m`, `q2m`, and `u10m`
evaluate it.
"""

import SurfaceFluxes.UniversalFunctions as UF

# WMO standard heights of screen-level temperature and humidity and of the
# anemometer, above the apparent sinks for heat and momentum
screen_height(::Type{FT}) where {FT} = FT(2)
anemometer_height(::Type{FT}) where {FT} = FT(10)

"""
    profile_shape(z, Δz_eff, ζ, z0, transport, surface_flux_params)

Return the dimensionless Monin-Obukhov profile `F̂(z) = log(z / z0) - ψ(z / L) +
ψ(z0 / L)` of SurfaceFluxes.jl for momentum or heat (`transport`,
`UF.MomentumTransport()` or `UF.HeatTransport()`) at the effective height `z`
above the displacement height, with the roughness length `z0` and the Obukhov
length `L = Δz_eff / ζ` of a flux solve at the effective forcing height
`Δz_eff` and stability parameter `ζ`. The height is clamped to the range from
`z0`, where the profile is zero, to the forcing height, so that the profile
does not extrapolate beyond the forcing.
"""
function profile_shape(z, Δz_eff, ζ, z0, transport, surface_flux_params)
    FT = typeof(z)
    κ = SurfaceFluxes.Parameters.von_karman_const(surface_flux_params)
    L = Δz_eff / ζ # Inf when neutral
    # With scale κ and zero surface value, the profile value is F̂ itself
    return SurfaceFluxes.compute_profile_value(
        surface_flux_params,
        L,
        z0,
        max(min(z, Δz_eff), z0),
        κ,
        FT(0),
        transport,
    )
end

"""
    screen_level_values(T_sfc, q_sfc, ustar, ζ, Δz_eff, displ, z0m, z0h, T_atmos,
                        q_atmos, z_screen, z_anemometer, earth_param_set)

Return the NamedTuple `(; T, q, u, g_h)` of the air temperature [K] and specific
humidity [kg/kg] at the height `z_screen` [m] above the apparent sink for heat
`displ + z0h`, the wind speed [m/s] at the height `z_anemometer` [m] above the
apparent sink for momentum `displ + z0m`, and the heat conductance `g_h` [m/s]
of a surface, reconstructed from the Monin-Obukhov profiles of its flux solve
(`turbulent_fluxes_at_a_point`). The inputs are the surface temperature
`T_sfc` [K] and specific humidity `q_sfc` at which the fluxes were evaluated,
the friction velocity `ustar` [m/s], the stability parameter `ζ` and the
effective forcing height `Δz_eff` [m], the height of the forcing above the
displacement height `displ` [m], of the solve, the roughness lengths for
momentum `z0m` and heat `z0h` [m], and the atmospheric temperature `T_atmos`
[K] and specific humidity `q_atmos` at the forcing height.

Between the surface and the forcing height, a quantity `X` carried by the heat
profile takes the value `X_sfc + (X_atmos - X_sfc) F̂_h(z) / F̂_h(Δz_eff)` at the
effective height `z`, with the dimensionless profile `F̂_h` of
[`profile_shape`](@ref) at the stability of the solve. The temperature follows
this relation in terms of the dry static energy, so it includes the adiabatic
temperature change `g / c_p` per meter between the screen and forcing heights.
The wind speed is `ustar F̂_m(z) / κ`, which at the forcing height is the wind
speed the solve used (including gustiness). Levels at or above the forcing
height take the forcing values.
"""
function screen_level_values(
    T_sfc::FT,
    q_sfc::FT,
    ustar::FT,
    ζ::FT,
    Δz_eff::FT,
    displ::FT,
    z0m::FT,
    z0h::FT,
    T_atmos::FT,
    q_atmos::FT,
    z_screen::FT,
    z_anemometer::FT,
    earth_param_set,
) where {FT}
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    κ = SurfaceFluxes.Parameters.von_karman_const(surface_flux_params)
    _grav = LP.grav(earth_param_set)
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    heat = UF.HeatTransport()
    momentum = UF.MomentumTransport()

    F̂_h_ref = profile_shape(Δz_eff, Δz_eff, ζ, z0h, heat, surface_flux_params)
    z_T = z0h + z_screen
    F̂_h = profile_shape(z_T, Δz_eff, ζ, z0h, heat, surface_flux_params)
    r = F̂_h_ref > 0 ? min(F̂_h / F̂_h_ref, FT(1)) : FT(1)
    # Heights above the surface of the forcing and of the screen level
    Δz_atmos = Δz_eff + displ
    Δz_T = min(z_T, Δz_eff) + displ
    T = T_sfc + (T_atmos - T_sfc) * r + _grav / cp_d * (r * Δz_atmos - Δz_T)
    q = q_sfc + (q_atmos - q_sfc) * r

    z_u = z0m + z_anemometer
    F̂_m = profile_shape(z_u, Δz_eff, ζ, z0m, momentum, surface_flux_params)
    u = ustar * max(F̂_m, FT(0)) / κ
    g_h = F̂_h_ref > 0 ? κ * ustar / F̂_h_ref : FT(0)
    return (; T, q, u, g_h)
end

"""
    screen_level_mean(::Val{name}, (w₁, s₁), (w₂, s₂), ...)

Return the mean of the screen-level value `name` (`:T`, `:q`, or `:u`) of the
NamedTuples `sᵢ` of [`screen_level_values`](@ref) of the surfaces of a land
model, weighted by the area fractions `wᵢ` of the surfaces times their heat
conductances `sᵢ.g_h`. The weighting treats the surfaces as parallel sources
to the atmosphere, as the land models do: the mean is the value a thermometer
would read if each surface's share of the exchange with the atmosphere is its
share of the total conductance. If the total conductance is zero, the value of
the first surface is returned.
"""
@inline function screen_level_mean(::Val{name}, pairs...) where {name}
    total = unrolled_sum(pairs) do (w, s)
        w * s.g_h
    end
    weighted = unrolled_sum(pairs) do (w, s)
        w * s.g_h * getproperty(s, name)
    end
    first_value = getproperty(pairs[1][2], name)
    return ifelse(total > 0, weighted / total, first_value)
end

# Sum over a tuple without allocating or looping at runtime
@inline unrolled_sum(f, t::Tuple{}) = false
@inline unrolled_sum(f, t::Tuple) = f(first(t)) + unrolled_sum(f, Base.tail(t))
