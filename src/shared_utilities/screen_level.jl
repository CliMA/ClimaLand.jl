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
    screen_level_values(T_sfc, q_sfc, ustar, ζ, Δz_eff, z0m, z0h, T_atmos, q_atmos,
                        z_screen, z_anemometer, earth_param_set)

Return the NamedTuple `(; T, q, u, g_h)` of the air temperature [K] and specific
humidity [kg/kg] at the height `z_screen` [m] above the apparent sink for heat
`displ + z0h`, the wind speed [m/s] at the height `z_anemometer` [m] above the
apparent sink for momentum `displ + z0m`, and the heat conductance `g_h` [m/s]
of a surface, reconstructed from the Monin-Obukhov profiles of its flux solve
(`turbulent_fluxes_at_a_point`) with `SurfaceFluxes.screen_level_values`. The
inputs are the surface temperature `T_sfc` [K] and specific humidity `q_sfc` at
which the fluxes were evaluated, the friction velocity `ustar` [m/s], the
stability parameter `ζ` and the effective forcing height `Δz_eff` [m], the
height of the forcing above the displacement height, of the solve, the roughness
lengths for momentum `z0m` and heat `z0h` [m], and the atmospheric temperature
`T_atmos` [K] and specific humidity `q_atmos` at the forcing height.

Between the surface and the forcing height, a quantity `X` carried by the heat
profile takes the value `X_sfc + (X_atmos - X_sfc) F̂_h(z) / F̂_h(Δz_eff)` at the
effective height `z`, with the dimensionless profile `F̂_h` of
`SurfaceFluxes.dimensionless_profile_value` at the stability of the solve. The
temperature follows this relation in terms of the dry static energy, with the
surface state at the displacement height (`SurfaceFluxes.surface_geopotential`),
so it includes the adiabatic temperature change `g / c_p` per meter between the
screen and forcing heights, and the displacement height itself does not enter.
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
    z0m::FT,
    z0h::FT,
    T_atmos::FT,
    q_atmos::FT,
    z_screen::FT,
    z_anemometer::FT,
    earth_param_set,
) where {FT}
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    κ = SurfaceFluxes.Parameters.von_karman_const(surface_flux_params)
    L_eff = Δz_eff / ζ # Inf when neutral
    rsl = SurfaceFluxes.NoRoughnessSubLayer()
    F̂_h_ref = SurfaceFluxes.dimensionless_profile_value(
        surface_flux_params,
        L_eff,
        z0h,
        Δz_eff,
        Δz_eff,
        UF.HeatTransport(),
        UF.PointValueScheme(),
        rsl,
    )
    g_h = ifelse(F̂_h_ref > 0, κ * ustar / max(F̂_h_ref, eps(FT)), FT(0))
    sc = SurfaceFluxes.SurfaceFluxConditions{FT}(
        FT(0),
        FT(0),
        FT(0),
        FT(0),
        FT(0),
        ustar,
        ζ,
        FT(0),
        g_h,
        T_sfc,
        q_sfc,
        L_eff,
        L_eff,
        ζ,
        true,
    )
    inputs = (;
        T_int = T_atmos,
        q_tot_int = q_atmos,
        q_liq_int = FT(0),
        q_ice_int = FT(0),
        Δz = Δz_eff,
        d = FT(0),
        roughness_model = SurfaceFluxes.ConstantRoughnessParams{FT}(z0m, z0h),
        roughness_inputs = nothing,
        rsl_model = rsl,
    )
    (; T, q, u) = SurfaceFluxes.screen_level_values(
        surface_flux_params,
        sc,
        inputs,
        z_screen,
        z_anemometer,
    )
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
