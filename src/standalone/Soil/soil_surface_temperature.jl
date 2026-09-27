#=
Soil skin (surface) temperature.

The soil surface ("skin") is treated as a layer with zero heat capacity that
exchanges radiation and turbulent heat/vapor with the air above, and is
connected to the center of the top soil layer by a thermal resistance

    r = Δz_top / κ_top + r_litter  (m² K / W),

where Δz_top is the distance between the surface and the center of the top
layer, κ_top is the thermal conductivity of the top layer, and r_litter is the
resistance of the litter, thatch, and standing dead material between the skin
and the mineral soil, supplied by the canopy biomass model
(`canopy.biomass.r_litter`) and zero for bare soil. The skin temperature T_sfc
satisfies the surface energy balance

    f(T_sfc) = SW_n + LW_n(T_sfc) + L(T_sfc) + H(T_sfc) - (T_top - T_sfc)/r = 0,

with all fluxes positive upward. It is solved for within the Monin-Obukhov
iterations of SurfaceFluxes.jl, in the same way as the snow surface temperature
(see `Snow.update_surf_temp!`), and the same solve yields the turbulent fluxes
at T_sfc, which are stored with it in `p.soil.turbulent_fluxes`. The surface
humidity and upwelling longwave radiation of the soil are evaluated at T_sfc
via `component_temperature`, and the soil receives the heat flux
(T_top - T_sfc)/r.
When the top cell is frozen (contains ice and is below the depressed freezing
temperature), T_sfc is capped at the depressed freezing temperature; the soil
then receives the atmospheric fluxes at that temperature, which exceed the
skin-top conduction, and the excess melts ice in the top cell.
With a `SlabLitter` surface layer, the node below the skin is the litter slab
instead of the top soil cell (see `litter_layer.jl`). Conduction between the
soil and a snowpack or lake sediment is unaffected and uses the top cell
temperature.

Under a canopy, the turbulent exchange between the ground and the atmosphere
is reduced. The conductance used for the soil turbulent fluxes is

    g_eff = W g_h + (1 - W) / (1/g_h + r_under),

where g_h is the conductance of bare soil computed by SurfaceFluxes.jl, W is
the canopy gap fraction (`p.soil.W_gap`), and r_under is the resistance between
the ground and the canopy air in a dense canopy (`p.soil.r_undercanopy`), in
series with the aerodynamic resistance to the atmosphere (see
`soil_conductance_ratio`). W and r_under are set by the integrated land models
with a canopy (see `src/integrated/undercanopy_conductance.jl`); without a
canopy, W = 1, r_under = 0, and g_eff = g_h. Because SurfaceFluxes.jl computes
the fluxes with g_h, the soil temperature and humidity passed to it are those
at an "interface" such that the fluxes with conductance g_h from the interface
equal the fluxes with conductance g_eff from the skin:

    T_i = T_a + (g_eff / g_h) (T_sfc - T_a),

where T_a = T_int + (Φ_int - Φ_sfc) / cp_d is the atmospheric temperature
brought dry-adiabatically to the surface, and analogously for the humidity (see
`soil_surface_vapor_weight`). This is the same approach as used for the canopy
(see `get_update_surface_temperature_function(::CanopyModel, Y, p)`).
=#

"""
    soil_surface_temperature(bc, p)

Return the soil surface temperature [K]: the skin temperature
`p.soil.turbulent_fluxes.T_sfc` for atmospherically driven soil, and the top
layer temperature otherwise.

Called from `component_temperature(::EnergyHydrology, Y, p)`.
"""
soil_surface_temperature(::AtmosDrivenFluxBC, p) = p.soil.turbulent_fluxes.T_sfc
soil_surface_temperature(_, p) =
    ClimaLand.Domains.top_center_to_surface(p.soil.T)

"""
    soil_surface_thermal_resistance(Δz_top, κ_top, r_litter)

Return the thermal resistance [m² K/W] between the soil skin and the center of
the top soil layer: the half-cell conduction resistance `Δz_top/κ_top` plus the
litter resistance `r_litter`.
"""
soil_surface_thermal_resistance(Δz_top, κ_top, r_litter) =
    Δz_top / κ_top + r_litter

"""
    soil_conductance_ratio(W::FT, r_under::FT, g_h::FT) where {FT}

Return the ratio `f_eff = g_eff / g_h` [-] of the conductance used for the
turbulent fluxes of the soil to the conductance of bare soil `g_h` [m/s], with

    g_eff = W g_h + (1 - W) / (1/g_h + r_under).

Here `W` is the canopy gap fraction [-] and `r_under` the resistance [s/m]
between the ground and the canopy air in a dense canopy, which is in series
with the aerodynamic resistance between the canopy air and the atmosphere. The
latter is approximated by the bare-soil resistance `1/g_h` computed by
SurfaceFluxes.jl, so that it includes the effects of atmospheric stability.
Then `W ≤ f_eff ≤ 1`, and `f_eff = 1` for bare soil (`W = 1` or `r_under = 0`).
"""
soil_conductance_ratio(W::FT, r_under::FT, g_h::FT) where {FT} =
    W + (1 - W) / (1 + g_h * r_under)

"""
    adiabatic_surface_air_temperature(inputs, param_set, thermo_params)

Return the temperature [K] of the air at the atmospheric reference height
brought dry-adiabatically to the surface, `T_int + (Φ_int - Φ_sfc) / cp_d`. The
sensible heat flux computed by SurfaceFluxes.jl is proportional to the
difference of the surface temperature from this value.
"""
function adiabatic_surface_air_temperature(inputs, param_set, thermo_params)
    Φ_sfc = SurfaceFluxes.surface_geopotential(inputs)
    Φ_int = SurfaceFluxes.interior_geopotential(param_set, inputs)
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    return inputs.T_int + (Φ_int - Φ_sfc) / cp_d
end

"""
    soil_interface_temperature(T_skin, T_a, f_eff)

Return the temperature [K] passed to SurfaceFluxes.jl, `T_a + f_eff (T_skin -
T_a)`, such that the sensible heat flux with the bare soil conductance `g_h`
equals that from the skin with conductance `g_eff = f_eff g_h`.
"""
soil_interface_temperature(T_skin, T_a, f_eff) = T_a + f_eff * (T_skin - T_a)

"""
    soil_skin_temperature(T_i, T_a, f_eff)

Return the skin temperature [K] with interface temperature `T_i`; the inverse
of [`soil_interface_temperature`](@ref).
"""
soil_skin_temperature(T_i::FT, T_a::FT, f_eff::FT) where {FT} =
    T_a + (T_i - T_a) / max(f_eff, eps(FT))

"""
    soil_surface_vapor_weight(q_air, q_src, g_liq, g_h, f_eff, β_ice, frozen)

Return the weight `w` [-] such that the soil surface specific humidity passed
to SurfaceFluxes.jl is `q_sfc = w * q_src + (1 - w) * q_air`, where `q_src` is
the specific humidity of the vapor source (see [`soil_vapor_source`](@ref)).
The vapor flux is then `ρ g_h (q_sfc - q_air) = ρ w g_h (q_src - q_air)`:
- For evaporation from unfrozen soil, `w g_h` is the conductance of the dry
  soil layer `g_liq` in series with the turbulent conductance
  `g_eff = f_eff g_h`, so `w = f_eff g_liq / (f_eff g_h + g_liq)`.
- For sublimation from frozen soil, `w = β_ice f_eff`.
- For dew or frost formation (`q_air > q_src`), which occurs at the surface
  and is not limited by the dry soil layer, `w = f_eff` (as in CLM5).
The flux vanishes at `q_air = q_src`, so it is continuous in `q_air`.

This is the same parameterization as
`get_update_surface_humidity_function(::EnergyHydrology)`, which the turbulent
flux solve uses.
"""
function soil_surface_vapor_weight(
    q_air::FT,
    q_src::FT,
    g_liq::FT,
    g_h::FT,
    f_eff::FT,
    β_ice::FT,
    frozen::Bool,
) where {FT}
    if q_air > q_src
        return f_eff
    elseif frozen
        return β_ice * f_eff
    else
        return f_eff * g_liq / (f_eff * g_h + g_liq)
    end
end

"""
    soil_skin_state(T_sfc, T_atmos, P_atmos, q_atmos, Δz, param_set, thermo_params,
                    ψ_sfc, Tf_depressed, earth_param_set) -> NamedTuple
    soil_skin_state(T_sfc, inputs, thermo_params, param_set, ψ_sfc, Tf_depressed,
                    earth_param_set) -> NamedTuple

Return `(; ρ_sfc, qsat, ∂qsat∂T, frozen, hr)`: the surface air density
[kg/m³], the saturation specific humidity at `T_sfc` adjusted for the soil
water potential [kg/kg], its derivative with respect to temperature [1/K],
whether the surface is frozen, and the Kelvin factor `hr` [-] (one over frozen
soil). The atmospheric state at the reference height `Δz` above the surface is
given either explicitly or through the SurfaceFluxes `inputs`.

Called from [`soil_skin_fluxes`](@ref), [`update_soil_q_vap_sfc_scheme`](@ref),
and [`solve_soil_surface_temperature_at_a_point`](@ref).
"""
function soil_skin_state(
    T_sfc::FT,
    T_atmos::FT,
    P_atmos::FT,
    q_atmos::FT,
    Δz::FT,
    param_set,
    thermo_params,
    ψ_sfc::FT,
    Tf_depressed::FT,
    earth_param_set,
) where {FT}
    ρ_sfc =
        ClimaLand.compute_ρ_sfc(param_set, T_atmos, P_atmos, q_atmos, Δz, T_sfc)
    qsat = soil_specific_humidity(
        T_sfc,
        ρ_sfc,
        ψ_sfc,
        Tf_depressed,
        earth_param_set,
    )
    # At T_sfc = Tf_depressed (the cap below), qsat is over ice, so the surface
    # is frozen there as well
    frozen = T_sfc <= Tf_depressed
    # qsat already includes the Kelvin factor for unfrozen soil
    hr = frozen ? FT(1) : soil_kelvin_factor(T_sfc, ψ_sfc, earth_param_set)
    LH =
        frozen ? Thermodynamics.latent_heat_sublim(thermo_params, T_sfc) :
        Thermodynamics.latent_heat_vapor(thermo_params, T_sfc)
    # The weak temperature dependence of the water potential factor is neglected
    ∂qsat∂T =
        Thermodynamics.∂q_vap_sat_∂T_from_L(thermo_params, qsat, LH, T_sfc)
    return (; ρ_sfc, qsat, ∂qsat∂T, frozen, hr)
end

function soil_skin_state(
    T_sfc::FT,
    inputs,
    thermo_params,
    param_set,
    ψ_sfc::FT,
    Tf_depressed::FT,
    earth_param_set,
) where {FT}
    P_atmos = Thermodynamics.air_pressure(
        thermo_params,
        inputs.T_int,
        inputs.ρ_int,
        inputs.q_tot_int,
    )
    return soil_skin_state(
        T_sfc,
        inputs.T_int,
        P_atmos,
        inputs.q_tot_int,
        inputs.Δz,
        param_set,
        thermo_params,
        ψ_sfc,
        Tf_depressed,
        earth_param_set,
    )
end

"""
    soil_vapor_source(q_air, qsat, ∂qsat∂T, hr, frozen)

Return the specific humidity [kg/kg] of the vapor source at the soil surface
and its derivative with respect to the skin temperature [1/K], consistent with
`effective_soil_vapor_humidity`: the soil-water-adjusted saturation value
`qsat` when the air is drier than that, `q_air` when the air is between `qsat`
and saturation over free water (so the Kelvin effect alone draws no vapor into
the soil), and saturation over free water when dew forms. For frozen soil, the
source is saturated with respect to ice.
"""
function soil_vapor_source(
    q_air::FT,
    qsat::FT,
    ∂qsat∂T::FT,
    hr::FT,
    frozen::Bool,
) where {FT}
    frozen && return (qsat, ∂qsat∂T)
    q_src = effective_soil_vapor_humidity(q_air, qsat, hr)
    ∂q_src∂T =
        q_air < qsat ? ∂qsat∂T :
        (q_air > qsat / max(hr, eps(FT)) ? ∂qsat∂T / max(hr, eps(FT)) : FT(0))
    return (q_src, ∂q_src∂T)
end

"""
    skin_flux_derivative(T_sfc, ρ_sfc, g_h, f_eff, w, ∂q_src∂T, ϵ, σ, thermo_params)

Return the derivative [W/m²/K] with respect to the skin temperature of the net
upward energy flux `LW_n + L + H` at the skin, given the surface air density
`ρ_sfc`, the bare-soil aerodynamic conductance for heat `g_h`, the conductance
ratio `f_eff` (see [`soil_conductance_ratio`](@ref)), the source weight `w` of
the surface humidity (see [`soil_surface_vapor_weight`](@ref)), and the
temperature derivative `∂q_src∂T` of the vapor source humidity.
"""
function skin_flux_derivative(
    T_sfc,
    ρ_sfc,
    g_h,
    f_eff,
    w,
    ∂q_src∂T,
    ϵ,
    σ,
    thermo_params,
)
    _LH_v0 = Thermodynamics.Parameters.LH_v0(thermo_params)
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    ∂L∂T = ρ_sfc * g_h * _LH_v0 * w * ∂q_src∂T
    ∂H∂T = ρ_sfc * g_h * f_eff * cp_d
    ∂LW_n∂T = 4 * ϵ * σ * T_sfc^3
    return ∂LW_n∂T + ∂L∂T + ∂H∂T
end

"""
    soil_skin_fluxes(T_sfc, T_a, f_eff, g_h, inputs, param_set, thermo_params,
                     ψ_sfc, Tf_depressed, g_liq, β_ice, ϵ, σ, earth_param_set)

Return `(; L, H, Λ)`: the latent and sensible heat fluxes [W/m²] from the soil
skin at temperature `T_sfc` as computed by SurfaceFluxes.jl with the interface
temperature and humidity (see the discussion at the top of this file), and the
derivative `Λ` [W/m²/K] of `LW_n + L + H` with respect to `T_sfc`. The air
density at the surface is evaluated as in SurfaceFluxes.jl, at the interface
temperature and humidity.

Called from [`update_soil_T_sfc_scheme`](@ref).
"""
function soil_skin_fluxes(
    T_sfc::FT,
    T_a::FT,
    f_eff::FT,
    g_h::FT,
    inputs,
    param_set,
    thermo_params,
    ψ_sfc::FT,
    Tf_depressed::FT,
    g_liq::FT,
    β_ice::FT,
    ϵ::FT,
    σ::FT,
    earth_param_set,
) where {FT}
    q_air = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
    (; qsat, ∂qsat∂T, frozen, hr) = soil_skin_state(
        T_sfc,
        inputs,
        thermo_params,
        param_set,
        ψ_sfc,
        Tf_depressed,
        earth_param_set,
    )
    (q_src, ∂q_src∂T) = soil_vapor_source(q_air, qsat, ∂qsat∂T, hr, frozen)
    w = soil_surface_vapor_weight(
        q_air,
        q_src,
        g_liq,
        g_h,
        f_eff,
        β_ice,
        frozen,
    )
    q_sfc = w * q_src + (1 - w) * q_air
    T_i = soil_interface_temperature(T_sfc, T_a, f_eff)
    ρ_sfc = SurfaceFluxes.surface_density(
        param_set,
        inputs.T_int,
        inputs.ρ_int,
        T_i,
        inputs.Δz,
        inputs.q_tot_int,
        inputs.q_liq_int,
        inputs.q_ice_int,
        q_sfc,
    )
    E = SurfaceFluxes.evaporation(
        param_set,
        inputs,
        g_h,
        inputs.q_tot_int,
        q_sfc,
        ρ_sfc,
        inputs.moisture_model,
    )
    L = SurfaceFluxes.latent_heat_flux(
        param_set,
        inputs,
        E,
        inputs.moisture_model,
    )
    H = SurfaceFluxes.sensible_heat_flux(
        param_set,
        inputs,
        g_h,
        inputs.T_int,
        T_i,
        ρ_sfc,
        E,
    )
    Λ = skin_flux_derivative(
        T_sfc,
        ρ_sfc,
        g_h,
        f_eff,
        w,
        ∂q_src∂T,
        ϵ,
        σ,
        thermo_params,
    )
    return (; L, H, Λ)
end

"""
    update_soil_T_sfc_scheme(ζ, param_set, thermo_params, inputs, scheme, u_star, z_0m, z_0b,
                             T_top, r, ϵ, σ, SW_n, LW_d, g_liq, β_ice, ψ_sfc, Tf_depressed,
                             W, r_under, earth_param_set)

Return the interface temperature [K] (see [`soil_interface_temperature`](@ref))
of the soil skin temperature that solves the skin energy balance
`r (SW_n + LW_n(T) + L(T) + H(T)) + (T - T_top) = 0` for the conductance `g_h`
of the current Monin-Obukhov state, with a fixed number of Newton iterations
starting from `T_top`, capped at the depressed freezing temperature
`Tf_depressed` when the top cell is frozen (contains ice, `β_ice > 0`, and
`T_top < Tf_depressed`). Used as the `update_T` callback of
`SurfaceFluxes.surface_fluxes`, which passes the first eight arguments. The
result depends only on the current Monin-Obukhov state, so that the mapping
between skin and interface temperature is consistent when the conductance
changes between iterations.

# Arguments
- `ζ`, `param_set`, `thermo_params`, `inputs`, `scheme`, `u_star`, `z_0m`,
  `z_0b`: Stability parameter [-], parameter sets, flux inputs, scheme,
  friction velocity [m/s], and roughness lengths [m], from SurfaceFluxes.
- `T_top`: Temperature of the node below the skin [K].
- `r`: Thermal resistance between the skin and that node [m² K/W].
- `ϵ`, `σ`: Surface emissivity [-] and Stefan-Boltzmann constant [W/m²/K⁴].
- `SW_n`, `LW_d`: Net shortwave (positive upward) and downwelling longwave
  radiation at the surface [W/m²].
- `g_liq`, `β_ice`, `ψ_sfc`, `Tf_depressed`: Dry-soil-layer vapor conductance
  [m/s], ice weight of the surface humidity [-], soil water potential [m], and
  depressed freezing temperature [K] at the surface.
- `W`, `r_under`: Canopy gap fraction [-] and under-canopy resistance [s/m]
  (see [`soil_conductance_ratio`](@ref)).
- `earth_param_set`: Land parameters.

Called from [`solve_soil_surface_temperature_at_a_point`](@ref).
"""
function update_soil_T_sfc_scheme(
    ζ,
    param_set,
    thermo_params,
    inputs,
    scheme,
    u_star,
    z_0m,
    z_0b,
    T_top,
    r,
    ϵ,
    σ,
    SW_n,
    LW_d,
    g_liq,
    β_ice,
    ψ_sfc,
    Tf_depressed,
    W,
    r_under,
    earth_param_set,
)
    g_h = SurfaceFluxes.heat_conductance(
        param_set,
        ζ,
        u_star,
        inputs,
        z_0m,
        z_0b,
        scheme,
    )
    f_eff = soil_conductance_ratio(W, r_under, g_h)
    T_a = adiabatic_surface_air_temperature(inputs, param_set, thermo_params)
    # Keyed on T_top, so trace ice in a cell above the melting point does not
    # pin the skin
    frozen_top = β_ice > 0 && T_top < Tf_depressed
    T_max = frozen_top ? Tf_depressed : oftype(T_top, Inf)
    T_sfc = T_top
    for _ in 1:3
        (; L, H, Λ) = soil_skin_fluxes(
            T_sfc,
            T_a,
            f_eff,
            g_h,
            inputs,
            param_set,
            thermo_params,
            ψ_sfc,
            Tf_depressed,
            g_liq,
            β_ice,
            ϵ,
            σ,
            earth_param_set,
        )
        LW_n = -ϵ * (LW_d - σ * T_sfc^4)
        ΔT = -(r * (SW_n + LW_n + L + H) + (T_sfc - T_top)) / (r * Λ + 1)
        T_sfc = min(T_sfc + ΔT, T_max)
    end
    return soil_interface_temperature(T_sfc, T_a, f_eff)
end

"""
    update_soil_q_vap_sfc_scheme(ζ, param_set, thermo_params, inputs, scheme, T_i, u_star,
                                 z_0m, z_0b, g_liq, β_ice, ψ_sfc, Tf_depressed, W, r_under,
                                 earth_param_set)

Return the surface specific humidity [kg/kg] of the soil passed to
SurfaceFluxes.jl, given the interface temperature `T_i` returned by
[`update_soil_T_sfc_scheme`](@ref) for the same Monin-Obukhov state. Used as
the `update_q` callback of `SurfaceFluxes.surface_fluxes` in the skin
temperature solve, which passes the first nine arguments.

Called from [`solve_soil_surface_temperature_at_a_point`](@ref).
"""
function update_soil_q_vap_sfc_scheme(
    ζ,
    param_set,
    thermo_params,
    inputs,
    scheme,
    T_i,
    u_star,
    z_0m,
    z_0b,
    g_liq,
    β_ice,
    ψ_sfc,
    Tf_depressed,
    W,
    r_under,
    earth_param_set,
)
    g_h = SurfaceFluxes.heat_conductance(
        param_set,
        ζ,
        u_star,
        inputs,
        z_0m,
        z_0b,
        scheme,
    )
    f_eff = soil_conductance_ratio(W, r_under, g_h)
    T_a = adiabatic_surface_air_temperature(inputs, param_set, thermo_params)
    T_sfc = soil_skin_temperature(T_i, T_a, f_eff)
    q_air = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
    (; qsat, ∂qsat∂T, frozen, hr) = soil_skin_state(
        T_sfc,
        inputs,
        thermo_params,
        param_set,
        ψ_sfc,
        Tf_depressed,
        earth_param_set,
    )
    (q_src, _) = soil_vapor_source(q_air, qsat, ∂qsat∂T, hr, frozen)
    w = soil_surface_vapor_weight(
        q_air,
        q_src,
        g_liq,
        g_h,
        f_eff,
        β_ice,
        frozen,
    )
    return w * q_src + (1 - w) * q_air
end

"""
    solve_soil_surface_temperature_at_a_point(return_extra_fluxes, return_sensitivity, T_below,
                                              r, ϵ, SW_n, LW_d, g_liq, β_ice, ψ_sfc,
                                              Tf_depressed, W, r_under, h_sfc, displ, P_atmos,
                                              T_atmos, q_atmos, u_atmos, roughness_model,
                                              atmos_h, gustiness, earth_param_set)

Solve the soil skin surface energy balance at a point within the Monin-Obukhov
iterations of SurfaceFluxes.jl, starting from `T_below`, and return the skin
temperature with the turbulent fluxes at it.

# Arguments
- `return_extra_fluxes`: `Val(true)` to also return the momentum and buoyancy
  fluxes, which a coupled atmosphere needs.
- `return_sensitivity`: `Val(true)` to also return the sensitivity `∂F∂T` of
  the atmospheric flux to `T_below`, which a [`SlabLitter`](@ref) needs.
- `T_below`: Temperature of the node below the skin (top soil layer or
  litter) [K].
- `r`: Thermal resistance between the skin and that node [m² K/W].
- `ϵ`: Surface emissivity [-].
- `SW_n`, `LW_d`: Net shortwave (positive upward) and downwelling longwave
  radiation at the surface [W/m²].
- `g_liq`: Conductance of the dry soil layer to water vapor [m/s].
- `β_ice`: Ice fraction weighting the surface humidity when frozen [-].
- `ψ_sfc`: Soil water potential at the surface [m].
- `Tf_depressed`: Depressed freezing temperature at the surface [K].
- `W`, `r_under`: Canopy gap fraction [-] and under-canopy resistance [s/m]
  (see [`soil_conductance_ratio`](@ref)).
- `h_sfc`, `displ`: Surface height and displacement height [m].
- `P_atmos`, `T_atmos`, `q_atmos`, `u_atmos`: Atmospheric pressure [Pa],
  temperature [K], specific humidity [kg/kg], and wind (a speed or a
  horizontal vector) [m/s] at the reference height `atmos_h` [m].
- `roughness_model`, `gustiness`: SurfaceFluxes roughness and gustiness
  specifications.
- `earth_param_set`: Land parameters.

# Returns
`(; lhf, shf, vapor_flux_liq, vapor_flux_ice, T_sfc)`, with `ρτxz`, `ρτyz`,
and `buoyancy_flux` before `T_sfc` if `return_extra_fluxes` is `Val(true)`: the
latent and sensible heat fluxes [W/m²], the vapor flux [m/s of liquid water]
from liquid water or from ice (sublimation at and below the freezing
temperature), and the skin temperature [K]; for `return_sensitivity` =
`Val(true)`, followed by `∂F∂T = Λ / (1 + r Λ)` [W/m²/K], the derivative of the
net upward atmospheric energy flux at the skin with respect to `T_below`, with
`Λ = ∂(LW_n + L + H)/∂T_sfc`. The fluxes are those computed by SurfaceFluxes.jl
from the interface temperature and humidity, which equal those from the skin
with the under-canopy conductance.

Called from [`update_soil_surface_temperature!`](@ref).
"""
function solve_soil_surface_temperature_at_a_point(
    return_extra_fluxes::Val,
    return_sensitivity::Val,
    T_below::FT,
    r::FT,
    ϵ::FT,
    SW_n::FT,
    LW_d::FT,
    g_liq::FT,
    β_ice::FT,
    ψ_sfc::FT,
    Tf_depressed::FT,
    W::FT,
    r_under::FT,
    h_sfc::FT,
    displ::FT,
    P_atmos::FT,
    T_atmos::FT,
    q_atmos::FT,
    u_atmos,
    roughness_model,
    atmos_h::FT,
    gustiness,
    earth_param_set,
) where {FT}
    config = SurfaceFluxes.SurfaceFluxConfig(roughness_model, gustiness)
    positional_default_args = (
        scheme = SurfaceFluxes.PointValueScheme(),
        solver_opts = nothing,
        flux_specs = nothing,
    )
    # A coupled atmosphere provides a wind vector; a prescribed one a speed
    u = u_atmos isa FT ? (u_atmos, FT(0)) : u_atmos
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    _grav = LP.grav(earth_param_set)
    _σ = LP.Stefan(earth_param_set)
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    Δz = atmos_h - h_sfc
    ρ_atmos =
        Thermodynamics.air_density(thermo_params, T_atmos, P_atmos, q_atmos)
    update_T(args...) = update_soil_T_sfc_scheme(
        args...,
        T_below,
        r,
        ϵ,
        _σ,
        SW_n,
        LW_d,
        g_liq,
        β_ice,
        ψ_sfc,
        Tf_depressed,
        W,
        r_under,
        earth_param_set,
    )
    update_q(args...) = update_soil_q_vap_sfc_scheme(
        args...,
        g_liq,
        β_ice,
        ψ_sfc,
        Tf_depressed,
        W,
        r_under,
        earth_param_set,
    )
    ρ_sfc = ClimaLand.compute_ρ_sfc(
        surface_flux_params,
        T_atmos,
        P_atmos,
        q_atmos,
        Δz,
        T_below,
    )
    q_sfc_guess = soil_specific_humidity(
        T_below,
        ρ_sfc,
        ψ_sfc,
        Tf_depressed,
        earth_param_set,
    )
    # The update_T callback restarts from T_below at each Monin-Obukhov
    # iteration, so the initial guess only seeds the first conductance
    output = SurfaceFluxes.surface_fluxes(
        surface_flux_params,
        T_atmos,
        q_atmos,
        FT(0),#phase_partition_atmos.liq,
        FT(0),#,phase_partition_atmos.ice,
        ρ_atmos,
        T_below,
        q_sfc_guess,
        _grav * h_sfc,
        Δz,
        displ,
        u,
        (FT(0), FT(0)), # u_sfc
        nothing, # roughness inputs
        config,
        positional_default_args...,
        update_T,
        update_q,
    )
    # Recover the skin temperature from the interface temperature; T_a is as
    # in `adiabatic_surface_air_temperature`
    g_h = output.g_h
    f_eff = soil_conductance_ratio(W, r_under, g_h)
    T_a = T_atmos + _grav * Δz / cp_d
    T_sfc = soil_skin_temperature(output.T_sfc, T_a, f_eff)
    # Vapor flux in volume of liquid water
    Ẽ = output.evaporation / LP.ρ_cloud_liq(earth_param_set)
    is_liquid = ClimaLand.heaviside(T_sfc, Tf_depressed)
    fluxes = (;
        lhf = output.lhf,
        shf = output.shf,
        vapor_flux_liq = Ẽ * is_liquid,
        vapor_flux_ice = Ẽ * (1 - is_liquid),
    )
    solution = skin_solution(
        return_extra_fluxes,
        fluxes,
        output,
        T_sfc,
        surface_flux_params,
        T_atmos,
        ρ_atmos,
        q_atmos,
        Δz,
    )
    return with_skin_sensitivity(
        return_sensitivity,
        solution,
        T_sfc,
        g_h,
        f_eff,
        r,
        ϵ,
        g_liq,
        β_ice,
        ψ_sfc,
        Tf_depressed,
        surface_flux_params,
        P_atmos,
        T_atmos,
        q_atmos,
        Δz,
        earth_param_set,
    )
end

"""
    skin_solution(return_extra_fluxes, fluxes, output, T_sfc, surface_flux_params, T_atmos,
                  ρ_atmos, q_atmos, Δz)

Return the NamedTuple stored in `p.soil.turbulent_fluxes`: the `fluxes` with
the skin temperature `T_sfc`, and, for `Val(true)`, the momentum fluxes and
the buoyancy flux of the SurfaceFluxes `output` before it.
"""
skin_solution(::Val{false}, fluxes, output, T_sfc, args...) =
    (; fluxes..., T_sfc)
function skin_solution(
    ::Val{true},
    fluxes,
    output,
    T_sfc,
    surface_flux_params,
    T_atmos,
    ρ_atmos,
    q_atmos,
    Δz,
)
    FT = typeof(T_sfc)
    ρ_sfc = SurfaceFluxes.surface_density(
        surface_flux_params,
        T_atmos,
        ρ_atmos,
        output.T_sfc,
        Δz,
        q_atmos,
        FT(0),
        FT(0),
        output.q_vap_sfc,
    )
    buoyancy_flux = SurfaceFluxes.buoyancy_flux(
        surface_flux_params,
        output.shf,
        output.lhf,
        output.T_sfc,
        ρ_sfc,
        output.q_vap_sfc,
        FT(0),
        FT(0),
    )
    return (;
        fluxes...,
        ρτxz = output.ρτxz,
        ρτyz = output.ρτyz,
        buoyancy_flux,
        T_sfc,
    )
end

"""
    with_skin_sensitivity(return_sensitivity, solution, T_sfc, g_h, f_eff, r, ϵ, g_liq, β_ice,
                          ψ_sfc, Tf_depressed, surface_flux_params, P_atmos, T_atmos,
                          q_atmos, Δz, earth_param_set)

Return the skin `solution`, followed for `Val(true)` by `∂F∂T = Λ / (1 + r Λ)`,
the derivative of the net upward atmospheric flux at the skin with respect to
the temperature of the node below it, with `Λ = ∂(LW_n + L + H)/∂T_sfc` at the
skin temperature `T_sfc`, bare-soil conductance `g_h`, and conductance ratio
`f_eff` (see [`soil_conductance_ratio`](@ref)).
"""
with_skin_sensitivity(::Val{false}, solution, args...) = solution
function with_skin_sensitivity(
    ::Val{true},
    solution,
    T_sfc,
    g_h,
    f_eff,
    r,
    ϵ,
    g_liq,
    β_ice,
    ψ_sfc,
    Tf_depressed,
    surface_flux_params,
    P_atmos,
    T_atmos,
    q_atmos,
    Δz,
    earth_param_set,
)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    (; ρ_sfc, qsat, ∂qsat∂T, frozen, hr) = soil_skin_state(
        T_sfc,
        T_atmos,
        P_atmos,
        q_atmos,
        Δz,
        surface_flux_params,
        thermo_params,
        ψ_sfc,
        Tf_depressed,
        earth_param_set,
    )
    (q_src, ∂q_src∂T) = soil_vapor_source(q_atmos, qsat, ∂qsat∂T, hr, frozen)
    w = soil_surface_vapor_weight(
        q_atmos,
        q_src,
        g_liq,
        g_h,
        f_eff,
        β_ice,
        frozen,
    )
    Λ = skin_flux_derivative(
        T_sfc,
        ρ_sfc,
        g_h,
        f_eff,
        w,
        ∂q_src∂T,
        ϵ,
        LP.Stefan(earth_param_set),
        thermo_params,
    )
    return (; solution..., ∂F∂T = Λ / (1 + r * Λ))
end

"""
    soil_surface_vapor_conductance!(g_soil_sfc, model::EnergyHydrology, Y, p)

Compute the conductance [m/s] of the dry soil layer (and any [`SlabLitter`](@ref)
surface layer in series) to water vapor at the soil surface into `g_soil_sfc`
and return it. The liquid water content is extrapolated to the surface from the
top two cells (capped by the top-cell value), and its effective saturation
within the ice-free pore space `ν - θ_i` (cf. the effective porosity of CLM5)
is clamped to `[0, 1]` without the floor above `θ_r` used elsewhere, so that
the conductance vanishes once the mobile water is exhausted.

Called from [`update_soil_surface_temperature!`](@ref) and
`get_update_surface_humidity_function(::EnergyHydrology)`.
"""
function soil_surface_vapor_conductance!(
    g_soil_sfc,
    model::EnergyHydrology,
    Y,
    p,
)
    FT = eltype(Y)
    (; ν, θ_r, d_ds, evap_p, evap_α, hydrology_cm, earth_param_set) =
        model.parameters
    hydrology_cm_sfc = ClimaLand.Domains.top_center_to_surface(hydrology_cm)
    S_c_sfc = hydrology_cm_sfc.S_c
    θ_r_sfc = ClimaLand.Domains.top_center_to_surface(θ_r)
    ν_top = ClimaLand.Domains.top_center_to_surface(ν)
    θ_i_top = ClimaLand.Domains.top_center_to_surface(Y.soil.θ_i)
    sqrt_ε = sqrt(eps(FT))
    ν_sfc = @. lazy(max(ν_top - θ_i_top, θ_r_sfc + sqrt_ε))
    ϑ_l_top = ClimaLand.Domains.top_center_to_surface(Y.soil.ϑ_l)
    θ_l_sfc = g_soil_sfc
    ClimaLand.Domains.linear_interpolation_to_surface!(
        θ_l_sfc,
        Y.soil.ϑ_l,
        model.domain.fields.z,
        model.domain.fields.Δz_top,
    )
    S_l_sfc = g_soil_sfc # currently set to θ_l_sfc
    S_min, S_max = FT(0), FT(1)
    # overwrite with S_l_sfc ∈ [0, 1], unclipped at θ_r
    @. S_l_sfc = clamp(
        (min(θ_l_sfc, ϑ_l_top) - θ_r_sfc) / (ν_sfc - θ_r_sfc),
        S_min,
        S_max,
    )
    _D_vapor = FT(LP.D_vapor(earth_param_set))
    # currently set to S_l_sfc; overwrite with the conductance
    g_soil_sfc .= soil_conductance.(
        S_l_sfc,
        S_c_sfc,
        d_ds,
        evap_p,
        evap_α,
        _D_vapor,
        ν_top,
        θ_r_sfc,
        θ_i_top,
    )
    add_litter_vapor_resistance!(g_soil_sfc, model.surface_layer, Y, _D_vapor)
    # Reusing g_soil_sfc for the intermediates keeps the kernel argument
    # count within the parameter memory limit of P100 GPUs
    return g_soil_sfc
end

"""
    update_soil_surface_temperature!(model::EnergyHydrology, SW_n, LW_d, r_litter, Y, p, t)
    update_soil_surface_temperature!(model::EnergyHydrology, Y, p, t)

Solve for the soil skin temperature from the surface energy balance and store
it, with the turbulent fluxes at it from the same Monin-Obukhov solve, in
`p.soil.turbulent_fluxes`; return `nothing`. A no-op unless the top boundary
condition is an `AtmosDrivenFluxBC`.

# Arguments
- `SW_n`: Net shortwave radiation at the soil surface, positive upward, i.e.
  minus the absorbed shortwave radiation [W/m²].
- `LW_d`: Downwelling longwave radiation at the soil surface [W/m²].
- `r_litter`: Litter thermal resistance used with a [`NoLitter`](@ref)
  surface layer [m² K/W].

`SW_n`, `LW_d`, and `r_litter` may be fields or lazy broadcasts; land models
with a canopy pass the radiation transmitted and emitted by the canopy and the
canopy's `biomass.r_litter`. The turbulent fluxes use the canopy gap fraction
`p.soil.W_gap` and the under-canopy resistance `p.soil.r_undercanopy`, which
must be up to date. The four-argument method is for bare soil exposed to the
sky: it sets `W_gap` and `r_undercanopy` to their bare-soil values (1 and 0)
and uses the downwelling radiation in `p.drivers`, the soil albedo, and no
litter resistance. The skin exchanges with the node below it given by
[`skin_lower_node`](@ref); with a [`SlabLitter`](@ref) surface layer, the
sensitivity of the atmospheric flux to the litter temperature is stored in
`p.soil.skin_solve` as well (see [`store_skin_solution!`](@ref)). The
atmospheric state at the reference height `atmos.h` is read from `p.drivers`,
whether prescribed or supplied by a coupler, and `p.soil.sfc_scratch` is
overwritten.

Called from the `soil_boundary_fluxes!` methods and, in integrated models,
from `lsm_radiant_energy_fluxes!`. See also
[`solve_soil_surface_temperature_at_a_point`](@ref).
"""
function update_soil_surface_temperature!(
    model::EnergyHydrology,
    SW_n,
    LW_d,
    r_litter,
    Y,
    p,
    t,
)
    bc = model.boundary_conditions.top
    bc isa AtmosDrivenFluxBC || return nothing
    atmos = bc.atmos
    earth_param_set = model.parameters.earth_param_set
    ν_sfc = ClimaLand.Domains.top_center_to_surface(model.parameters.ν)
    θ_i_sfc = ClimaLand.Domains.top_center_to_surface(Y.soil.θ_i)
    ψ_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.ψ)
    Tf_depressed_sfc =
        ClimaLand.Domains.top_center_to_surface(p.soil.Tf_depressed)
    ϵ = ClimaLand.surface_emissivity(model, Y, p)
    T_below, r = skin_lower_node(model.surface_layer, model, Y, p, r_litter)
    g_soil_sfc =
        soil_surface_vapor_conductance!(p.soil.sfc_scratch, model, Y, p)
    h_sfc = ClimaLand.surface_height(model, Y, p)
    roughness_model = ClimaLand.surface_roughness_model(model, Y, p)
    displ = ClimaLand.surface_displacement_height(model, Y, p)
    gustiness = SurfaceFluxes.ConstantGustinessSpec(atmos.gustiness)
    return_extra_fluxes = Val(ClimaLand.return_momentum_fluxes(atmos))
    return_sensitivity = Val(model.surface_layer isa SlabLitter)
    skin = @. lazy(
        solve_soil_surface_temperature_at_a_point(
            return_extra_fluxes,
            return_sensitivity,
            T_below,
            r,
            ϵ,
            SW_n,
            LW_d,
            g_soil_sfc,
            (θ_i_sfc ./ ν_sfc) .^ 4, # β_ice
            ψ_sfc,
            Tf_depressed_sfc,
            p.soil.W_gap,
            p.soil.r_undercanopy,
            h_sfc,
            displ,
            p.drivers.P,
            p.drivers.T,
            p.drivers.q,
            p.drivers.u,
            roughness_model,
            atmos.h,
            gustiness,
            earth_param_set,
        ),
    )
    store_skin_solution!(model.surface_layer, p, skin)
    return nothing
end

function update_soil_surface_temperature!(model::EnergyHydrology, Y, p, t)
    model.boundary_conditions.top isa AtmosDrivenFluxBC || return nothing
    FT = eltype(Y)
    p.soil.W_gap .= 1
    p.soil.r_undercanopy .= 0
    α_sfc = ClimaLand.surface_albedo(model, Y, p)
    SW_n = @. lazy(-(1 - α_sfc) * p.drivers.SW_d)
    return update_soil_surface_temperature!(
        model,
        SW_n,
        p.drivers.LW_d,
        FT(0),
        Y,
        p,
        t,
    )
end
