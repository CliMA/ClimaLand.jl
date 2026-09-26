#=
Soil skin (surface) temperature.

The soil surface ("skin") is treated as a layer with zero heat capacity that
exchanges radiation and turbulent heat/vapor with the air above, and is
connected to the center of the top soil layer by a thermal resistance

    r = Δz_top / κ_top + r_litter  (m² K / W),

where Δz_top is the distance between the surface and the center of the top
layer, κ_top is the thermal conductivity of the top layer, and r_litter is the
(additional) thermal resistance of a litter/thatch layer on top of the mineral
soil. The skin temperature T_sfc satisfies the surface energy balance

    f(T_sfc) = SW_n + LW_n(T_sfc) + L(T_sfc) + H(T_sfc) - (T_top - T_sfc) / r = 0,

with all fluxes positive upward. It is solved for within the Monin-Obukhov
iterations of SurfaceFluxes.jl, in the same way as the snow surface temperature
(see `Snow.update_surf_temp!`). The turbulent fluxes, surface humidity, and
upwelling longwave radiation of the soil are then evaluated at T_sfc via
`component_temperature`, and the soil receives the heat flux (T_top - T_sfc)/r.

Under a canopy, the turbulent exchange between the ground and the atmosphere
is reduced. The conductance used for the soil turbulent fluxes is

    g_eff = W g_h + (1 - W) / (1/g_h + r_under),

where `g_h` is the conductance of bare soil computed by SurfaceFluxes.jl, `W`
is the canopy gap fraction (`p.soil.W_gap`), and `r_under` is the resistance
between the ground and the canopy air in a dense canopy
(`p.soil.r_undercanopy`), in series with the aerodynamic resistance to the
atmosphere (see `soil_conductance_ratio`). `W` and `r_under` are set by the
integrated land models with a canopy (see
`src/integrated/undercanopy_conductance.jl`). Without a canopy, `W = 1`,
`r_under = 0`, and `g_eff = g_h`. Because SurfaceFluxes.jl computes the fluxes with
`g_h`, the soil temperature and humidity passed to SurfaceFluxes.jl are those
at an "interface" such that the fluxes with conductance `g_h` from the interface
equal the fluxes with conductance `g_eff` from the skin:

    T_i = T_a + (g_eff / g_h) (T_sfc - T_a),

where `T_a = T_int + (Φ_int - Φ_sfc) / cp_d` is the atmospheric temperature
brought dry-adiabatically to the surface, and analogously for the humidity (see
`soil_surface_vapor_weight`). This is the same approach as used for the canopy
(see `get_update_surface_temperature_function(::CanopyModel, Y, p)`).
=#

"""
    soil_surface_temperature(bc, p)

Returns the soil surface temperature: the skin temperature `p.soil.T_sfc` for
atmospherically driven soil, and the top layer temperature otherwise.
"""
soil_surface_temperature(::AtmosDrivenFluxBC, p) = p.soil.T_sfc
soil_surface_temperature(_, p) =
    ClimaLand.Domains.top_center_to_surface(p.soil.T)

"""
    initialize_soil_surface_temperature!(bc, p)

Sets the skin temperature (if present) to the top layer temperature; it is
subsequently updated from the surface energy balance when the boundary fluxes
are computed. Also sets the canopy gap fraction `p.soil.W_gap` and the
under-canopy resistance `p.soil.r_undercanopy` to their bare soil values (1
and 0); integrated models with a canopy update these before the soil fluxes are
computed.
"""
function initialize_soil_surface_temperature!(::AtmosDrivenFluxBC, p)
    p.soil.T_sfc .= ClimaLand.Domains.top_center_to_surface(p.soil.T)
    p.soil.W_gap .= 1
    p.soil.r_undercanopy .= 0
    return nothing
end
initialize_soil_surface_temperature!(_, p) = nothing

"""
    soil_surface_thermal_resistance(Δz_top::FT, κ_top::FT, r_litter::FT) where {FT}

Returns the thermal resistance (m² K/W) between the soil skin and the center of
the top soil layer: the half-cell conduction resistance `Δz_top/κ_top` plus the
resistance of the litter layer `r_litter`.
"""
soil_surface_thermal_resistance(Δz_top::FT, κ_top::FT, r_litter::FT) where {FT} =
    Δz_top / max(κ_top, eps(FT)) + max(r_litter, FT(0))

"""
    soil_conductance_ratio(W::FT, r_under::FT, g_h::FT) where {FT}

Returns the ratio `f_eff = g_eff / g_h` of the conductance used for the
turbulent fluxes of the soil to the conductance of bare soil `g_h`, with

    g_eff = W g_h + (1 - W) / (1/g_h + r_under).

Here `W` is the canopy gap fraction and `r_under` the resistance (s/m) between
the ground and the canopy air in a dense canopy, which is in series with the
aerodynamic resistance between the canopy air and the atmosphere. The latter is
approximated by the bare-soil resistance `1/g_h` computed by SurfaceFluxes.jl,
so that it includes the effects of atmospheric stability. Then
`W ≤ f_eff ≤ 1`, and `f_eff = 1` for bare soil (`W = 1` or `r_under = 0`).
"""
soil_conductance_ratio(W::FT, r_under::FT, g_h::FT) where {FT} =
    W + (1 - W) / (1 + g_h * r_under)

"""
    adiabatic_surface_air_temperature(inputs, param_set, thermo_params)

Returns the temperature of the air at the atmospheric reference height brought
dry-adiabatically to the surface, `T_int + (Φ_int - Φ_sfc) / cp_d`. The sensible
heat flux computed by SurfaceFluxes.jl is proportional to the difference of the
surface temperature from this value.
"""
function adiabatic_surface_air_temperature(inputs, param_set, thermo_params)
    Φ_sfc = SurfaceFluxes.surface_geopotential(inputs)
    Φ_int = SurfaceFluxes.interior_geopotential(param_set, inputs)
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    return inputs.T_int + (Φ_int - Φ_sfc) / cp_d
end

"""
    soil_interface_temperature(T_skin, T_a, f_eff)

Returns the temperature passed to SurfaceFluxes.jl, `T_a + f_eff (T_skin - T_a)`,
such that the sensible heat flux with the bare soil conductance `g_h` equals
that from the skin with conductance `g_eff = f_eff g_h`.
"""
soil_interface_temperature(T_skin, T_a, f_eff) = T_a + f_eff * (T_skin - T_a)

"""
    soil_skin_temperature(T_i, T_a, f_eff)

Inverse of `soil_interface_temperature`.
"""
soil_skin_temperature(T_i::FT, T_a::FT, f_eff::FT) where {FT} =
    T_a + (T_i - T_a) / max(f_eff, eps(FT))

"""
    soil_surface_vapor_weight(q_air, q_src, g_liq, g_h, f_eff, β_ice, frozen)

Returns the weight `w` such that the soil surface specific humidity passed to
SurfaceFluxes.jl is `q_sfc = w * q_src + (1 - w) * q_air`, where `q_src` is the
specific humidity of the vapor source (see `soil_vapor_source`). The vapor flux
is then `ρ g_h (q_sfc - q_air) = ρ w g_h (q_src - q_air)`:
- For evaporation from unfrozen soil, `w g_h` is the conductance of the dry
  soil layer `g_liq` in series with the turbulent conductance `g_eff = f_eff g_h`,
  so `w = f_eff g_liq / (f_eff g_h + g_liq)`.
- For sublimation from frozen soil, `w = β_ice f_eff`.
- For dew or frost formation (`q_air > q_src`), which occurs at the surface
  and is not limited by the dry soil layer, `w = f_eff` (as in CLM5).
The flux vanishes at `q_air = q_src`, so it is continuous in `q_air`.

This is used both in the skin temperature solve and in
`ClimaLand.get_update_surface_humidity_function(::EnergyHydrology, Y, p)`.
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
    soil_skin_state(T_sfc, inputs, thermo_params, param_set, ψ_sfc, Tf_depressed, earth_param_set)

Helper returning the surface air density, the (soil water potential adjusted)
saturation specific humidity at `T_sfc`, its derivative with respect to
temperature, whether the surface is frozen, and the Kelvin factor.
"""
function soil_skin_state(
    T_sfc::FT,
    inputs,
    thermo_params,
    param_set,
    ψ_sfc::FT,
    Tf_depressed::FT,
    earth_param_set,
) where {FT}
    T_atmos = inputs.T_int
    ρ_atmos = inputs.ρ_int
    q_atmos = inputs.q_tot_int
    P_atmos =
        Thermodynamics.air_pressure(thermo_params, T_atmos, ρ_atmos, q_atmos)
    ρ_sfc = ClimaLand.compute_ρ_sfc(
        param_set,
        T_atmos,
        P_atmos,
        q_atmos,
        inputs.Δz,
        T_sfc,
    )
    qsat = soil_specific_humidity(
        T_sfc,
        ρ_sfc,
        ψ_sfc,
        Tf_depressed,
        earth_param_set,
    )
    frozen = T_sfc < Tf_depressed
    # Kelvin factor (qsat above already includes it for unfrozen soil)
    hr = frozen ? FT(1) : soil_kelvin_factor(T_sfc, ψ_sfc, earth_param_set)
    LH =
        frozen ? Thermodynamics.latent_heat_sublim(thermo_params, T_sfc) :
        Thermodynamics.latent_heat_vapor(thermo_params, T_sfc)
    # The (weak) temperature dependence of the soil water potential factor is neglected
    ∂qsat∂T = Thermodynamics.∂q_vap_sat_∂T_from_L(thermo_params, qsat, LH, T_sfc)
    return (; ρ_sfc, qsat, ∂qsat∂T, frozen, hr)
end

"""
    soil_vapor_source(q_air, qsat, ∂qsat∂T, hr, frozen)

Returns the specific humidity of the vapor source at the soil surface and its
derivative with respect to the skin temperature, consistent with
`effective_soil_vapor_humidity`. For frozen soil, the source is saturated with
respect to ice.
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
    soil_skin_fluxes(T_sfc, T_a, f_eff, g_h, inputs, param_set, thermo_params,
                     ψ_sfc, Tf_depressed, g_liq, β_ice, earth_param_set)

Returns the latent and sensible heat fluxes from the soil skin at temperature
`T_sfc` as computed by SurfaceFluxes.jl with the interface temperature and
humidity (see the discussion at the top of this file), their derivatives with
respect to `T_sfc`, and the surface specific humidity passed to SurfaceFluxes.jl.
The air density at the surface is evaluated as in SurfaceFluxes.jl, at the
interface temperature and humidity.
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
    w = soil_surface_vapor_weight(q_air, q_src, g_liq, g_h, f_eff, β_ice, frozen)
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
    _LH_v0 = Thermodynamics.Parameters.LH_v0(thermo_params)
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    ∂L∂T = ρ_sfc * g_h * _LH_v0 * w * ∂q_src∂T
    ∂H∂T = ρ_sfc * g_h * f_eff * cp_d
    return (; L, H, ∂L∂T, ∂H∂T, q_sfc)
end

"""
    update_soil_T_sfc_scheme(ζ, param_set, thermo_params, inputs, scheme, u_star, z_0m, z_0b,
                             T_top, r, ϵ, σ, SW_n, LW_d, g_liq, β_ice, ψ_sfc, Tf_depressed,
                             W, r_under, earth_param_set)

The `update_T` callback of `SurfaceFluxes.surface_fluxes` in the soil skin
temperature solve. For the conductance `g_h` given by the current Monin-Obukhov
state, the skin energy balance
`r (SW_n + LW_n(T) + L(T) + H(T)) + (T - T_top) = 0` is solved with a fixed
number of Newton iterations starting from `T_top`, and the corresponding
interface temperature (see the discussion at the top of this file) is returned.
The result depends only on the current Monin-Obukhov state, not on the
previous iterate, so that the mapping between skin and interface temperature is
consistent even when the conductance changes between iterations.
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
    T_sfc = T_top
    for _ in 1:3
        (; L, H, ∂L∂T, ∂H∂T) = soil_skin_fluxes(
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
            earth_param_set,
        )
        LW_n = -ϵ * (LW_d - σ * T_sfc^4)
        ∂LW_n∂T = 4 * ϵ * σ * T_sfc^3
        T_sfc -=
            (r * (SW_n + LW_n + L + H) + (T_sfc - T_top)) /
            (r * (∂LW_n∂T + ∂L∂T + ∂H∂T) + 1)
    end
    return soil_interface_temperature(T_sfc, T_a, f_eff)
end

"""
    update_soil_q_vap_sfc_scheme(ζ, param_set, thermo_params, inputs, scheme, T_i, u_star,
                                 z_0m, z_0b, g_liq, β_ice, ψ_sfc, Tf_depressed, W, r_under,
                                 earth_param_set)

Surface specific humidity of the soil passed to SurfaceFluxes.jl, given the
interface temperature `T_i` returned by `update_soil_T_sfc_scheme` for the same
Monin-Obukhov state; used as the `update_q` callback of
`SurfaceFluxes.surface_fluxes` in the skin temperature solve.
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
    w = soil_surface_vapor_weight(q_air, q_src, g_liq, g_h, f_eff, β_ice, frozen)
    return w * q_src + (1 - w) * q_air
end

"""
    solve_soil_surface_temperature_at_a_point(T_top, r, ϵ, SW_n, LW_d, g_liq, β_ice, ψ_sfc,
                                              Tf_depressed, W, r_under, h_sfc, displ,
                                              P_atmos, T_atmos, q_atmos, u_atmos,
                                              roughness_model, atmos_h, gustiness,
                                              earth_param_set)

Solves the soil skin surface energy balance at a point and returns the skin
temperature. The initial guess is the top layer temperature `T_top`.
"""
function solve_soil_surface_temperature_at_a_point(
    T_top::FT,
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
)::FT where {FT}
    config = SurfaceFluxes.SurfaceFluxConfig(roughness_model, gustiness)
    positional_default_args = (
        scheme = SurfaceFluxes.PointValueScheme(),
        solver_opts = nothing,
        flux_specs = nothing,
    )
    # u is already a vector when we get it from a coupled atmosphere, otherwise we need to make it one
    u = u_atmos isa FT ? (u_atmos, FT(0)) : u_atmos
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    _grav = LP.grav(earth_param_set)
    _σ = LP.Stefan(earth_param_set)
    ρ_atmos =
        Thermodynamics.air_density(thermo_params, T_atmos, P_atmos, q_atmos)
    update_T(args...) = update_soil_T_sfc_scheme(
        args...,
        T_top,
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
        atmos_h - h_sfc,
        T_top,
    )
    q_sfc_guess =
        soil_specific_humidity(T_top, ρ_sfc, ψ_sfc, Tf_depressed, earth_param_set)
    # The skin temperature callback does not depend on the guess (see
    # `update_soil_T_sfc_scheme`); T_a is as in `adiabatic_surface_air_temperature`.
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    T_a = T_atmos + _grav * (atmos_h - h_sfc) / cp_d
    T_guess = T_top
    output = SurfaceFluxes.surface_fluxes(
        surface_flux_params,
        T_atmos,
        q_atmos,
        FT(0),#phase_partition_atmos.liq,
        FT(0),#,phase_partition_atmos.ice,
        ρ_atmos,
        T_guess,
        q_sfc_guess,
        _grav * h_sfc,
        atmos_h - h_sfc,
        displ,
        u,
        (FT(0), FT(0)), # u_sfc
        nothing, # roughness inputs
        config,
        positional_default_args...,
        update_T,
        update_q,
    )
    # Recover the skin temperature from the interface temperature
    f_eff = soil_conductance_ratio(W, r_under, output.g_h)
    T_sfc = soil_skin_temperature(output.T_sfc, T_a, f_eff)
    return isfinite(T_sfc) ? T_sfc : T_top
end

"""
    soil_surface_vapor_conductance!(g_soil_sfc, model::EnergyHydrology, Y, p)

Computes the conductance (m/s) of the dry soil layer to water vapor at the soil
surface, in place, and returns it. This is used in the surface humidity
parameterization of the soil.
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
    # Ice-free pore space (cf. the effective porosity used in CLM5)
    ν_top = ClimaLand.Domains.top_center_to_surface(ν)
    θ_i_top = ClimaLand.Domains.top_center_to_surface(Y.soil.θ_i)
    ν_sfc = @. lazy(max(ν_top - θ_i_top, θ_r_sfc + sqrt(eps(FT))))
    θ_l_sfc = g_soil_sfc
    ClimaLand.Domains.linear_interpolation_to_surface!(
        θ_l_sfc,
        p.soil.θ_l,
        model.domain.fields.z,
        model.domain.fields.Δz_top,
    )
    @. θ_l_sfc = clamp(θ_l_sfc, θ_r_sfc + eps(FT), ν_sfc)
    S_l_sfc = g_soil_sfc # currently set to θ_l_sfc
    @. S_l_sfc = effective_saturation(ν_sfc, θ_l_sfc, θ_r_sfc) # overwrite with S_l_sfc
    _D_vapor = FT(LP.D_vapor(earth_param_set))
    # currently set to S_l_sfc; overwrite with the conductance
    g_soil_sfc .=
        soil_conductance.(
            S_l_sfc,
            S_c_sfc,
            d_ds,
            evap_p,
            evap_α,
            _D_vapor,
            ν_sfc,
            θ_r_sfc,
        )
    # the above is jumping through hoops so that we dont hit the parameter memory limit on P100...
    return g_soil_sfc
end

"""
    update_soil_surface_temperature!(model::EnergyHydrology, SW_n, LW_d, r_litter, Y, p, t)

Solves for the soil skin temperature `p.soil.T_sfc` from the surface energy
balance, given the net shortwave radiation at the soil surface `SW_n`
(positive upward, i.e., minus the absorbed shortwave radiation), the downwelling
longwave radiation at the soil surface `LW_d`, and the litter thermal
resistance `r_litter` (m² K/W). The arguments may be fields or lazy broadcasted
objects. The turbulent fluxes use the canopy gap fraction `p.soil.W_gap` and
the under-canopy resistance `p.soil.r_undercanopy`, which must be up to date.

The skin temperature is only solved for when the soil is driven by a
`PrescribedAtmosphere`; otherwise it is left equal to the top layer temperature.
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
    return update_soil_surface_temperature!(
        bc.atmos,
        model,
        SW_n,
        LW_d,
        r_litter,
        Y,
        p,
        t,
    )
end

update_soil_surface_temperature!(atmos, model, SW_n, LW_d, r_litter, Y, p, t) =
    nothing

function update_soil_surface_temperature!(
    atmos::ClimaLand.PrescribedAtmosphere,
    model::EnergyHydrology,
    SW_n,
    LW_d,
    r_litter,
    Y,
    p,
    t,
)
    earth_param_set = model.parameters.earth_param_set
    ν_sfc = ClimaLand.Domains.top_center_to_surface(model.parameters.ν)
    θ_i_sfc = ClimaLand.Domains.top_center_to_surface(Y.soil.θ_i)
    T_top = ClimaLand.Domains.top_center_to_surface(p.soil.T)
    κ_top = ClimaLand.Domains.top_center_to_surface(p.soil.κ)
    ψ_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.ψ)
    Tf_depressed_sfc =
        ClimaLand.Domains.top_center_to_surface(p.soil.Tf_depressed)
    Δz_top = model.domain.fields.Δz_top
    ϵ = ClimaLand.surface_emissivity(model, Y, p)
    # Note: this uses p.soil.sfc_scratch, which is also used (and overwritten) when
    # computing the turbulent fluxes of the soil. The skin solve below completes first.
    g_soil_sfc = soil_surface_vapor_conductance!(p.soil.sfc_scratch, model, Y, p)
    h_sfc = ClimaLand.surface_height(model, Y, p)
    roughness_model = ClimaLand.surface_roughness_model(model, Y, p)
    displ = ClimaLand.surface_displacement_height(model, Y, p)
    gustiness = SurfaceFluxes.ConstantGustinessSpec(atmos.gustiness)
    r = @. lazy(soil_surface_thermal_resistance(Δz_top, κ_top, r_litter))
    p.soil.T_sfc .=
        solve_soil_surface_temperature_at_a_point.(
            T_top,
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
        )
    return nothing
end
