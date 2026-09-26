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
are computed.
"""
initialize_soil_surface_temperature!(::AtmosDrivenFluxBC, p) =
    p.soil.T_sfc .= ClimaLand.Domains.top_center_to_surface(p.soil.T)
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
    soil_surface_vapor_weight(q_air, qsat, g_liq, g_h, β_ice, frozen)

Returns the weight `w` such that the surface specific humidity of the soil
is `q_sfc = w * qsat + (1 - w) * q_air`. This must be consistent with the
surface humidity parameterization in
`ClimaLand.get_update_surface_humidity_function(::EnergyHydrology, Y, p)`.
"""
function soil_surface_vapor_weight(
    q_air::FT,
    qsat::FT,
    g_liq::FT,
    g_h::FT,
    β_ice::FT,
    frozen::Bool,
) where {FT}
    if frozen
        return q_air < qsat ? β_ice : FT(1)
    else
        return (g_liq / g_h) / (1 + g_liq / g_h)
    end
end

"""
    soil_skin_state(T_sfc, inputs, thermo_params, param_set, ψ_sfc, Tf_depressed, earth_param_set)

Helper returning the surface air density, the (soil water potential adjusted)
saturation specific humidity at `T_sfc`, its derivative with respect to
temperature, and whether the surface is frozen.
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
    LH =
        frozen ? Thermodynamics.latent_heat_sublim(thermo_params, T_sfc) :
        Thermodynamics.latent_heat_vapor(thermo_params, T_sfc)
    # The (weak) temperature dependence of the soil water potential factor is neglected
    ∂qsat∂T = Thermodynamics.∂q_vap_sat_∂T_from_L(thermo_params, qsat, LH, T_sfc)
    return (; ρ_sfc, qsat, ∂qsat∂T, frozen)
end

"""
    update_soil_T_sfc_scheme(ζ, param_set, thermo_params, inputs, scheme, u_star, z_0m, z_0b,
                             T_top, r, ϵ, σ, SW_n, LW_d, g_liq, β_ice, ψ_sfc, Tf_depressed,
                             earth_param_set)

Newton update of the soil skin temperature used as the `update_T` callback of
`SurfaceFluxes.surface_fluxes`: returns `T + ΔT` with `ΔT = -f(T)/f'(T)`, where
`r f(T) = r (SW_n + LW_n(T) + L(T) + H(T)) + (T - T_top)`.
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
    earth_param_set,
)
    T_sfc = inputs.T_sfc_guess
    q_air = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
    (; ρ_sfc, qsat, ∂qsat∂T, frozen) = soil_skin_state(
        T_sfc,
        inputs,
        thermo_params,
        param_set,
        ψ_sfc,
        Tf_depressed,
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
    w = soil_surface_vapor_weight(q_air, qsat, g_liq, g_h, β_ice, frozen)
    q_sfc = w * qsat + (1 - w) * q_air
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
        T_sfc,
        ρ_sfc,
        E,
    )
    _LH_v0 = Thermodynamics.Parameters.LH_v0(thermo_params)
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    ∂L∂T = ρ_sfc * g_h * _LH_v0 * w * ∂qsat∂T
    ∂H∂T = ρ_sfc * g_h * cp_d
    LW_n = -ϵ * (LW_d - σ * T_sfc^4)
    ∂LW_n∂T = 4 * ϵ * σ * T_sfc^3
    ΔT =
        -(r * (SW_n + LW_n + L + H) + (T_sfc - T_top)) /
        (r * (∂LW_n∂T + ∂L∂T + ∂H∂T) + 1)
    return T_sfc + ΔT
end

"""
    update_soil_q_vap_sfc_scheme(ζ, param_set, thermo_params, inputs, scheme, T_sfc, u_star,
                                 z_0m, z_0b, g_liq, β_ice, ψ_sfc, Tf_depressed, earth_param_set)

Surface specific humidity of the soil at the skin temperature `T_sfc`, used as
the `update_q` callback of `SurfaceFluxes.surface_fluxes` in the skin
temperature solve.
"""
function update_soil_q_vap_sfc_scheme(
    ζ,
    param_set,
    thermo_params,
    inputs,
    scheme,
    T_sfc,
    u_star,
    z_0m,
    z_0b,
    g_liq,
    β_ice,
    ψ_sfc,
    Tf_depressed,
    earth_param_set,
)
    q_air = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
    (; qsat, frozen) = soil_skin_state(
        T_sfc,
        inputs,
        thermo_params,
        param_set,
        ψ_sfc,
        Tf_depressed,
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
    w = soil_surface_vapor_weight(q_air, qsat, g_liq, g_h, β_ice, frozen)
    return w * qsat + (1 - w) * q_air
end

"""
    solve_soil_surface_temperature_at_a_point(T_top, r, ϵ, SW_n, LW_d, g_liq, β_ice, ψ_sfc,
                                              Tf_depressed, h_sfc, displ, P_atmos, T_atmos,
                                              q_atmos, u_atmos, roughness_model, atmos_h,
                                              gustiness, earth_param_set)

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
        earth_param_set,
    )
    update_q(args...) = update_soil_q_vap_sfc_scheme(
        args...,
        g_liq,
        β_ice,
        ψ_sfc,
        Tf_depressed,
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
    output = SurfaceFluxes.surface_fluxes(
        surface_flux_params,
        T_atmos,
        q_atmos,
        FT(0),#phase_partition_atmos.liq,
        FT(0),#,phase_partition_atmos.ice,
        ρ_atmos,
        T_top,
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
    T_sfc = output.T_sfc
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
    ν_sfc = ClimaLand.Domains.top_center_to_surface(ν)
    θ_r_sfc = ClimaLand.Domains.top_center_to_surface(θ_r)
    θ_l_sfc = g_soil_sfc
    ClimaLand.Domains.linear_interpolation_to_surface!(
        θ_l_sfc,
        p.soil.θ_l,
        model.domain.fields.z,
        model.domain.fields.Δz_top,
    )
    @. θ_l_sfc = max(θ_l_sfc, θ_r_sfc + eps(FT))
    S_l_sfc = g_soil_sfc # currently set to θ_l_sfc
    @. S_l_sfc = effective_saturation(ν_sfc, θ_l_sfc, θ_r_sfc) # overwrite with S_l_sfc
    _D_vapor = FT(LP.D_vapor(earth_param_set))
    # currently set to S_l_sfc; overwrite with the conductance
    g_soil_sfc .=
        soil_conductance.(S_l_sfc, S_c_sfc, d_ds, evap_p, evap_α, _D_vapor)
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
objects.

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
