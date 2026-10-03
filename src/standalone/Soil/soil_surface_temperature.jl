#=
Soil skin (surface) temperature.

The soil surface ("skin") is treated as a layer with zero heat capacity that
exchanges radiation and turbulent heat/vapor with the air above, and is
connected to the center of the top soil layer by the half-cell conduction
resistance

    r = Δz_top / κ_top  (m² K / W),

where Δz_top is the distance between the surface and the center of the top
layer and κ_top is the thermal conductivity of the top layer. The skin
temperature T_sfc satisfies the surface energy balance

    f(T_sfc) = SW_n + LW_n(T_sfc) + L(T_sfc) + H(T_sfc) - (T_top - T_sfc)/r = 0,

with all fluxes positive upward. It is solved for within the Monin-Obukhov
iterations of SurfaceFluxes.jl, in the same way as the snow surface temperature
(see `Snow.update_surf_temp!`), and the same solve yields the turbulent fluxes
at T_sfc, which are stored with it in `p.soil.turbulent_fluxes`. The surface
humidity and upwelling longwave radiation of the soil are evaluated at T_sfc
via `component_temperature`. The soil receives the net radiation and turbulent
fluxes at T_sfc, which equal the conduction (T_top - T_sfc)/r to the tolerance
of the solve.
When the top cell is frozen (contains ice and is below the depressed freezing
temperature), T_sfc is capped at the depressed freezing temperature. The
atmospheric fluxes at the cap then exceed the conduction, and the excess energy
warms the top cell; its ice melts through the soil phase change once the cell
reaches the freezing temperature, which also releases the cap.
Conduction between the soil and a snowpack or lake sediment is unaffected and
uses the top cell temperature.
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
    frozen_soil_vapor_weight(θ_i, ν)

Return the weight `β_ice = (θ_i / ν)^4` [-] of the saturation specific
humidity in the surface specific humidity of a frozen soil losing water by
sublimation, given the volumetric ice content `θ_i` and the porosity `ν`
[m³/m³]. This is a heuristic without a published source: it suppresses
sublimation from a frozen surface whose pores hold little ice, in the same
way that the dry-soil-layer conductance suppresses evaporation from an
unfrozen surface holding little liquid water, and reaches one when the pores
are filled with ice. The exponent is chosen so that the suppression is strong
at intermediate ice contents.

Called from `update_soil_surface_temperature!` and
`get_update_surface_humidity_function(::EnergyHydrology)`.
"""
frozen_soil_vapor_weight(θ_i, ν) = (θ_i / ν)^4

"""
    soil_surface_vapor_weight(q_air, qsat, g_liq, g_h, β_ice, frozen)

Return the weight `w` [-] such that the surface specific humidity of the soil
is `q_sfc = w * qsat + (1 - w) * q_air`: the ice fraction `β_ice` over a frozen
surface losing water, one over a frozen surface gaining water, and the
conductance ratio `g_liq / (g_liq + g_h)` otherwise. This is the same
parameterization as `get_update_surface_humidity_function(::EnergyHydrology)`,
which `turbulent_fluxes!(dest, atmos, ::EnergyHydrology, Y, p, t)` uses.
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
    soil_skin_state(T_sfc, inputs, thermo_params, param_set, ψ_sfc, Tf_depressed, earth_param_set) -> NamedTuple

Return `(; ρ_sfc, qsat, ∂qsat∂T, frozen)`: the surface air density
[kg/m³], the saturation specific humidity at `T_sfc` adjusted for the soil
water potential [kg/kg], its derivative with respect to temperature [1/K],
and whether the surface is frozen.

Called from the SurfaceFluxes callbacks [`update_soil_T_sfc_scheme`](@ref) and
[`update_soil_q_vap_sfc_scheme`](@ref).
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
    # At T_sfc = Tf_depressed (the cap below), qsat is over ice, so the surface
    # is frozen there as well
    frozen = T_sfc <= Tf_depressed
    LH =
        frozen ? Thermodynamics.latent_heat_sublim(thermo_params, T_sfc) :
        Thermodynamics.latent_heat_vapor(thermo_params, T_sfc)
    # The weak temperature dependence of the water potential factor is neglected
    ∂qsat∂T =
        Thermodynamics.∂q_vap_sat_∂T_from_L(thermo_params, qsat, LH, T_sfc)
    return (; ρ_sfc, qsat, ∂qsat∂T, frozen)
end

"""
    update_soil_T_sfc_scheme(ζ, param_set, thermo_params, inputs, scheme, u_star, z_0m, z_0b,
                             T_top, r, ϵ, σ, SW_n, LW_d, g_liq, β_ice, ψ_sfc, Tf_depressed,
                             earth_param_set)

Return the Newton update of the soil skin temperature [K] from
`ClimaLand.surface_temperature_newton_update`, with the layer below at `T_top`,
capped at the depressed freezing temperature `Tf_depressed` when the top cell
is frozen (contains ice, `β_ice > 0`, and `T_top < Tf_depressed`). Used as the
`update_T` callback of `SurfaceFluxes.surface_fluxes`, which passes the first
eight arguments.

# Arguments
- `ζ`, `param_set`, `thermo_params`, `inputs`, `scheme`, `u_star`, `z_0m`,
  `z_0b`: Stability parameter [-], parameter sets, flux inputs, scheme,
  friction velocity [m/s], and roughness lengths [m], from SurfaceFluxes.
- `T_top`: Temperature of the top soil cell [K].
- `r`: Thermal resistance between the skin and the top cell [m² K/W].
- `ϵ`, `σ`: Surface emissivity [-] and Stefan-Boltzmann constant [W/m²/K⁴].
- `SW_n`, `LW_d`: Net shortwave (positive upward) and downwelling longwave
  radiation at the surface [W/m²].
- `g_liq`, `β_ice`, `ψ_sfc`, `Tf_depressed`: Dry-soil-layer vapor conductance
  [m/s], ice weight of the surface humidity [-], soil water potential [m], and
  depressed freezing temperature [K] at the surface.
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
    # Keyed on T_top, so trace ice in a cell above the melting point does not
    # pin the skin
    frozen_top = β_ice > 0 && T_top < Tf_depressed
    T_max = frozen_top ? Tf_depressed : oftype(T_sfc, Inf)
    return ClimaLand.surface_temperature_newton_update(
        inputs,
        param_set,
        thermo_params,
        g_h,
        ρ_sfc,
        w * qsat + (1 - w) * q_air, # q_sfc
        w * ∂qsat∂T, # ∂q_sfc∂T
        T_top,
        r,
        T_max,
        ϵ,
        σ,
        SW_n,
        LW_d,
    )
end

"""
    update_soil_q_vap_sfc_scheme(ζ, param_set, thermo_params, inputs, scheme, T_sfc, u_star,
                                 z_0m, z_0b, g_liq, β_ice, ψ_sfc, Tf_depressed, earth_param_set)

Return the surface specific humidity [kg/kg] of the soil at the skin
temperature `T_sfc`. Used as the `update_q` callback of
`SurfaceFluxes.surface_fluxes` in the skin temperature solve, which passes the
first nine arguments.

Called from [`solve_soil_surface_temperature_at_a_point`](@ref).
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
    solve_soil_surface_temperature_at_a_point(return_extra_fluxes, T_top, r, ϵ, SW_n, LW_d,
                                              g_liq, β_ice, ψ_sfc, Tf_depressed, h_sfc, displ,
                                              P_atmos, T_atmos, q_atmos, u_atmos,
                                              roughness_model, atmos_h, gustiness,
                                              update_∂T_sfc∂T, update_∂q_sfc∂T,
                                              earth_param_set)

Solve the soil skin surface energy balance at a point within the Monin-Obukhov
iterations of SurfaceFluxes.jl, starting from `T_top`, and return the skin
temperature with the turbulent fluxes at it.

# Arguments
- `return_extra_fluxes`: `Val(true)` to also return the momentum and buoyancy
  fluxes, which a coupled atmosphere needs.
- `T_top`: Temperature of the top soil cell [K].
- `r`: Thermal resistance between the skin and the top cell [m² K/W].
- `ϵ`: Surface emissivity [-].
- `SW_n`, `LW_d`: Net shortwave (positive upward) and downwelling longwave
  radiation at the surface [W/m²].
- `g_liq`: Conductance of the dry soil layer to water vapor [m/s].
- `β_ice`: Ice fraction weighting the surface humidity when frozen [-].
- `ψ_sfc`: Soil water potential at the surface [m].
- `Tf_depressed`: Depressed freezing temperature at the surface [K].
- `h_sfc`, `displ`: Surface height and displacement height [m].
- `P_atmos`, `T_atmos`, `q_atmos`, `u_atmos`: Atmospheric pressure [Pa],
  temperature [K], specific humidity [kg/kg], and wind (a speed or a
  horizontal vector) [m/s] at the reference height `atmos_h` [m].
- `roughness_model`, `gustiness`: SurfaceFluxes roughness and gustiness
  specifications.
- `update_∂T_sfc∂T`, `update_∂q_sfc∂T`: Functions for the derivatives of the
  surface temperature and humidity with respect to the soil temperature, as in
  `turbulent_fluxes!`.
- `earth_param_set`: Land parameters.

# Returns
`(; lhf, shf, vapor_flux_liq, vapor_flux_ice, T_sfc)`, with `ρτxz`, `ρτyz`,
and `buoyancy_flux` before `T_sfc` if `return_extra_fluxes` is `Val(true)`: the
latent and sensible heat fluxes [W/m²], the vapor flux [m/s of liquid water]
from liquid water or from ice (sublimation at and below the freezing
temperature), and the skin temperature [K] (see `soil_turbulent_fluxes`).

Called from [`update_soil_surface_temperature!`](@ref).
"""
@inline function solve_soil_surface_temperature_at_a_point(
    return_extra_fluxes::Val,
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
    update_∂T_sfc∂T::UDT,
    update_∂q_sfc∂T::UDQ,
    earth_param_set,
) where {FT, UDT, UDQ}
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    _σ = LP.Stefan(earth_param_set)
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
    q_sfc_guess = soil_specific_humidity(
        T_top,
        ρ_sfc,
        ψ_sfc,
        Tf_depressed,
        earth_param_set,
    )
    output = ClimaLand.surface_fluxes_at_a_point(
        T_top,
        q_sfc_guess,
        update_T,
        update_q,
        P_atmos,
        T_atmos,
        q_atmos,
        u_atmos,
        atmos_h,
        h_sfc,
        displ,
        roughness_model,
        gustiness,
        earth_param_set,
    )
    fluxes = ClimaLand.turbulent_fluxes_from_output(
        return_extra_fluxes,
        output,
        output.T_sfc,
        output.q_vap_sfc,
        update_∂T_sfc∂T,
        update_∂q_sfc∂T,
        P_atmos,
        T_atmos,
        q_atmos,
        atmos_h - h_sfc,
        displ,
        earth_param_set,
    )
    return soil_turbulent_fluxes(fluxes, output.T_sfc, Tf_depressed)
end

"""
    soil_surface_vapor_conductance!(g_soil_sfc, model::EnergyHydrology, Y, p)

Compute the conductance [m/s] of the dry soil layer to water vapor at the soil
surface into `g_soil_sfc` and return it. The liquid water content is
extrapolated to the surface from the top two cells.

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
    ν_sfc = ClimaLand.Domains.top_center_to_surface(ν)
    θ_r_sfc = ClimaLand.Domains.top_center_to_surface(θ_r)
    θ_i_sfc = ClimaLand.Domains.top_center_to_surface(Y.soil.θ_i)
    θ_l_sfc = g_soil_sfc
    ClimaLand.Domains.linear_interpolation_to_surface!(
        θ_l_sfc,
        p.soil.θ_l,
        model.domain.fields.z,
        model.domain.fields.Δz_top,
    )
    # The extrapolation can leave the physical range [θ_r, ν]
    ε = eps(FT)
    @. θ_l_sfc = clamp(θ_l_sfc, θ_r_sfc + ε, ν_sfc)
    S_l_sfc = g_soil_sfc # currently set to θ_l_sfc
    @. S_l_sfc = effective_saturation(ν_sfc, θ_l_sfc, θ_r_sfc) # overwrite with S_l_sfc
    _D_vapor = FT(LP.D_vapor(earth_param_set))
    # currently set to S_l_sfc; overwrite with the conductance
    g_soil_sfc .= soil_conductance.(
        S_l_sfc,
        S_c_sfc,
        d_ds,
        evap_p,
        evap_α,
        _D_vapor,
        ν_sfc,
        θ_r_sfc,
        θ_i_sfc,
    )
    # Reusing g_soil_sfc for the intermediates keeps the kernel argument
    # count within the parameter memory limit of P100 GPUs
    return g_soil_sfc
end

"""
    update_soil_surface_temperature!(model::EnergyHydrology, SW_n, LW_d, Y, p, t;
                                     h_atmos = nothing, u_atmos = p.drivers.u,
                                     T_atmos = p.drivers.T, q_atmos = p.drivers.q,
                                     gustiness = nothing)
    update_soil_surface_temperature!(model::EnergyHydrology, Y, p, t)

Solve for the soil skin temperature from the surface energy balance and store
it, with the turbulent fluxes at it from the same Monin-Obukhov solve, in
`p.soil.turbulent_fluxes`, and update the surface specific humidity
`p.soil.q_sfc` at it; return `nothing`. A no-op unless the top boundary
condition is an `AtmosDrivenFluxBC`.

# Arguments
- `SW_n`: Net shortwave radiation at the soil surface, positive upward, i.e.
  minus the absorbed shortwave radiation [W/m²].
- `LW_d`: Downwelling longwave radiation at the soil surface [W/m²].

`SW_n` and `LW_d` may be fields or lazy broadcasts; land models with a canopy
pass the radiation transmitted and emitted by the canopy. The four-argument
method is for soil exposed to the sky and uses the downwelling radiation in
`p.drivers` and the soil albedo. The atmospheric pressure is read from
`p.drivers`, whether prescribed or supplied by a coupler, and
`p.soil.sfc_scratch` is overwritten. The reference height `h_atmos`, the wind
`u_atmos`, temperature `T_atmos`, and specific humidity `q_atmos` at it, and
the gustiness model default to those of the atmospheric forcing; land models with a canopy pass the sub-canopy
reference height, attenuated wind, canopy-air temperature and humidity, and
zero gustiness (see `Canopy.subcanopy_forcing`).

Called from the `soil_boundary_fluxes!` methods and, in integrated models,
from `lsm_radiant_energy_fluxes!`. See also
[`solve_soil_surface_temperature_at_a_point`](@ref).
"""
function update_soil_surface_temperature!(
    model::EnergyHydrology,
    SW_n,
    LW_d,
    Y,
    p,
    t;
    h_atmos = nothing,
    u_atmos = p.drivers.u,
    T_atmos = p.drivers.T,
    q_atmos = p.drivers.q,
    gustiness = nothing,
)
    bc = model.boundary_conditions.top
    bc isa AtmosDrivenFluxBC || return nothing
    atmos = bc.atmos
    h_atmos = isnothing(h_atmos) ? atmos.h : h_atmos
    gustiness =
        isnothing(gustiness) ?
        SurfaceFluxes.ConstantGustinessSpec(atmos.gustiness) : gustiness
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
    g_soil_sfc =
        soil_surface_vapor_conductance!(p.soil.sfc_scratch, model, Y, p)
    h_sfc = ClimaLand.surface_height(model, Y, p)
    roughness_model = ClimaLand.surface_roughness_model(model, Y, p)
    displ = ClimaLand.surface_displacement_height(model, Y, p)
    r = @. lazy(Δz_top / κ_top)
    β_ice = @. lazy(frozen_soil_vapor_weight(θ_i_sfc, ν_sfc))
    return_extra_fluxes = Val(ClimaLand.return_momentum_fluxes(atmos))
    update_∂T_sfc∂T = ClimaLand.get_∂T_sfc∂T_function(model, Y, p)
    update_∂q_sfc∂T = ClimaLand.get_∂q_sfc∂T_function(model, Y, p)
    p.soil.turbulent_fluxes .= solve_soil_surface_temperature_at_a_point.(
        return_extra_fluxes,
        T_top,
        r,
        ϵ,
        SW_n,
        LW_d,
        g_soil_sfc,
        β_ice,
        ψ_sfc,
        Tf_depressed_sfc,
        h_sfc,
        displ,
        p.drivers.P,
        T_atmos,
        q_atmos,
        u_atmos,
        roughness_model,
        h_atmos,
        gustiness,
        update_∂T_sfc∂T,
        update_∂q_sfc∂T,
        earth_param_set,
    )
    # Updates the cached surface humidity to match the new skin temperature
    ClimaLand.component_specific_humidity(model, Y, p)
    return nothing
end

function update_soil_surface_temperature!(model::EnergyHydrology, Y, p, t)
    α_sfc = ClimaLand.surface_albedo(model, Y, p)
    SW_n = @. lazy(-(1 - α_sfc) * p.drivers.SW_d)
    return update_soil_surface_temperature!(
        model,
        SW_n,
        p.drivers.LW_d,
        Y,
        p,
        t,
    )
end
