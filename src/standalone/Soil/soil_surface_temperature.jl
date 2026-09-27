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
via `component_temperature`, and the soil receives the heat flux
(T_top - T_sfc)/r.
When the top cell is frozen (contains ice and is below the depressed freezing
temperature), T_sfc is capped at the depressed freezing temperature; the soil
then receives the atmospheric fluxes at that temperature, which exceed the
skin-top conduction, and the excess melts ice in the top cell.
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

Return the Newton update `T + ΔT` of the soil skin temperature [K], with
`ΔT = -f(T)/f'(T)` and
`r f(T) = r (SW_n + LW_n(T) + L(T) + H(T)) + (T - T_top)`, capped at the
depressed freezing temperature `Tf_depressed` when the top cell is frozen
(contains ice, `β_ice > 0`, and `T_top < Tf_depressed`). Used as the
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
    # Keyed on T_top, so trace ice in a cell above the melting point does not
    # pin the skin
    frozen_top = β_ice > 0 && T_top < Tf_depressed
    T_max = frozen_top ? Tf_depressed : oftype(T_sfc, Inf)
    return min(T_sfc + ΔT, T_max)
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
- `earth_param_set`: Land parameters.

# Returns
`(; lhf, shf, vapor_flux_liq, vapor_flux_ice, T_sfc)`, with `ρτxz`, `ρτyz`,
and `buoyancy_flux` before `T_sfc` if `return_extra_fluxes` is `Val(true)`: the
latent and sensible heat fluxes [W/m²], the vapor flux [m/s of liquid water]
from liquid water or from ice (sublimation at and below the freezing
temperature), and the skin temperature [K].

Called from [`update_soil_surface_temperature!`](@ref).
"""
function solve_soil_surface_temperature_at_a_point(
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
    q_sfc_guess = soil_specific_humidity(
        T_top,
        ρ_sfc,
        ψ_sfc,
        Tf_depressed,
        earth_param_set,
    )
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
    # Vapor flux in volume of liquid water
    Ẽ = output.evaporation / LP.ρ_cloud_liq(earth_param_set)
    is_liquid = ClimaLand.heaviside(T_sfc, Tf_depressed)
    fluxes = (;
        lhf = output.lhf,
        shf = output.shf,
        vapor_flux_liq = Ẽ * is_liquid,
        vapor_flux_ice = Ẽ * (1 - is_liquid),
    )
    return skin_solution(
        return_extra_fluxes,
        fluxes,
        output,
        surface_flux_params,
        T_atmos,
        ρ_atmos,
        q_atmos,
        atmos_h - h_sfc,
    )
end

"""
    skin_solution(return_extra_fluxes, fluxes, output, surface_flux_params, T_atmos, ρ_atmos,
                  q_atmos, Δz)

Return the NamedTuple stored in `p.soil.turbulent_fluxes`: the `fluxes` with
the skin temperature `output.T_sfc`, and, for `Val(true)`, the momentum fluxes
and the buoyancy flux of the SurfaceFluxes `output` before it.
"""
skin_solution(::Val{false}, fluxes, output, args...) =
    (; fluxes..., T_sfc = output.T_sfc)
function skin_solution(
    ::Val{true},
    fluxes,
    output,
    surface_flux_params,
    T_atmos,
    ρ_atmos,
    q_atmos,
    Δz,
)
    FT = typeof(output.T_sfc)
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
        T_sfc = output.T_sfc,
    )
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
    update_soil_surface_temperature!(model::EnergyHydrology, SW_n, LW_d, Y, p, t)
    update_soil_surface_temperature!(model::EnergyHydrology, Y, p, t)

Solve for the soil skin temperature from the surface energy balance and store
it, with the turbulent fluxes at it from the same Monin-Obukhov solve, in
`p.soil.turbulent_fluxes`; return `nothing`. A no-op unless the top boundary
condition is an `AtmosDrivenFluxBC`.

# Arguments
- `SW_n`: Net shortwave radiation at the soil surface, positive upward, i.e.
  minus the absorbed shortwave radiation [W/m²].
- `LW_d`: Downwelling longwave radiation at the soil surface [W/m²].

`SW_n` and `LW_d` may be fields or lazy broadcasts; land models with a canopy
pass the radiation transmitted and emitted by the canopy. The four-argument
method is for soil exposed to the sky and uses the downwelling radiation in
`p.drivers` and the soil albedo. The atmospheric state at the reference height
`atmos.h` is read from `p.drivers`, whether prescribed or supplied by a
coupler, and `p.soil.sfc_scratch` is overwritten.

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
    t,
)
    bc = model.boundary_conditions.top
    bc isa AtmosDrivenFluxBC || return nothing
    atmos = bc.atmos
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
    gustiness = SurfaceFluxes.ConstantGustinessSpec(atmos.gustiness)
    r = @. lazy(Δz_top / κ_top)
    return_extra_fluxes = Val(ClimaLand.return_momentum_fluxes(atmos))
    p.soil.turbulent_fluxes .= solve_soil_surface_temperature_at_a_point.(
        return_extra_fluxes,
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
