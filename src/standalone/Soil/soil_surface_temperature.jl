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
(see `Snow.update_surf_temp!`), refined by bisection where those iterations do
not close the balance, and the same solve yields the turbulent fluxes at
T_sfc, which are stored with it in `p.soil.turbulent_fluxes`. The surface
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
[m³/m³]. The flux at given air humidity and surface temperature does not
otherwise depend on the ice present, so this weight makes the sublimation
vanish as the ice does, faster than `θ_i` itself, which keeps an explicit
time step from drawing on ice that is not there; it plays the role of the
dry-soil-layer conductance for an unfrozen surface. The weight is a
heuristic without a published source; the exponent gives strong suppression
at intermediate ice contents.

Called from `update_soil_surface_temperature!` and
`get_update_surface_humidity_function(::EnergyHydrology)`.
"""
frozen_soil_vapor_weight(θ_i, ν) = (θ_i / ν)^4

"""
    soil_evaporation_beta(q_air, qsat, g_liq, g_h, β_ice, frozen)

Return the β factor [-] of the β approach to soil evaporation (Mahfouf and
Noilhan, 1991): the surface specific humidity of the soil is
`q_sfc = β * qsat + (1 - β) * q_air`, so that the vapor flux is
`β g_h (qsat - q_air)`. β is the ice weight `β_ice` over a frozen surface
losing water (`frozen_soil_vapor_weight`), one over a frozen surface gaining
water, and the conductance ratio `g_liq / (g_liq + g_h)` of the dry soil
layer otherwise.
"""
function soil_evaporation_beta(q_air, qsat, g_liq, g_h, β_ice, frozen::Bool)
    β_frozen = ifelse(q_air < qsat, β_ice, one(β_ice))
    β_unfrozen = (g_liq / g_h) / (1 + g_liq / g_h)
    return ifelse(frozen, promote(β_frozen, β_unfrozen)...)
end

"""
    soil_skin_state(T_sfc, inputs, thermo_params, param_set, ψ_sfc, Tf_depressed, earth_param_set) -> NamedTuple

Return `(; ρ_sfc, qsat, ∂qsat∂T, frozen)`: the surface air density
[kg/m³], extrapolated over the effective height with the current surface
humidity guess of the solve, the saturation specific humidity at `T_sfc` adjusted for the soil
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
    ρ_sfc = SurfaceFluxes.surface_density(
        param_set,
        inputs,
        T_sfc,
        inputs.q_vap_sfc_guess,
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
    q_air = SurfaceFluxes.interior_vapor_specific_humidity(inputs)
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
    w = soil_evaporation_beta(q_air, qsat, g_liq, g_h, β_ice, frozen)
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
    q_air = SurfaceFluxes.interior_vapor_specific_humidity(inputs)
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
    w = soil_evaporation_beta(q_air, qsat, g_liq, g_h, β_ice, frozen)
    return w * qsat + (1 - w) * q_air
end

"""
    soil_flux_inputs(return_extra_fluxes::Val)

Return, as a `Val`, the names of the quantities of
`ClimaLand.turbulent_fluxes_at_a_point` from which `soil_turbulent_fluxes`
forms `p.soil.turbulent_fluxes`: the energy and vapor fluxes, the momentum and
buoyancy fluxes if `return_extra_fluxes` is `Val(true)` (coupled runs), the
effective surface humidity `q_sfc`, and the similarity scales `ustar`, `ζ`,
`Δz_eff` of the solve.
"""
soil_flux_inputs(::Val{false}) =
    Val((:lhf, :shf, :vapor_flux, :q_sfc, :ustar, :ζ, :Δz_eff))
soil_flux_inputs(::Val{true}) = Val((
    :lhf,
    :shf,
    :vapor_flux,
    :ρτxz,
    :ρτyz,
    :buoyancy_flux,
    :q_sfc,
    :ustar,
    :ζ,
    :Δz_eff,
))


"""
    soil_turbulent_fluxes(fluxes, T_sfc, Tf_depressed)

Return the NamedTuple stored in `p.soil.turbulent_fluxes` from the NamedTuple
`fluxes` of `ClimaLand.turbulent_fluxes_at_a_point` at the surface temperature
`T_sfc` [K]: the vapor flux is attributed to liquid water evaporating above the
depressed freezing temperature `Tf_depressed` [K] and to ice sublimating at
and below it, the temperature derivatives are dropped, and `T_sfc`, the
effective surface humidity `q_sfc` of the solve (which includes the dry soil
layer resistance, see `soil_evaporation_beta`), and the similarity scales
`ustar`, `ζ`, `Δz_eff` of the solve are appended after the momentum and
buoyancy fluxes, if present.
"""
function soil_turbulent_fluxes(fluxes, T_sfc, Tf_depressed)
    is_liquid = ClimaLand.heaviside(T_sfc, Tf_depressed)
    extra_fluxes = Base.structdiff(
        fluxes,
        NamedTuple{(
            :lhf,
            :shf,
            :vapor_flux,
            :∂lhf∂T,
            :∂shf∂T,
            :T_sfc,
            :q_sfc,
            :ustar,
            :ζ,
            :Δz_eff,
        )},
    )
    return (;
        lhf = fluxes.lhf,
        shf = fluxes.shf,
        vapor_flux_liq = fluxes.vapor_flux * is_liquid,
        vapor_flux_ice = fluxes.vapor_flux * (1 - is_liquid),
        extra_fluxes...,
        T_sfc,
        q_sfc = fluxes.q_sfc,
        ustar = fluxes.ustar,
        ζ = fluxes.ζ,
        Δz_eff = fluxes.Δz_eff,
    )
end

"""
    soil_skin_energy_residual(output, T_top, r, ϵ, σ, SW_n, LW_d)

Return the residual `r F + (T_sfc - T_top)` [K] of the soil skin energy balance
at the skin temperature `T_sfc` of the SurfaceFluxes.jl `output`, with
`F = SW_n - ϵ (LW_d - σ T_sfc⁴) + lhf + shf` [W/m²] the fluxes leaving the skin
(positive upward) and `(T_top - T_sfc)/r` the conduction from the top cell.
The residual increases with `T_sfc`.

Called from [`solve_soil_surface_temperature_at_a_point`](@ref).
"""
@inline function soil_skin_energy_residual(output, T_top, r, ϵ, σ, SW_n, LW_d)
    T_sfc = output.T_sfc
    F = SW_n - ϵ * (LW_d - σ * T_sfc^4) + output.lhf + output.shf
    return r * F + (T_sfc - T_top)
end

"""
    blend_surface_flux_conditions(a, b, weight)

Return the SurfaceFluxes.jl `SurfaceFluxConditions` whose fields are
`(1 - weight)` times those of `a` plus `weight` times those of `b`, converged
where both are.

Called from [`solve_soil_surface_temperature_at_a_point`](@ref).
"""
@inline function blend_surface_flux_conditions(
    a::SurfaceFluxes.SurfaceFluxConditions{FT},
    b::SurfaceFluxes.SurfaceFluxConditions{FT},
    weight::FT,
) where {FT}
    mix(x, y) = (1 - weight) * x + weight * y
    return SurfaceFluxes.SurfaceFluxConditions{FT}(
        mix(a.shf, b.shf),
        mix(a.lhf, b.lhf),
        mix(a.evaporation, b.evaporation),
        mix(a.ρτxz, b.ρτxz),
        mix(a.ρτyz, b.ρτyz),
        mix(a.ustar, b.ustar),
        mix(a.ζ, b.ζ),
        mix(a.Cd, b.Cd),
        mix(a.g_h, b.g_h),
        mix(a.T_sfc, b.T_sfc),
        mix(a.q_vap_sfc, b.q_vap_sfc),
        mix(a.L_MO, b.L_MO),
        mix(a.L_eff, b.L_eff),
        mix(a.ζ_eff, b.ζ_eff),
        a.converged && b.converged,
    )
end

"""
    SKIN_ENERGY_BALANCE_TOLERANCE

Imbalance [W/m²] of the soil skin energy balance above which
[`solve_soil_surface_temperature_at_a_point`](@ref) refines the skin
temperature by bisection.
"""
const SKIN_ENERGY_BALANCE_TOLERANCE = 5

"""
    SKIN_BISECTION_ITERATIONS

Number of bisection steps of the skin temperature refinement in
[`solve_soil_surface_temperature_at_a_point`](@ref); each halves the bracket and
solves the Monin-Obukhov equations at fixed skin temperature.
"""
const SKIN_BISECTION_ITERATIONS = 12

"""
    solve_soil_surface_temperature_at_a_point(return_extra_fluxes, T_top, r, ϵ, SW_n, LW_d,
                                              g_liq, β_ice, ψ_sfc, Tf_depressed, z_sfc, displ,
                                              P_atmos, T_atmos, q_atmos, u_atmos,
                                              roughness_model, atmos_h, gustiness, reference,
                                              update_∂T_sfc∂T, update_∂q_sfc∂T,
                                              earth_param_set)

Solve the soil skin surface energy balance at a point within the Monin-Obukhov
iterations of SurfaceFluxes.jl, starting from `T_top`, and return the skin
temperature with the turbulent fluxes at it.

The Newton updates within the Monin-Obukhov iterations neglect the dependence
of the heat conductance on the skin temperature. At low wind speed over a wet
surface, the conductance jumps between its capped stable value and a much
larger unstable value where the virtual temperatures of the surface and the air
cross, and these updates alternate between the two regimes. Where the energy
balance is then left open by more than `SKIN_ENERGY_BALANCE_TOLERANCE`, and the
skin is not held at the freezing temperature, the skin temperature is refined
by bisection between the skin temperature of the solve, `T_sfc`, and
`T_sfc - R`, which bracket the root of the residual `R`
(`soil_skin_energy_residual`), with the Monin-Obukhov equations solved at each
fixed skin temperature. At a jump of the conductance the residual changes sign
without a root, and the solutions at the two ends of the final bracket are
blended in the proportion that closes the balance.

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
- `z_sfc`, `displ`: Surface elevation and displacement height [m].
- `P_atmos`, `T_atmos`, `q_atmos`, `u_atmos`: Atmospheric pressure [Pa],
  temperature [K], specific humidity [kg/kg], and wind (a speed or a
  horizontal vector) [m/s] at the reference height `atmos_h` [m].
- `roughness_model`, `gustiness`: SurfaceFluxes roughness and gustiness
  specifications.
- `update_∂T_sfc∂T`, `update_∂q_sfc∂T`: Functions for the derivatives of the
  surface temperature and humidity with respect to the soil temperature,
- `earth_param_set`: Land parameters.

# Returns
`(; lhf, shf, vapor_flux_liq, vapor_flux_ice, T_sfc)`, with `ρτxz`, `ρτyz`,
and `buoyancy_flux` before `T_sfc` if `return_extra_fluxes` is `Val(true)`: the
latent and sensible heat fluxes [W/m²], the vapor flux [m/s of liquid water]
from liquid water or from ice (sublimation at and below the freezing
temperature), and the skin temperature in K (see `soil_turbulent_fluxes`).

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
    z_sfc::FT,
    displ::FT,
    P_atmos::FT,
    T_atmos::FT,
    q_atmos::FT,
    u_atmos,
    roughness_model,
    atmos_h::FT,
    gustiness,
    reference,
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
    z_0m = SurfaceFluxes.momentum_roughness(
        roughness_model,
        nothing,
        surface_flux_params,
        nothing,
    )
    Δz_eff = ClimaLand.compute_Δz_eff(atmos_h, z_0m, displ, reference)
    ρ_sfc = ClimaLand.compute_ρ_sfc(
        surface_flux_params,
        T_atmos,
        P_atmos,
        q_atmos,
        Δz_eff,
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
        z_sfc,
        displ,
        roughness_model,
        gustiness,
        reference,
        earth_param_set,
    )
    fluxes_at(T_sfc) = ClimaLand.surface_fluxes_at_a_point(
        T_sfc,
        q_sfc_guess,
        nothing,
        update_q,
        P_atmos,
        T_atmos,
        q_atmos,
        u_atmos,
        atmos_h,
        z_sfc,
        displ,
        roughness_model,
        gustiness,
        reference,
        earth_param_set,
    )
    residual(output) =
        soil_skin_energy_residual(output, T_top, r, ϵ, _σ, SW_n, LW_d)
    R = residual(output)
    # At the freezing cap the energy balance is not closed by design
    frozen_top = β_ice > 0 && T_top < Tf_depressed
    if !frozen_top && abs(R) > r * FT(SKIN_ENERGY_BALANCE_TOLERANCE)
        # The residual increases with the skin temperature, so the root lies
        # between T_sfc and T_sfc - R
        output_b = fluxes_at(output.T_sfc - R)
        R_b = residual(output_b)
        above_root = R > 0
        output_lo = ifelse(above_root, output_b, output)
        output_hi = ifelse(above_root, output, output_b)
        R_lo = ifelse(above_root, R_b, R)
        R_hi = ifelse(above_root, R, R_b)
        for _ in 1:SKIN_BISECTION_ITERATIONS
            output_mid = fluxes_at((output_lo.T_sfc + output_hi.T_sfc) / 2)
            R_mid = residual(output_mid)
            above_root = R_mid > 0
            output_lo = ifelse(above_root, output_lo, output_mid)
            output_hi = ifelse(above_root, output_mid, output_hi)
            R_lo = ifelse(above_root, R_lo, R_mid)
            R_hi = ifelse(above_root, R_mid, R_hi)
        end
        # Where the conductance jumps, no single skin temperature closes the
        # balance; the two sides are combined in the proportion that does
        weight = clamp(R_lo / min(R_lo - R_hi, -eps(FT)), FT(0), FT(1))
        output = blend_surface_flux_conditions(output_lo, output_hi, weight)
    end
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
        atmos_h,
        displ,
        roughness_model,
        reference,
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
                                     h_atmos = atmos.h, u_atmos = p.drivers.u,
                                     T_atmos = p.drivers.T, q_atmos = p.drivers.q,
                                     gustiness = ClimaLand.gustiness_spec(atmos))
    update_soil_surface_temperature!(model::EnergyHydrology, Y, p, t)

Solve for the soil skin temperature from the surface energy balance and store
it, with the turbulent fluxes at it from the same Monin-Obukhov solve, in
`p.soil.turbulent_fluxes`, and update the surface specific humidity
`p.soil.q_sfc` at it; return `nothing`. The top boundary condition of the
soil must be an `AtmosDrivenFluxBC` with the atmospheric driver `atmos`; an
`ArgumentError` is thrown otherwise.

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
the gustiness model default to those of the atmospheric forcing; land models
with a canopy pass the sub-canopy reference height, attenuated wind,
canopy-air temperature and humidity, and gustiness model (see
`Canopy.subcanopy_forcing`).

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
    h_atmos = model.boundary_conditions.top.atmos.h,
    u_atmos = p.drivers.u,
    T_atmos = p.drivers.T,
    q_atmos = p.drivers.q,
    gustiness = ClimaLand.gustiness_spec(model.boundary_conditions.top.atmos),
    reference = model.boundary_conditions.top.atmos.reference,
)
    bc = model.boundary_conditions.top
    bc isa AtmosDrivenFluxBC || throw(
        ArgumentError(
            "The soil skin temperature solve requires an AtmosDrivenFluxBC at the top of the soil",
        ),
    )
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
    z_sfc = ClimaLand.surface_elevation(model, Y, p)
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
        z_sfc,
        displ,
        p.drivers.P,
        T_atmos,
        q_atmos,
        u_atmos,
        roughness_model,
        h_atmos,
        gustiness,
        reference,
        update_∂T_sfc∂T,
        update_∂q_sfc∂T,
        earth_param_set,
    )
    # Updates the cached surface humidity to match the new skin temperature
    ClimaLand.component_specific_humidity(
        model,
        Y,
        p,
        (; h_atmos, T_atmos, q_atmos, reference),
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
