export AbstractCanopyInterceptionModel, NoInterception, CLM5Interception

"""
    AbstractCanopyInterceptionModel{FT} <: AbstractCanopyComponent{FT}

Abstract type for models of the interception of precipitation by the canopy
(leaves and stems) and the evaporation of the intercepted water.
"""
abstract type AbstractCanopyInterceptionModel{FT} <: AbstractCanopyComponent{FT} end

ClimaLand.name(::AbstractCanopyInterceptionModel) = :interception
Base.broadcastable(m::AbstractCanopyInterceptionModel) = tuple(m)

"""
    NoInterception{FT} <: AbstractCanopyInterceptionModel{FT}

No canopy interception: all precipitation reaches the ground, and the canopy
vapor flux is entirely transpiration.
"""
struct NoInterception{FT} <: AbstractCanopyInterceptionModel{FT} end

"""
    CLM5Interception{FT} <: AbstractCanopyInterceptionModel{FT}

Interception of liquid precipitation by leaves and stems and evaporation of
the intercepted water from the wetted plant area, following CLM5 (Lawrence
et al., 2019, CLM5 technical note, "Canopy Hydrology"):

- intercepted rain  I = α tanh(PAI) P_liq,
- maximum storage   W_max = p_liq PAI,
- drip              D = max(W + I Δt - W_max, 0) / Δt,
- wetted fraction   f_wet = min((W* / W_max)^(2/3), f_wet_max),

where PAI = LAI + SAI is the plant area index, W the canopy liquid water
store (m of liquid water per unit ground area), and W* = W + (I - D) Δt the
store after interception and drip. The wetted plant area
evaporates through the leaf boundary-layer conductance only, in parallel with
transpiration from the dry leaf area, which uses the stomatal and
boundary-layer conductances in series (see `canopy_vapor_conductances`).
Condensation (dew) onto the canopy is added to the store rather than to the
plant water. Wet-canopy evaporation is limited smoothly so that it cannot
exceed the available water within one timestep.

Snowfall currently passes through the canopy unchanged, and the intercepted
water is treated as liquid. The timestep `Δt` must equal the timestep of the
simulation. The internal energy of the intercepted water, ρ_l e_l(T_air)
(I - D), is added to the canopy energy, following the convention used for the
energy of the root water uptake (see `interception_energy_flux`).

$(DocStringExtensions.FIELDS)
"""
struct CLM5Interception{FT} <: AbstractCanopyInterceptionModel{FT}
    "Interception efficiency coefficient α in I = α tanh(PAI) P_liq (unitless)"
    α_liq::FT
    "Maximum liquid water storage per unit plant area (m)"
    p_liq::FT
    "Maximum wetted fraction of the plant area (unitless)"
    f_wet_max::FT
    "Timestep of the simulation (s), used for drip and the evaporation limit"
    Δt::FT
end

"""
    CLM5Interception{FT}(
        toml_dict::CP.ParamDict,
        Δt;
        α_liq = toml_dict["canopy_interception_efficiency"],
        p_liq = toml_dict["canopy_liquid_storage_per_area"],
        f_wet_max = toml_dict["canopy_maximum_wetted_fraction"],
    ) where {FT}

Constructs a `CLM5Interception` model from the parameters in `toml_dict` and
the simulation timestep `Δt` (s).
"""
function CLM5Interception{FT}(
    toml_dict::CP.ParamDict,
    Δt;
    α_liq = toml_dict["canopy_interception_efficiency"],
    p_liq = toml_dict["canopy_liquid_storage_per_area"],
    f_wet_max = toml_dict["canopy_maximum_wetted_fraction"],
) where {FT}
    return CLM5Interception{FT}(FT(α_liq), FT(p_liq), FT(f_wet_max), FT(Δt))
end

ClimaLand.prognostic_vars(::CLM5Interception) = (:W,)
ClimaLand.prognostic_types(::CLM5Interception{FT}) where {FT} = (FT,)
ClimaLand.prognostic_domain_names(::CLM5Interception) = (:surface,)

ClimaLand.auxiliary_vars(::CLM5Interception) = (
    :f_wet,
    :intercepted_liq,
    :drip,
    :throughfall_liq,
    :E_max,
    :transpiration,
    :fluxes,
)
ClimaLand.auxiliary_types(::CLM5Interception{FT}) where {FT} = (
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    NamedTuple{
        (:lhf, :shf, :vapor_flux, :∂lhf∂T, :∂shf∂T, :transpiration),
        Tuple{FT, FT, FT, FT, FT, FT},
    },
)
ClimaLand.auxiliary_domain_names(::CLM5Interception) =
    (:surface, :surface, :surface, :surface, :surface, :surface, :surface)

"""
    update_interception!(p, Y, t, model::AbstractCanopyInterceptionModel, canopy)

Updates the interception rate, drip, throughfall, wetted fraction and the
limit on wet-canopy evaporation in the cache. This must be called after the
area indices are updated and before the ground boundary fluxes, which use the
throughfall.
"""
update_interception!(p, Y, t, model::NoInterception, canopy) = nothing

function update_interception!(p, Y, t, model::CLM5Interception{FT}, canopy) where {FT}
    (; α_liq, p_liq, f_wet_max, Δt) = model
    area_index = p.canopy.biomass.area_index
    PAI = @. lazy(area_index.leaf + area_index.stem)
    W = Y.canopy.interception.W
    ρ_liq = LP.ρ_cloud_liq(canopy.earth_param_set)
    # P_liq is negative (downward); the interception, drip and E_max are
    # positive magnitudes, the throughfall has the same sign as P_liq.
    @. p.canopy.interception.intercepted_liq =
        α_liq * tanh(max(PAI, FT(0))) * max(-p.drivers.P_liq, FT(0))
    @. p.canopy.interception.drip =
        max(
            W + p.canopy.interception.intercepted_liq * Δt -
            p_liq * max(PAI, FT(0)),
            FT(0),
        ) / Δt
    @. p.canopy.interception.throughfall_liq =
        p.drivers.P_liq + p.canopy.interception.intercepted_liq -
        p.canopy.interception.drip
    # As in CLM5, the wetted fraction and the water available for evaporation
    # are evaluated after the store is updated with interception and drip
    W_star = @. lazy(
        max(
            W +
            (p.canopy.interception.intercepted_liq - p.canopy.interception.drip) *
            Δt,
            FT(0),
        ),
    )
    @. p.canopy.interception.f_wet =
        wetted_fraction(W_star, PAI, p_liq, f_wet_max)
    @. p.canopy.interception.E_max = ρ_liq * W_star / Δt
    return nothing
end

"""
    wetted_fraction(W::FT, PAI::FT, p_liq::FT, f_wet_max::FT) where {FT}

Returns the wetted fraction of the plant area, (W / W_max)^(2/3) with
W_max = p_liq PAI, limited to `f_wet_max`. The canopy is treated as dry when
the plant area index is below 0.05, the threshold below which the canopy
turbulent fluxes are neglected.
"""
function wetted_fraction(W::FT, PAI::FT, p_liq::FT, f_wet_max::FT) where {FT}
    W_max = p_liq * PAI
    f = min(
        (max(W, FT(0)) / max(W_max, eps(FT)))^(FT(2) / FT(3)),
        f_wet_max,
    )
    return ifelse(PAI < FT(0.05), FT(0), f)
end

"""
    canopy_vapor_conductances(
        u_star::FT,
        leaf_Cd::FT,
        LAI::FT,
        PAI::FT,
        r_stomata_canopy::FT,
        f_wet::FT,
        E_max::FT,
        q_canopy::FT,
        q_air::FT,
        ρ_air::FT,
        dew_to_storage::Bool,
    ) where {FT}

Returns the conductances (m/s) for the evaporation of intercepted water
`g_wet` and for transpiration `g_tr`, which act in parallel between the
canopy (at saturation specific humidity `q_canopy`) and the canopy surface
air.

The leaf boundary-layer conductance per unit plant area is `leaf_Cd u_star`.
Transpiration occurs from the dry leaf area (1 - f_wet) LAI through the
stomatal and boundary-layer conductances in series; the wetted plant area
f_wet PAI evaporates through the boundary layer only. The wet-canopy
conductance is combined harmonically with the conductance
E_max / (ρ_air (q_canopy - q_air)) at which the evaporation would use up the
available water `E_max` (kg/m²/s) within one timestep, which keeps the store
non-negative without clipping. If `dew_to_storage` is true and the air is
more humid than the canopy (condensation), all of the condensation is
directed to the canopy store through the boundary layer of the whole plant
area, and there is no transpiration.

With `f_wet = 0`, `E_max = 0` and `dew_to_storage = false`, `g_wet = 0` and
`g_tr` is the conductance of the canopy without interception.
"""
function canopy_vapor_conductances(
    u_star::FT,
    leaf_Cd::FT,
    LAI::FT,
    PAI::FT,
    r_stomata_canopy::FT,
    f_wet::FT,
    E_max::FT,
    q_canopy::FT,
    q_air::FT,
    ρ_air::FT,
    dew_to_storage::Bool,
) where {FT}
    g_leaf = leaf_Cd * u_star * LAI
    g_stomata = 1 / r_stomata_canopy
    g_tr = (1 - f_wet) * (g_stomata * g_leaf / (g_leaf + g_stomata))
    g_b_plant = leaf_Cd * u_star * PAI
    g_wet = f_wet * g_b_plant
    Δq = q_canopy - q_air
    g_cap = E_max / (ρ_air * max(Δq, eps(FT)))
    g_wet_limited = g_wet * g_cap / (g_wet + g_cap + eps(FT))
    condensing = dew_to_storage & (Δq < 0)
    g_wet_eff = ifelse(condensing, g_b_plant, g_wet_limited)
    g_tr_eff = ifelse(condensing, zero(FT), g_tr)
    return (g_wet_eff, g_tr_eff)
end

"""
    interception_flux_args(model::AbstractCanopyInterceptionModel, p)

Returns the wetted fraction, the evaporation limit (kg/m²/s), and whether
condensation is routed to the canopy store, as used in the canopy vapor
conductances. Without interception, these are `0`, `0`, and `false`.
"""
interception_flux_args(model::NoInterception{FT}, p) where {FT} =
    (FT(0), FT(0), false)
interception_flux_args(model::CLM5Interception, p) =
    (p.canopy.interception.f_wet, p.canopy.interception.E_max, true)

"""
    canopy_transpiration(model::AbstractCanopyInterceptionModel, p)

Returns the transpiration (m/s of liquid water), the part of the canopy
vapor flux that is drawn from the plant water.
"""
canopy_transpiration(model::NoInterception, p) =
    p.canopy.turbulent_fluxes.vapor_flux
canopy_transpiration(model::CLM5Interception, p) =
    p.canopy.interception.transpiration

"""
    liquid_throughfall(p)

Returns the liquid water flux reaching the ground below the canopy
(negative downward, m/s); see `ClimaLand.liquid_throughfall`.
"""
liquid_throughfall(p) = ClimaLand.liquid_throughfall(p)

"""
    interception_energy_flux(model::AbstractCanopyInterceptionModel, p, canopy)

Returns the internal energy flux (W/m², positive into the canopy) carried by
the liquid water added to the canopy store, ρ_l e_l(T_air) (I - D). The ground
receives the throughfall at the same energy per unit volume, so the total
energy of the precipitation is conserved. Like the energy of the root water
uptake, this energy enters the canopy energy balance, while the heat capacity
of the intercepted water is neglected.
"""
interception_energy_flux(model::NoInterception{FT}, p, canopy) where {FT} =
    FT(0)
function interception_energy_flux(model::CLM5Interception, p, canopy)
    earth_param_set = canopy.earth_param_set
    return @. lazy(
        Soil.volumetric_internal_energy_liq(p.drivers.T, earth_param_set) *
        (p.canopy.interception.intercepted_liq - p.canopy.interception.drip),
    )
end

function ClimaLand.make_compute_exp_tendency(model::CLM5Interception, canopy)
    function compute_exp_tendency!(dY, Y, p, t)
        # dW/dt = I - D - E_wet, where the wet-canopy evaporation E_wet
        # (negative for condensation) is the canopy vapor flux minus transpiration
        @. dY.canopy.interception.W =
            p.canopy.interception.intercepted_liq - p.canopy.interception.drip -
            (
                p.canopy.turbulent_fluxes.vapor_flux -
                p.canopy.interception.transpiration
            )
    end
    return compute_exp_tendency!
end

"""
    ClimaLand.total_liq_water_vol_per_area!(
        surface_field,
        model::AbstractCanopyInterceptionModel,
        canopy,
        Y,
        p,
        t,
    )

Adds the liquid water volume per unit ground area stored on the canopy to
`surface_field`.
"""
ClimaLand.total_liq_water_vol_per_area!(
    surface_field,
    model::NoInterception,
    canopy,
    Y,
    p,
    t,
) = nothing
function ClimaLand.total_liq_water_vol_per_area!(
    surface_field,
    model::CLM5Interception,
    canopy,
    Y,
    p,
    t,
)
    @. surface_field += Y.canopy.interception.W
    return nothing
end

"""
    check_interception_forcing(model::AbstractCanopyInterceptionModel, atmos)

Checks that the interception model is compatible with the atmospheric
forcing. The partitioning of the canopy vapor flux into transpiration and
wet-canopy evaporation is currently computed only with prescribed
atmospheric forcing.
"""
check_interception_forcing(model::NoInterception, atmos) = nothing
function check_interception_forcing(model::CLM5Interception, atmos)
    atmos isa ClimaLand.PrescribedAtmosphere || throw(
        ArgumentError(
            "CLM5Interception currently requires a PrescribedAtmosphere",
        ),
    )
    return nothing
end
