using SurfaceFluxes
using Thermodynamics
using StaticArrays
import SurfaceFluxes.Parameters as SFP

import ClimaLand: turbulent_fluxes!, AbstractBC, get_earth_param_set

function get_earth_param_set(model::CanopyModel)
    return model.earth_param_set
end

"""
    AbstractCanopyBC <: ClimaLand.AbstractBC

An abstract type for boundary conditions for the canopy model.
"""
abstract type AbstractCanopyBC <: ClimaLand.AbstractBC end
"""
    AtmosDrivenCanopyBC{
        A <: AbstractAtmosphericDrivers,
        B <: AbstractRadiativeDrivers,
        G <: AbstractGroundConditions,
        R <: AbstractCanopyFluxParameterization,
        C::Tuple
    } <: AbstractCanopyBC

A struct used to specify the canopy fluxes, referred
to as "boundary conditions", at the surface and
bottom of the canopy, for water and energy.

These fluxes include turbulent surface fluxes
computed with Monin-Obukhov theory, radiative fluxes,
and root extraction.
$(DocStringExtensions.FIELDS)
"""
struct AtmosDrivenCanopyBC{
    A <: AbstractAtmosphericDrivers,
    B <: AbstractRadiativeDrivers,
    G <: AbstractGroundConditions,
    R <: AbstractCanopyFluxParameterization,
    C <: Tuple,
} <: AbstractCanopyBC
    "The atmospheric conditions driving the model"
    atmos::A
    "The radiative fluxes driving the model"
    radiation::B
    "Ground conditions"
    ground::G
    "Turbulent flux (latent, sensible, vapor, and momentum) parameterization"
    turbulent_flux_parameterization::R
    "Prognostic land components present"
    prognostic_land_components::C
end

"""
    AtmosDrivenCanopyBC(
        atmos,
        radiation,
        ground,
        turbulent_flux_parameterization;
        prognostic_land_components = (:canopy,),
    )

An outer constructor for `AtmosDrivenCanopyBC` which is
intended for use as a default when running canopy
models.

This is also checks the logic that:
- If the `ground` field is Prescribed, :soil should not be a prognostic_land_component
- If the `ground` field is not Prescribed, :soil should be modeled prognostically.
"""
function AtmosDrivenCanopyBC(
    atmos,
    radiation,
    ground,
    turbulent_flux_parameterization;
    prognostic_land_components = (:canopy,),
)
    if typeof(ground) <: PrescribedGroundConditions
        @assert !(:soil ∈ prognostic_land_components)
    else
        @assert :soil ∈ prognostic_land_components
    end

    args = (
        atmos,
        radiation,
        ground,
        turbulent_flux_parameterization,
        prognostic_land_components,
    )
    return AtmosDrivenCanopyBC(args...)
end

function ClimaLand.get_drivers(bc::AtmosDrivenCanopyBC)
    if typeof(bc.ground) <: PrescribedGroundConditions
        return (bc.atmos, bc.radiation, bc.ground)
    else
        return (bc.atmos, bc.radiation)
    end
end


function make_update_boundary_fluxes(canopy::CanopyModel)
    function update_boundary_fluxes!(p, Y, t)
        canopy_boundary_fluxes!(p, canopy, Y, t)
    end
    return update_boundary_fluxes!
end

"""
    canopy_boundary_fluxes!(p::NamedTuple,
                            canopy::CanopyModel,
                            Y::ClimaCore.Fields.FieldVector,
                            t,
                            )

Computes the boundary fluxes for the canopy prognostic
equations; updates the specific fields in the auxiliary
state `p` which hold these variables. This function is called
within the explicit tendency of the canopy model.

- `p.canopy.turbulent_fluxes`: Canopy SHF, LHF, transpiration, derivatives of these with respect to T,q
- `p.canopy.hydraulics.fa[end]`: Transpiration
- `p.canopy.hydraulics.fa_roots`: Root water flux
- `p.canopy.radiative_transfer.LW_n`: net long wave radiation
- `p.canopy.radiative_transfer.SW_n`: net short wave radiation
"""
function canopy_boundary_fluxes!(
    p::NamedTuple,
    canopy::CanopyModel,
    Y::ClimaCore.Fields.FieldVector,
    t,
)
    # Note that in three functions below,
    # we dispatch off of the ground conditions `bc.ground`
    # to handle standalone canopy simulations vs integrated ones

    bc = canopy.boundary_conditions

    # Update the canopy radiation
    canopy_radiant_energy_fluxes!(
        p,
        bc.ground,
        canopy,
        bc.radiation,
        canopy.earth_param_set,
        Y,
        t,
    )

    # Compute transpiration, SHF, LHF
    ClimaLand.turbulent_fluxes!(
        p.canopy.turbulent_fluxes,
        bc.atmos,
        canopy,
        Y,
        p,
        t,
    )
    # Due to roundoff problem when multiplying and dividing by cp_d, set
    # SHF to zero if LAI < 0.01
    clip_flux_zero_lai(shf::FT, lai::FT) where {FT} =
        lai < FT(0.01) ? FT(0) : shf
    @. p.canopy.turbulent_fluxes.shf = clip_flux_zero_lai(
        p.canopy.turbulent_fluxes.shf,
        p.canopy.biomass.area_index.leaf + p.canopy.biomass.area_index.stem,
    )

    # Update the root flux of water per unit ground area in place
    root_water_flux_per_ground_area!(
        p.canopy.hydraulics.fa_roots,
        bc.ground,
        canopy.hydraulics,
        canopy,
        Y,
        p,
        t,
    )
    # Update the root flux of energy per unit ground area in place
    root_energy_flux_per_ground_area!(
        p.canopy.energy.fa_energy_roots,
        bc.ground,
        canopy.energy,
        canopy,
        Y,
        p,
        t,
    )
end

"""
    ClimaLand.component_temperature(model::CanopyModel, Y, p)

a helper function which returns the component temperature for the canopy
model, which is stored in the aux state.
"""
function ClimaLand.component_temperature(model::CanopyModel, Y, p)
    return canopy_temperature(model.energy, model, Y, p)
end

"""
    ClimaLand.component_specific_humidity(model::CanopyModel, Y, p)

a helper function which returns the surface specific humidity for the canopy
model.
"""
function ClimaLand.component_specific_humidity(model::CanopyModel, Y, p)
    earth_param_set = get_earth_param_set(model)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    T_sfc = component_temperature(model, Y, p)
    ρ_sfc = @. lazy(
        ClimaLand.compute_ρ_sfc(thermo_params, p.drivers.thermal_state, T_sfc),
    )
    q_sfc = @. lazy(
        Thermodynamics.q_vap_saturation_generic(
            thermo_params,
            T_sfc,
            ρ_sfc,
            Thermodynamics.Liquid(),
        ),
    )
    return q_sfc
end

"""
    ClimaLand.surface_displacement_height(model::CanopyModel, Y, p)

a helper function which returns the displacement height for the canopy
model.
"""
function ClimaLand.surface_displacement_height(
    model::CanopyModel{FT},
    Y,
    p,
) where {FT}
    sfp = model.boundary_conditions.turbulent_flux_parameterization
    return sfp.displ
end

"""
    ClimaLand.surface_roughness_model(model::CanopyModel, Y, p)

a helper function which returns the surface roughness model for the canopy
model.
"""
function ClimaLand.surface_roughness_model(
    model::CanopyModel{FT},
    Y,
    p,
) where {FT}
    sfp = model.boundary_conditions.turbulent_flux_parameterization
    return @. lazy(
        SurfaceFluxes.ConstantRoughnessParams{FT}(sfp.z_0m, sfp.z_0b),
    )
end

"""
    ClimaLand.get_update_surface_humidity_function(model::CanopyModel, Y, p)

a helper function which computes and returns the function which updates the guess 
for surface specific humidity to the actual value, for the canopy model.
"""
function ClimaLand.get_update_surface_humidity_function(
    model::CanopyModel,
    Y,
    p,
)
    sfp = model.boundary_conditions.turbulent_flux_parameterization
    Cd = sfp.Cd
    LAI = p.canopy.biomass.area_index.leaf
    r_stomata_canopy = p.canopy.conductance.r_stomata_canopy
    function update_q_vap_sfc_at_a_point(
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        T_sfc,
        u_star,
        z_0m,
        z_0b,
        leaf_Cd,
        LAI,
        r_stomata_canopy,
    )
        FT = eltype(param_set)
        g_leaf = leaf_Cd * max(u_star, FT(1)) * LAI
        g_stomata = 1 / r_stomata_canopy
        g_land = g_stomata * g_leaf / (g_leaf + g_stomata)
        g_h = SurfaceFluxes.heat_conductance(
            param_set,
            ζ,
            u_star,
            inputs,
            z_0m,
            z_0b,
            scheme,
        )

        q_vap_int = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
        q_canopy = inputs.q_vap_sfc_guess

        # Solve for q_sfc analytically to satisfy balance of fluxes:
        # Flux_aero = ρ * g_h * (q_sfc - q_atm)
        # Flux_stom = ρ * (q_canopy - q_sfc) / r_land
        # Equating fluxes: g_h * (q_sfc - q_atm) = (q_canopy - q_sfc) / r_land
        # q_sfc * (g_h + 1/r_land) = q_canopy/r_land + g_h * q_atm
        # q_sfc = (q_canopy + g_h * r_land * q_atm) / (1 + g_h * r_land)

        q_new = (g_land / g_h * q_canopy + q_vap_int) / (1 + g_land / g_h)
        return q_new
    end
    # Closure
    update_q_vap_sfc_field(LAI_val, r_val, leaf_Cd) =
        (args...) ->
            update_q_vap_sfc_at_a_point(args..., leaf_Cd, LAI_val, r_val)
    return @. lazy(update_q_vap_sfc_field(LAI, r_stomata_canopy, Cd))
end

"""
    ClimaLand.get_update_surface_temperature_function(model::CanopyModel, Y, p)

a helper function which computes and returns the function which updates the guess 
for surface temperature to the actual value, for the canopy model.
"""
function ClimaLand.get_update_surface_temperature_function(
    model::CanopyModel,
    Y,
    p,
)
    sfp = model.boundary_conditions.turbulent_flux_parameterization
    Cd = sfp.Cd
    AI = @. lazy(
        p.canopy.biomass.area_index.leaf + p.canopy.biomass.area_index.stem,
    )
    function update_T_sfc_at_a_point(
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        u_star,
        z_0m,
        z_0b,
        leaf_Cd,
        AI,
    )
        FT = eltype(param_set)
        Φ_sfc = SurfaceFluxes.surface_geopotential(inputs)
        Φ_int = SurfaceFluxes.interior_geopotential(param_set, inputs)
        T_int = inputs.T_int
        T_canopy = inputs.T_sfc_guess
        g_h = SurfaceFluxes.heat_conductance(
            param_set,
            ζ,
            u_star,
            inputs,
            z_0m,
            z_0b,
            scheme,
        )
        ws = SurfaceFluxes.windspeed(param_set, ζ, u_star, inputs)
        g_land = leaf_Cd * max(u_star, FT(1)) * AI

        ΔΦ = Φ_int - Φ_sfc
        cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
        cp_d = cp_d

        #        T_0 = Thermodynamics.Parameters.T_0(thermo_params)
        #        DSE_int = Thermodynamics.dry_static_energy(thermo_params, T_int, Φ_int)
        #       DSE_sfc = Thermodynamics.dry_static_energy(thermo_params, T_int + ΔΦ / cp_d, Φ_sfc)
        #       ΔDSE = DSE_int - DSE_sfc
        #       @show ΔDSE

        #       DSE_int2 = cp_d*(T_int-T_0) + Φ_int
        #       DSE_sfc2 = cp_d*(T_int+ - T_0)+ (Φ_int - Φ_sfc) + Φ_sfc
        #       @show DSE_int2 - DSE_sfc2
        T_sfc =
            (T_int + T_canopy * g_land / g_h + ΔΦ / cp_d) / (1 + g_land / g_h)
        return FT(T_sfc)
    end
    # Closure
    update_T_sfc_field(AI_val, leaf_Cd) =
        (args...) -> update_T_sfc_at_a_point(args..., leaf_Cd, AI_val)
    return @. lazy(update_T_sfc_field(AI, Cd))
end

"""
    ClimaLand.get_∂q_sfc∂T_function(model::CanopyModel, Y, p)

a helper function which creates and returns the function which computes
the partial derivative of the surface specific humididity with respect to
the canopy temperature.
"""
function ClimaLand.get_∂q_sfc∂T_function(model::CanopyModel, Y, p)
    sfp = model.boundary_conditions.turbulent_flux_parameterization
    Cd = sfp.Cd
    LAI = p.canopy.biomass.area_index.leaf
    r_stomata_canopy = p.canopy.conductance.r_stomata_canopy
    function update_∂q_sfc∂T_at_a_point(
        u_star,
        g_h,
        T_sfc,
        P_sfc,
        earth_param_set,
        leaf_Cd,
        LAI,
        r_stomata_canopy,
    )
        FT = eltype(earth_param_set)
        _T_freeze = LP.T_freeze(earth_param_set)
        u_star_safe = max(u_star, FT(1))
        r_land =
            1 / (leaf_Cd * u_star_safe) / max(LAI, eps(FT)) + r_stomata_canopy
        ∂q_sfc∂q = 1 / (1 + g_h * r_land)
        return ∂q_sfc∂q *
               ClimaLand.partial_q_sat_partial_T_liq(P_sfc, T_sfc - _T_freeze)
    end
    # Closure
    update_∂q_sfc∂T_field(LAI_val, r_val, leaf_Cd) =
        (args...) ->
            update_∂q_sfc∂T_at_a_point(args..., leaf_Cd, LAI_val, r_val)
    return @. lazy(update_∂q_sfc∂T_field(LAI, r_stomata_canopy, Cd))
end

"""
    ClimaLand.get_∂T_sfc∂T_function(model::CanopyModel, Y, p)

a helper function which creates and returns the function which computes
the partial derivative of the surface temperature with respect to
the canopy temperature.
"""
function ClimaLand.get_∂T_sfc∂T_function(model::CanopyModel, Y, p)
    sfp = model.boundary_conditions.turbulent_flux_parameterization
    Cd = sfp.Cd
    AI = @. lazy(
        p.canopy.biomass.area_index.leaf + p.canopy.biomass.area_index.stem,
    )
    function update_∂T_sfc∂T_at_a_point(
        u_star,
        g_h,
        earth_param_set,
        leaf_Cd,
        AI,
    )
        FT = eltype(earth_param_set)
        g_land = leaf_Cd * max(u_star, FT(1)) * AI
        ∂T_sfc∂T = (g_land / g_h) / (1 + g_land / g_h)
        return ∂T_sfc∂T
    end
    # Closure
    update_∂T_sfc∂T_field(AI_val, leaf_Cd) =
        (args...) -> update_∂T_sfc∂T_at_a_point(args..., leaf_Cd, AI_val)
    return @. lazy(update_∂T_sfc∂T_field(AI, Cd))
end

"""
    boundary_vars(bc, ::ClimaLand.TopBoundary)
    boundary_var_domain_names(bc, ::ClimaLand.TopBoundary)
    boundary_var_types(::AbstractCanopyEnergyModel, bc, ::ClimaLand.TopBoundary)

Fallbacks for the boundary conditions methods which add the turbulent
fluxes to the auxiliary variables.
"""
boundary_vars(bc, ::ClimaLand.TopBoundary) = (:turbulent_fluxes,)
boundary_var_domain_names(bc, ::ClimaLand.TopBoundary) = (:surface,)
boundary_var_types(::CanopyModel{FT}, bc, ::ClimaLand.TopBoundary) where {FT} =
    (
        NamedTuple{
            (:lhf, :shf, :vapor_flux, :∂lhf∂T, :∂shf∂T),
            Tuple{FT, FT, FT, FT, FT},
        },
    )

"""
    boundary_var_types(
        ::CanopyModel{FT},
        ::AtmosDrivenCanopyBC{<:CoupledAtmosphere, <:CoupledRadiativeFluxes},
        ::ClimaLand.TopBoundary,
    ) where {FT}

An extension of the `boundary_var_types` method for AtmosDrivenCanopyBC. This
specifies the type of the additional variables.

This method includes additional flux-related properties needed by the atmosphere:
momentum fluxes (`ρτxz`, `ρτyz`) and the buoyancy flux (`buoy_flux`).
These are updated in place when the coupler computes turbulent fluxes,
rather than in `canopy_boundary_fluxes!`.

Note that we currently store these in the land model because the coupler
computes turbulent land/atmosphere fluxes using ClimaLand functions, and
the land model needs to be able to store the fluxes as an intermediary.
Once we compute fluxes entirely within the coupler, we can remove this.
"""
boundary_var_types(
    ::CanopyModel{FT},
    ::AtmosDrivenCanopyBC{<:CoupledAtmosphere, <:CoupledRadiativeFluxes},
    ::ClimaLand.TopBoundary,
) where {FT} = (
    NamedTuple{
        (
            :lhf,
            :shf,
            :vapor_flux,
            :∂lhf∂qT,
            :∂shf∂T,
            :ρτxz,
            :ρτyz,
            :buoyancy_flux,
        ),
        Tuple{FT, FT, FT, FT, FT, FT, FT, FT},
    },
)
