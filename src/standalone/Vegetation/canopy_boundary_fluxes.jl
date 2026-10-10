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
    # Monin-Obukhov similarity needs the forcing above the roughness sublayer
    if atmos isa PrescribedAtmosphere &&
       turbulent_flux_parameterization isa MoninObukhovCanopyFluxes
        (; displ, z_0m) = turbulent_flux_parameterization
        clearance = minimum(@. atmos.h - displ - z_0m)
        clearance > 0 || throw(
            ArgumentError(
                "The atmospheric reference height `atmos.h` must exceed the canopy displacement height plus the momentum roughness length everywhere; the minimum clearance is $clearance m.",
            ),
        )
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

- `p.canopy.turbulent_fluxes`: Canopy SHF, LHF, transpiration, derivatives of these with respect to T,
  the canopy temperature and effective surface humidity of the flux solve, and its similarity scales
- `p.canopy.hydraulics.fa[end]`: Transpiration
- `p.canopy.hydraulics.fa_roots`: Root water flux
- `p.canopy.radiative_transfer.LW_n`: net long wave radiation
- `p.canopy.radiative_transfer.SW_n`: net short wave radiation
"""
NVTX.@annotate function canopy_boundary_fluxes!(
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
    canopy_turbulent_fluxes!(p, canopy, Y, t)
    canopy_root_fluxes!(p, canopy, Y, t)
end

"""
    canopy_turbulent_fluxes!(p, canopy::CanopyModel, Y, t)

Compute the canopy turbulent fluxes `p.canopy.turbulent_fluxes` (sensible and
latent heat, transpiration, their temperature derivatives, and the canopy-air
state and similarity scales of the solve) with `ClimaLand.turbulent_fluxes!`,
and zero them where plants are absent
([`zero_canopy_fluxes_without_plants!`](@ref)); return `nothing`.

Called from [`canopy_boundary_fluxes!`](@ref) and, before the ground skin
solves that read the canopy-air state, from the integrated models.
"""
function canopy_turbulent_fluxes!(p, canopy::CanopyModel, Y, t)
    bc = canopy.boundary_conditions
    ClimaLand.turbulent_fluxes!(
        p.canopy.turbulent_fluxes,
        bc.atmos,
        canopy,
        Y,
        p,
        t,
    )
    zero_canopy_fluxes_without_plants!(p.canopy.turbulent_fluxes, p)
    return nothing
end

"""
    canopy_root_fluxes!(p, canopy::CanopyModel, Y, t)

Update the root fluxes of water `p.canopy.hydraulics.fa_roots` and of energy
`p.canopy.energy.fa_energy_roots` per unit ground area in place; return
`nothing`. Called from [`canopy_boundary_fluxes!`](@ref) and from the
integrated models.
"""
function canopy_root_fluxes!(p, canopy::CanopyModel, Y, t)
    bc = canopy.boundary_conditions
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
    return nothing
end

"""
    plants_present(LAI, SAI)

Whether the plant area index `LAI + SAI` reaches the threshold 0.05 below
which the canopy exchanges no heat with the air and the ground below it is
forced as a standalone surface. Called from
[`zero_canopy_fluxes_without_plants!`](@ref) and [`subcanopy_forcing`](@ref).
"""
plants_present(LAI, SAI) = LAI + SAI >= 0.05

"""
    zero_canopy_fluxes_without_plants!(turbulent_fluxes, p)

Set the canopy sensible heat flux and its temperature derivative to zero
where plants are absent ([`plants_present`](@ref)), and the
latent heat flux, transpiration, and the latent heat flux derivative to
zero where the leaf area index is below 0.05 or the transpiration is not
positive (the canopy does not condense water); modify `turbulent_fluxes`
in place and return `nothing`. The area indices are read from
`p.canopy.biomass.area_index`.

Below these area indices the Monin-Obukhov solve exchanges heat with a
canopy of negligible conductance, and the fluxes it returns are roundoff.

Called from `canopy_boundary_fluxes!` and from the implicit boundary flux
update of the canopy energy model, so that the explicit and implicit
fluxes agree.
"""
function zero_canopy_fluxes_without_plants!(turbulent_fluxes, p)
    zero_without_plants(X, lai, sai) =
        ifelse(plants_present(lai, sai), X, zero(X))
    zero_below_lai(X, vapor_flux, lai) =
        ifelse(lai < 0.05 || vapor_flux <= 0, zero(X), X)
    area_index = p.canopy.biomass.area_index
    @. turbulent_fluxes.shf = zero_without_plants(
        turbulent_fluxes.shf,
        area_index.leaf,
        area_index.stem,
    )
    @. turbulent_fluxes.∂shf∂T = zero_without_plants(
        turbulent_fluxes.∂shf∂T,
        area_index.leaf,
        area_index.stem,
    )
    @. turbulent_fluxes.∂lhf∂T = zero_below_lai(
        turbulent_fluxes.∂lhf∂T,
        turbulent_fluxes.vapor_flux,
        area_index.leaf,
    )
    @. turbulent_fluxes.lhf = zero_below_lai(
        turbulent_fluxes.lhf,
        turbulent_fluxes.vapor_flux,
        area_index.leaf,
    )
    @. turbulent_fluxes.vapor_flux = zero_below_lai(
        turbulent_fluxes.vapor_flux,
        turbulent_fluxes.vapor_flux,
        area_index.leaf,
    )
    return nothing
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
    T_air = p.drivers.T
    P_air = p.drivers.P
    q_air = p.drivers.q
    # Below we approximate the surface air density with the
    # atmospheric density to make it independent of T
    # This makes our estimate of the derivatives more exact later on
    q_sfc = @. lazy(
        Thermodynamics.q_vap_saturation(
            thermo_params,
            T_sfc,
            Thermodynamics.air_density(thermo_params, T_air, P_air, q_air),
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
    q_canopy = component_specific_humidity(model, Y, p)
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
        q_canopy,
    )
        g_leaf = leaf_Cd * u_star * LAI
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

        # Solve for q_sfc analytically to satisfy balance of fluxes:
        # Flux_aero = ρ * g_h * (q_sfc - q_atm)
        # Flux_stom = ρ * (q_canopy - q_sfc) / r_land
        # Equating fluxes: g_h * (q_sfc - q_atm) = (q_canopy - q_sfc) / r_land
        # q_sfc * (g_h + 1/r_land) = q_canopy/r_land + g_h * q_atm
        # q_sfc = (q_canopy + g_h * r_land * q_atm) / (1 + g_h * r_land)

        q_new = (g_land / g_h * q_canopy + q_vap_int) / (1 + g_land / g_h)
        # Condensation is not modelled, so q_new >= q_vap_int
        return max(q_new, q_vap_int)
    end
    # Closure
    update_q_vap_sfc_field(Cd, LAI, r, qc) =
        (args...) -> update_q_vap_sfc_at_a_point(args..., Cd, LAI, r, qc)
    return @. lazy(update_q_vap_sfc_field(Cd, LAI, r_stomata_canopy, q_canopy))
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
    # Sensible heat is exchanged by all plant surfaces, leaves and stems
    # (plant area index); transpiration (humidity callback) passes through
    # stomata and uses the leaf area only.
    area_index = p.canopy.biomass.area_index
    T_canopy = canopy_temperature(model.energy, model, Y, p)
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
        area_index_pt,
        T_canopy,
    )
        Φ_sfc = SurfaceFluxes.surface_geopotential(param_set, inputs)
        Φ_int = SurfaceFluxes.interior_geopotential(param_set, inputs)
        T_int = inputs.T_int
        g_h = SurfaceFluxes.heat_conductance(
            param_set,
            ζ,
            u_star,
            inputs,
            z_0m,
            z_0b,
            scheme,
        )
        AI = area_index_pt.leaf + area_index_pt.stem
        g_land = leaf_Cd * u_star * AI

        ΔΦ = Φ_int - Φ_sfc
        cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
        T_sfc =
            (T_int + T_canopy * g_land / g_h + ΔΦ / cp_d) / (1 + g_land / g_h)
        return T_sfc
    end
    # Closure
    update_T_sfc_field(Cd, AI, T_c) =
        (args...) -> update_T_sfc_at_a_point(args..., Cd, AI, T_c)
    return @. lazy(update_T_sfc_field(Cd, area_index, T_canopy))
end


"""
    ClimaLand.get_∂q_sfc∂T_function(model::CanopyModel, Y, p)

a helper function which creates and returns the function which computes
the partial derivative of the surface specific humididity with respect to
the canopy temperature.

The derivative is zero where the surface humidity is held at the
atmospheric value by `get_update_surface_humidity_function`, i.e. where
the saturation specific humidity of the canopy `q_sat` does not exceed
the atmospheric specific humidity.
"""
function ClimaLand.get_∂q_sfc∂T_function(model::CanopyModel, Y, p)
    sfp = model.boundary_conditions.turbulent_flux_parameterization
    Cd = sfp.Cd
    LAI = p.canopy.biomass.area_index.leaf
    r_stomata_canopy = p.canopy.conductance.r_stomata_canopy
    q_air = p.drivers.q
    function update_∂q_sfc∂T_at_a_point(
        u_star,
        g_h,
        q_sat,
        T_sfc,
        earth_param_set,
        leaf_Cd,
        LAI,
        r_stomata_canopy,
        q_air,
    )
        g_leaf = leaf_Cd * u_star * LAI
        g_stomata = 1 / r_stomata_canopy
        g_land = g_stomata * g_leaf / (g_leaf + g_stomata)
        ∂q_sfc∂q = (g_land / g_h) / (1 + g_land / g_h)
        ∂q_sfc∂T =
            ∂q_sfc∂q * ClimaLand.partial_q_sat_partial_T(
                q_sat,
                T_sfc,
                Thermodynamics.Liquid(),
                earth_param_set,
            )
        return ifelse(q_sat > q_air, ∂q_sfc∂T, zero(∂q_sfc∂T))
    end
    # Closure
    update_∂q_sfc∂T_field(LAI_val, r_val, leaf_Cd, q_air_val) =
        (args...) -> update_∂q_sfc∂T_at_a_point(
            args...,
            leaf_Cd,
            LAI_val,
            r_val,
            q_air_val,
        )
    return @. lazy(update_∂q_sfc∂T_field(LAI, r_stomata_canopy, Cd, q_air))
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
    area_index = p.canopy.biomass.area_index
    function update_∂T_sfc∂T_at_a_point(
        u_star,
        g_h,
        earth_param_set,
        leaf_Cd,
        area_index_pt,
    )
        AI = area_index_pt.leaf + area_index_pt.stem
        g_land = leaf_Cd * u_star * AI
        ∂T_sfc∂T = (g_land / g_h) / (1 + g_land / g_h)
        return ∂T_sfc∂T
    end
    # Closure
    update_∂T_sfc∂T_field(AI_val, leaf_Cd) =
        (args...) -> update_∂T_sfc∂T_at_a_point(args..., leaf_Cd, AI_val)
    return @. lazy(update_∂T_sfc∂T_field(area_index, Cd))
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
            (
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
            ),
            NTuple{10, FT},
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
            :∂lhf∂T,
            :∂shf∂T,
            :ρτxz,
            :ρτyz,
            :buoyancy_flux,
            :T_sfc,
            :q_sfc,
            :ustar,
            :ζ,
            :Δz_eff,
        ),
        NTuple{13, FT},
    },
)

"""
    subcanopy_forcing(canopy::CanopyModel, p, h_sfc)

Return the NamedTuple `(; h_atmos, u_atmos, T_atmos, q_atmos, gustiness)` of
lazy broadcasts of the reference height [m], the wind [m/s], the air
temperature [K], the specific humidity [kg/kg], and the gustiness model for
the turbulent fluxes of a ground surface (soil or snow) at height `h_sfc`
below the canopy, from the atmospheric forcing in `p.drivers`,
the canopy height and plant area index in `canopy.biomass` and
`p.canopy.biomass.area_index`, the canopy-air state in
`p.canopy.turbulent_fluxes`, and the canopy turbulent flux parameterization;
see [`subcanopy_reference_height`](@ref), [`subcanopy_wind`](@ref), and
[`ground_gustiness`](@ref).

Where plants are present ([`plants_present`](@ref)), the reference height is
the sub-canopy height, the wind the attenuated wind there, and the gustiness the
model of the forcing without its floor, which is folded into that wind;
`T_atmos` and `q_atmos` are the canopy-air temperature
`p.canopy.turbulent_fluxes.T_sfc` and specific humidity
`p.canopy.turbulent_fluxes.q_sfc` at the scalar roughness height `d + z_0b`
from the canopy Monin-Obukhov solve, which the integrated models run first,
so that the ground exchanges heat and water vapor with the canopy air space
rather than directly with the atmosphere above the canopy (Shuttleworth and
Wallace, 1985). Where plants are absent, the ground is forced as a standalone
surface: the forcing height, wind, temperature, humidity, and gustiness model
of `p.drivers` and the atmospheric driver.
"""
function subcanopy_forcing(canopy::CanopyModel, p, h_sfc)
    atmos = canopy.boundary_conditions.atmos
    sfp = canopy.boundary_conditions.turbulent_flux_parameterization
    (; displ, z_0m, subcanopy_min_reference_height, subcanopy_wind_extinction) =
        sfp
    h_forcing = atmos.h
    atmos_gustiness = ClimaLand.gustiness_spec(atmos)
    surface_flux_params =
        LP.surface_fluxes_parameters(ClimaLand.get_earth_param_set(canopy))
    floor =
        SurfaceFluxes.minimum_wind_speed(atmos_gustiness, surface_flux_params)
    height = canopy.biomass.height
    area_index = p.canopy.biomass.area_index
    u = p.drivers.u
    plants = @. lazy(plants_present(area_index.leaf, area_index.stem))
    h_subcanopy = @. lazy(
        subcanopy_reference_height(
            h_sfc,
            displ,
            z_0m,
            subcanopy_min_reference_height,
            h_forcing,
        ),
    )
    h_atmos = @. lazy(ifelse(plants, h_subcanopy, h_forcing))
    u_atmos = @. lazy(
        ifelse(
            plants,
            subcanopy_wind(
                u,
                floor,
                h_forcing,
                h_subcanopy,
                height,
                displ,
                z_0m,
                area_index.leaf + area_index.stem,
                subcanopy_wind_extinction,
            ),
            u,
        ),
    )
    T_ca = p.canopy.turbulent_fluxes.T_sfc
    q_ca = p.canopy.turbulent_fluxes.q_sfc
    T_atmos = @. lazy(ifelse(plants, T_ca, p.drivers.T))
    q_atmos = @. lazy(ifelse(plants, q_ca, p.drivers.q))
    gustiness = ground_gustiness(atmos_gustiness, plants)
    return (; h_atmos, u_atmos, T_atmos, q_atmos, gustiness)
end
