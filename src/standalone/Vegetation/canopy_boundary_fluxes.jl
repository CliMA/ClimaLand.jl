using SurfaceFluxes
using Thermodynamics
using StaticArrays
import SurfaceFluxes.Parameters as SFP

import ClimaLand:
    surface_temperature,
    surface_evaporative_scaling,
    surface_height,
    surface_resistance,
    displacement_height,
    turbulent_fluxes,
    AbstractBC

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
        C::Tuple
    } <: AbstractCanopyBC

A struct used to specify the canopy fluxes, referred
to as ``boundary conditions", at the surface and
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
    C <: Tuple,
} <: AbstractCanopyBC
    "The atmospheric conditions driving the model"
    atmos::A
    "The radiative fluxes driving the model"
    radiation::B
    "Ground conditions"
    ground::G
    "Prognostic land components present"
    prognostic_land_components::C
end

"""
    AtmosDrivenCanopyBC(
        atmos,
        radiation,
        ground;
        prognostic_land_components = (:canopy,),
    )

An outer constructor for `AtmosDrivenCanopyBC` which is
intended for use as a default when running standlone canopy
models. 

This is also checks the logic that:
- If the `ground` field is Prescribed, :soil should not be a prognostic_land_component
- If the `ground` field is not Prescribed, :soil should be modeled prognostically.
"""
function AtmosDrivenCanopyBC(
    atmos,
    radiation,
    ground;
    prognostic_land_components = (:canopy,),
)
    if typeof(ground) <: PrescribedGroundConditions
        @assert !(:soil ∈ prognostic_land_components)
    else
        @assert :soil ∈ prognostic_land_components
    end

    args = (atmos, radiation, ground, prognostic_land_components)
    return AtmosDrivenCanopyBC(args...)
end

"""
    ClimaLand.displacment_height(model::CanopyModel, Y, p)

A helper function which returns the displacement height for the canopy
model.

See Cowan 1968; Brutsaert 1982, pp. 113–116; Campbell and Norman 1998, p. 71; Shuttleworth 2012, p. 343; Monteith and Unsworth 2013, p. 304.
"""
function ClimaLand.displacement_height(model::CanopyModel{FT}, Y, p) where {FT}
    return FT(0.67) * model.hydraulics.compartment_surfaces[end]
end

"""
    ClimaLand.surface_resistance(
        model::CanopyModel{FT},
        Y,
        p,
        t,
    ) where {FT}
Returns the stomatal resistance field of the
`CanopyModel` canopy.
"""
function ClimaLand.surface_resistance(
    model::CanopyModel{FT},
    Y,
    p,
    t,
) where {FT}
    return p.canopy.conductance.r_stomata_canopy
end

"""
    ClimaLand.surface_temperature(model::CanopyModel, Y, p, t)

A helper function which returns the temperature for the canopy
model.
"""
function ClimaLand.surface_temperature(model::CanopyModel, Y, p, t)
    return canopy_temperature(model.energy, model, Y, p, t)
end

"""
    ClimaLand.surface_height(model::CanopyModel, Y, _...)

A helper function which returns the surface height for the canopy
model, which is stored in the parameter struct.
"""
function ClimaLand.surface_height(model::CanopyModel{FT}, _...) where {FT}
    return FT(0)#model.domain.fields.z_sfc FIX!
end

function make_update_boundary_fluxes(canopy::CanopyModel)
    function update_boundary_fluxes!(p, Y, t)
        canopy_boundary_fluxes!(p, canopy, Y, t)
    end
    return update_boundary_fluxes!
end

"""
    canopy_boundary_fluxes!(p::NamedTuple,
                            canopy::CanopyModel{
                                FT,
                                <:AutotrophicRespirationModel,
                                <:Union{BeerLambertModel, TwoStreamModel},
                                <:Union{FarquharModel,OptimalityFarquharModel},
                                <:MedlynConductanceModel,
                                <:PlantHydraulicsModel,
                                <:Union{PrescribedCanopyTempModel,BigLeafEnergyModel}
                            },
                            Y::ClimaCore.Fields.FieldVector,
                            t,
                            ) where {FT}

Computes the boundary fluxes for the canopy prognostic
equations; updates the specific fields in the auxiliary
state `p` which hold these variables. This function is called
within the explicit tendency of the canopy model.

- `p.canopy.energy.turbulent_fluxes`: Canopy SHF, LHF, transpiration, derivatives of these with respect to T,q
- `p.canopy.hydraulics.fa[end]`: Transpiration
- `p.canopy.hydraulics.fa_roots`: Root water flux
- `p.canopy.radiative_transfer.LW_n`: net long wave radiation
- `p.canopy.radiative_transfer.SW_n`: net short wave radiation
"""
function canopy_boundary_fluxes!(
    p::NamedTuple,
    canopy::CanopyModel{
        FT,
        <:AutotrophicRespirationModel,
        <:Union{BeerLambertModel, TwoStreamModel},
        <:Union{FarquharModel, OptimalityFarquharModel},
        <:MedlynConductanceModel,
        <:PlantHydraulicsModel,
        <:Union{PrescribedCanopyTempModel, BigLeafEnergyModel},
    },
    Y::ClimaCore.Fields.FieldVector,
    t,
) where {FT}
    bc = canopy.boundary_conditions
    radiation = bc.radiation
    atmos = bc.atmos
    root_water_flux = p.canopy.hydraulics.fa_roots
    root_energy_flux = p.canopy.energy.fa_energy_roots
    fa = p.canopy.hydraulics.fa
    LAI = p.canopy.hydraulics.area_index.leaf
    SAI = p.canopy.hydraulics.area_index.stem
    canopy_tf = p.canopy.energy.turbulent_fluxes
    i_end = canopy.hydraulics.n_stem + canopy.hydraulics.n_leaf
    # Compute transpiration, SHF, LHF
    canopy_tf .= ClimaLand.turbulent_fluxes(atmos, canopy, Y, p, t)
    # Transpiration is per unit ground area, not leaf area (mult by LAI)
    fa.:($i_end) .= PlantHydraulics.transpiration_per_ground_area(
        canopy.hydraulics.transpiration,
        Y,
        p,
        t,
    )
    # Note that in the three functions below,
    # we dispatch off of the ground conditions `bc.ground`
    # to handle standalone canopy simulations vs integrated ones

    # Update the root flux of water per unit ground area in place
    root_water_flux_per_ground_area!(
        root_water_flux,
        bc.ground,
        canopy.hydraulics,
        Y,
        p,
        t,
    )
    # Update the root flux of energy per unit ground area in place
    root_energy_flux_per_ground_area!(
        root_energy_flux,
        bc.ground,
        canopy.energy,
        Y,
        p,
        t,
    )

    # Update the canopy radiation
    canopy_radiant_energy_fluxes!(
        p,
        bc.ground,
        canopy,
        bc.radiation,
        canopy.parameters.earth_param_set,
        Y,
        t,
    )

end

"""
    function ClimaLand.turbulent_fluxes(
        atmos::PrescribedAtmosphere,
        model::CanopyModel,
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t,
    )

A canopy specific function for compute turbulent fluxes with the atmosphere;
returns the latent heat flux, sensible heat flux, vapor flux, and aerodynamic resistance.

 We cannot use the default version in src/shared_utilities/drivers.jl
because the canopy requires a different resistance for vapor and sensible heat
fluxes, and the resistances depend on ustar, which we must compute using
SurfaceFluxes before adjusting to account for these resistances.
"""
function ClimaLand.turbulent_fluxes(
    atmos::PrescribedAtmosphere,
    model::CanopyModel,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    T_sfc = ClimaLand.surface_temperature(model, Y, p, t)
    h_sfc = ClimaLand.surface_height(model, Y, p)
    r_stomata_canopy = ClimaLand.surface_resistance(model, Y, p, t)
    d_sfc = ClimaLand.displacement_height(model, Y, p)
    u_air = p.drivers.u
    h_air = atmos.h
    return canopy_turbulent_fluxes_at_a_point.(
        T_sfc,
        h_sfc,
        r_stomata_canopy,
        d_sfc,
        p.drivers.thermal_state,
        u_air,
        h_air,
        p.canopy.hydraulics.area_index.leaf,
        p.canopy.hydraulics.area_index.stem,
        atmos.gustiness,
        model.parameters.z_0m,
        model.parameters.z_0b,
        Ref(model.parameters.earth_param_set),
    )
end

"""
    function canopy_turbulent_fluxes_at_a_point(
        T_sfc::FT,
        h_sfc::FT,
        leaf_r_stomata::FT,
        d_sfc::FT,
        ts_in,
        u::FT,
        h::FT,
        LAI::FT,
        SAI::FT,
        gustiness::FT,
        z_0m::FT,
        z_0b::FT,
        earth_param_set::EP,
    ) where {FT <: AbstractFloat, EP}

Computes the turbulent surface fluxes for the canopy at a point
and returns the fluxes in a named tuple.

Note that an additiontal resistance is used in computing both
evaporation and sensible heat flux, and this modifies the output
of `SurfaceFluxes.surface_conditions`.
"""
function canopy_turbulent_fluxes_at_a_point(
    T_sfc::FT,
    h_sfc::FT,
    leaf_r_stomata::FT,
    d_sfc::FT,
    ts_in,
    u::FT,
    h::FT,
    LAI::FT,
    SAI::FT,
    gustiness::FT,
    z_0m::FT,
    z_0b::FT,
    earth_param_set::EP,
) where {FT <: AbstractFloat, EP}
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    # The following will not run on GPU
    #    h - d_sfc - h_sfc < 0 &&
    #        @error("Surface height is larger than atmos height in surface fluxes")
    state_in = SurfaceFluxes.StateValues(
        h - d_sfc - h_sfc,
        SVector{2, FT}(u, 0),
        ts_in,
    )

    ρ_sfc = compute_ρ_sfc(thermo_params, ts_in, T_sfc)
    q_sfc = Thermodynamics.q_vap_saturation_generic(
        thermo_params,
        T_sfc,
        ρ_sfc,
        Thermodynamics.Liquid(),
    )
    ts_sfc = Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_sfc, T_sfc, q_sfc)

    state_sfc = SurfaceFluxes.StateValues(FT(0), SVector{2, FT}(0, 0), ts_sfc)

    # State containers
    states = SurfaceFluxes.ValuesOnly(
        state_in,
        state_sfc,
        z_0m,
        z_0b,
        beta = FT(1),
        gustiness = gustiness,
    )
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    scheme = SurfaceFluxes.PointValueScheme()
    conditions =
        SurfaceFluxes.surface_conditions(surface_flux_params, states, scheme)
    _LH_v0::FT = LP.LH_v0(earth_param_set)
    _ρ_liq::FT = LP.ρ_cloud_liq(earth_param_set)
    cp_m_sfc::FT = Thermodynamics.cp_m(thermo_params, ts_sfc)
    r_ae::FT = 1 / (conditions.Ch * SurfaceFluxes.windspeed(states))
    ustar::FT = conditions.ustar
    ρ_sfc = Thermodynamics.air_density(thermo_params, ts_sfc)
    T_int = Thermodynamics.air_temperature(thermo_params, ts_in)
    Rm_int = Thermodynamics.gas_constant_air(thermo_params, ts_in)
    ρ_air = Thermodynamics.air_density(thermo_params, ts_in)

    r_b::FT = FT(1 / 0.01 * (ustar / 0.04)^(-1 / 2)) # CLM 5, tech note Equation 5.122
    leaf_r_b = r_b / LAI
    canopy_r_b = r_b / (LAI + SAI)

    E0::FT =
        SurfaceFluxes.evaporation(surface_flux_params, states, conditions.Ch)
    E = E0 * r_ae / (leaf_r_b + leaf_r_stomata + r_ae) # CLM 5, tech note Equation 5.101, and fig 5.2b, assuming all sunlit, f_wet = 0
    Ẽ = E / _ρ_liq

    SH =
        SurfaceFluxes.sensible_heat_flux(
            surface_flux_params,
            conditions.Ch,
            states,
            scheme,
        ) * r_ae / (canopy_r_b + r_ae)
    # The above follows from CLM 5, tech note Equation 5.88, setting H_v = SH and solving to remove T_s, ignoring difference between cp in atmos and above canopy. 
    LH = _LH_v0 * E

    # Derivatives
    # We ignore ∂r_ae/∂T_sfc, ∂u*/∂T_sfc, ∂r_stomata∂Tc
    ∂ρsfc∂Tc =
        ρ_air *
        (Thermodynamics.cv_m(thermo_params, ts_in) / Rm_int) *
        (T_sfc / T_int)^(Thermodynamics.cv_m(thermo_params, ts_in) / Rm_int - 1) /
        T_int
    ∂cp_m_sfc∂Tc = 0 # Possibly can address at a later date

    ∂LHF∂qc =
        ρ_sfc * _LH_v0 / (leaf_r_b + leaf_r_stomata + r_ae) +
        LH / ρ_sfc * ∂ρsfc∂Tc

    ∂SHF∂Tc =
        ρ_sfc * cp_m_sfc / (canopy_r_b + r_ae) +
        SH / ρ_sfc * ∂ρsfc∂Tc +
        SH / cp_m_sfc * ∂cp_m_sfc∂Tc
    return (
        lhf = LH,
        shf = SH,
        transpiration = Ẽ,
        r_ae = r_ae,
        ∂LHF∂qc = ∂LHF∂qc,
        ∂SHF∂Tc = ∂SHF∂Tc,
    )

end
