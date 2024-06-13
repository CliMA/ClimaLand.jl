using SurfaceFluxes
using Thermodynamics
using StaticArrays
import SurfaceFluxes.Parameters as SFP

import ClimaLand:
    surface_temperature,
    surface_specific_humidity,
    surface_evaporative_scaling,
    surface_height,
    surface_resistance,
    displacement_height


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
    earth_param_set = model.parameters.earth_param_set
    R = FT(LP.gas_constant(earth_param_set))
    ρ_liq = FT(LP.ρ_cloud_liq(earth_param_set))
    P_air = p.drivers.P
    T_air = p.drivers.T
    leaf_conductance = p.canopy.conductance.gs
    canopy_conductance =
        upscale_leaf_conductance.(
            leaf_conductance,
            p.canopy.hydraulics.area_index.leaf,
            T_air,
            R,
            P_air,
        )

    return 1 ./ canopy_conductance # [s/m]
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
function ClimaLand.surface_height(model::CanopyModel, _...)
    return model.hydraulics.compartment_surfaces[1]
end

"""
    ClimaLand.surface_specific_humidity(model::CanopyModel, Y, p)

A helper function which returns the surface specific humidity for the canopy
model, which is stored in the aux state.
"""
function ClimaLand.surface_specific_humidity(
    model::CanopyModel,
    Y,
    p,
    T_canopy,
    ρ_canopy,
)
    thermo_params =
        LP.thermodynamic_parameters(model.parameters.earth_param_set)
    return Thermodynamics.q_vap_saturation_generic.(
        Ref(thermo_params),
        T_canopy,
        ρ_canopy,
        Ref(Thermodynamics.Liquid()),
    )
end

function make_update_boundary_fluxes(canopy::CanopyModel)
    function update_boundary_fluxes!(p, Y, t)
        canopy_boundary_fluxes!(p, canopy, canopy.radiation, canopy.atmos, Y, t)
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
                            radiation::PrescribedRadiativeFluxes,
                            atmos::PrescribedAtmosphere,
                            Y::ClimaCore.Fields.FieldVector,
                            t,
                            ) where {FT}

Computes the boundary fluxes for the canopy prognostic
equations; updates the specific fields in the auxiliary
state `p` which hold these variables. This function is called
within the explicit tendency of the canopy model.

- `p.canopy.energy.shf`: Canopy SHF
- `p.canopy.energy.lhf`: Canopy LHF
- `p.canopy.hydraulics.fa[end]`: Transpiration
- `p.canopy.conductance.transpiration`: Transpiration (stored twice; to be addressed in a future PR)
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
    radiation::PrescribedRadiativeFluxes,
    atmos::PrescribedAtmosphere,
    Y::ClimaCore.Fields.FieldVector,
    t,
) where {FT}

    root_water_flux = p.canopy.hydraulics.fa_roots
    root_energy_flux = p.canopy.energy.fa_energy_roots
    fa = p.canopy.hydraulics.fa
    LAI = p.canopy.hydraulics.area_index.leaf
    SAI = p.canopy.hydraulics.area_index.stem
    transpiration = p.canopy.conductance.transpiration
    shf = p.canopy.energy.shf
    lhf = p.canopy.energy.lhf
    r_ae = p.canopy.energy.r_ae
    i_end = canopy.hydraulics.n_stem + canopy.hydraulics.n_leaf

    # Compute transpiration, SHF, LHF
    canopy_tf = canopy_turbulent_fluxes(atmos, canopy, Y, p, t)
    transpiration .= canopy_tf.vapor_flux
    shf .= canopy_tf.shf
    lhf .= canopy_tf.lhf
    r_ae .= canopy_tf.r_ae

    # Transpiration is per unit ground area, not leaf area (mult by LAI)
    fa.:($i_end) .= PlantHydraulics.transpiration_per_ground_area(
        canopy.hydraulics.transpiration,
        Y,
        p,
        t,
    )
    # Update the root flux of water per unit ground area in place
    root_water_flux_per_ground_area!(
        root_water_flux,
        canopy.soil_driver,
        canopy.hydraulics,
        Y,
        p,
        t,
    )

    root_energy_flux_per_ground_area!(
        root_energy_flux,
        canopy.soil_driver,
        canopy.energy,
        Y,
        p,
        t,
    )

    canopy_radiant_energy_fluxes!(
        p,
        canopy.soil_driver,
        canopy,
        canopy.radiation,
        canopy.parameters.earth_param_set,
        Y,
        t,
    )

end

"""
    function canopy_turbulent_fluxes(
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
function canopy_turbulent_fluxes(
    atmos::PrescribedAtmosphere,
    model::CanopyModel,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    T_sfc = ClimaLand.surface_temperature(model, Y, p, t)
    ρ_sfc = ClimaLand.surface_air_density(atmos, model, Y, p, t, T_sfc)
    q_sfc = ClimaLand.surface_specific_humidity(model, Y, p, T_sfc, ρ_sfc)
    h_sfc = ClimaLand.surface_height(model, Y, p)
    leaf_r_stomata = ClimaLand.surface_resistance(model, Y, p, t)
    d_sfc = ClimaLand.displacement_height(model, Y, p)
    u_air = p.drivers.u
    h_air = atmos.h
    return canopy_turbulent_fluxes_at_a_point.(
        T_sfc,
        q_sfc,
        ρ_sfc,
        h_sfc,
        leaf_r_stomata,
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
        q_sfc::FT,
        ρ_sfc::FT,
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
"""
function canopy_turbulent_fluxes_at_a_point(
    T_sfc::FT,
    q_sfc::FT,
    ρ_sfc::FT,
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
    ts_sfc = Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_sfc, T_sfc, q_sfc)

    state_sfc = SurfaceFluxes.StateValues(FT(0), SVector{2, FT}(0, 0), ts_sfc)
    state_in = SurfaceFluxes.StateValues(
        h - d_sfc - h_sfc,
        SVector{2, FT}(u, 0),
        ts_in,
    )

    # State containers
    sc = SurfaceFluxes.ValuesOnly(
        state_in,
        state_sfc,
        z_0m,
        z_0b,
        beta = FT(1),
        gustiness = gustiness,
    )
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    conditions = SurfaceFluxes.surface_conditions(
        surface_flux_params,
        sc;
        tol_neutral = SFP.cp_d(surface_flux_params) / 100000,
    )
    _LH_v0::FT = LP.LH_v0(earth_param_set)
    _ρ_liq::FT = LP.ρ_cloud_liq(earth_param_set)

    cp_m::FT = Thermodynamics.cp_m(thermo_params, ts_in)
    T_in::FT = Thermodynamics.air_temperature(thermo_params, ts_in)
    ΔT = T_in - T_sfc
    r_ae::FT = 1 / (conditions.Ch * SurfaceFluxes.windspeed(sc))
    ρ_air::FT = Thermodynamics.air_density(thermo_params, ts_in)
    ustar::FT = conditions.ustar
    r_b::FT = FT(1 / 0.01 * (ustar / 0.04)^(-1 / 2)) # CLM 5, tech note Equation 5.122
    leaf_r_b = r_b / LAI
    canopy_r_b = r_b / (LAI + SAI)
    E0::FT = SurfaceFluxes.evaporation(surface_flux_params, sc, conditions.Ch)
    E = E0 * r_ae / (leaf_r_b + leaf_r_stomata + r_ae) # CLM 5, tech note Equation 5.101, and fig 5.2b, assuming all sunlit, f_wet = 0
    Ẽ = E / _ρ_liq
    H = -ρ_air * cp_m * ΔT / (canopy_r_b + r_ae) # CLM 5, tech note Equation 5.88, setting H_v = H and solving to remove T_s
    LH = _LH_v0 * E
    return (lhf = LH, shf = H, vapor_flux = Ẽ, r_ae = r_ae)
end
