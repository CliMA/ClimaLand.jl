import ClimaLSM:
    surface_temperature,
    surface_specific_humidity,
    surface_air_density,
    surface_evaporative_scaling,
    surface_height,
    surface_resistance

export canopy_turbulent_surface_fluxes

"""
    canopy_turbulent_surface_fluxes(atmos::PrescribedAtmosphere{FT},
                          model::CanopyModel,
                          Y,
                          p,
                          t::FT) where {FT}

Computes canopy transpiration using Monin-Obukhov Surface Theory,
the prescribed atmospheric conditions, and the canopy conductance.

Please note that in the future the SurfaceFluxes.jl code will compute
fluxes taking into account the canopy conductance, so that
what is returned by `surface_fluxes` is correct. At present, it does not,
so we are adjusting for it after the fact here in both ET, LHF, and SHF.
"""
function canopy_turbulent_surface_fluxes(
    atmos::PrescribedAtmosphere{FT},
    model::CanopyModel,
    Y,
    p,
    t::FT,
) where {FT}
    conditions = surface_fluxes(atmos, model, Y, p, t) # per unit m^2 of leaf
    return conditions.vapor_flux, conditions.shf, conditions.lhf
end

"""
    ClimaLSM.surface_resistance(
        model::CanopyModel{FT},
        Y,
        p,
        t,
    ) where {FT}
Returns the surface resistance field of the
`CanopyModel` canopy.
"""
function ClimaLSM.surface_resistance(model::CanopyModel{FT}, Y, p, t) where {FT}
    earth_param_set = model.parameters.earth_param_set
    R = FT(LSMP.gas_constant(earth_param_set))
    ρ_liq = FT(LSMP.ρ_cloud_liq(earth_param_set))
    P_air::FT = model.atmos.P(t)
    T_air::FT = model.atmos.T(t)
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
    ClimaLSM.surface_temperature(model::CanopyModel, Y, p, t)

A helper function which returns the temperature for the canopy
model.
"""
function ClimaLSM.surface_temperature(model::CanopyModel, Y, p, t)
    return canopy_temperature(model.energy, model, Y, p, t)
end

"""
    ClimaLSM.surface_height(model::CanopyModel, Y, _...)

A helper function which returns the surface height for the canopy
model, which is stored in the parameter struct.
"""
function ClimaLSM.surface_height(model::CanopyModel, _...)
    return model.hydraulics.compartment_surfaces[end]
end

"""
    ClimaLSM.surface_specific_humidity(model::CanopyModel, Y, p)

A helper function which returns the surface specific humidity for the canopy
model, which is stored in the aux state.
"""
function ClimaLSM.surface_specific_humidity(
    model::CanopyModel,
    Y,
    p,
    T_canopy,
    ρ_canopy,
)
    thermo_params =
        LSMP.thermodynamic_parameters(model.parameters.earth_param_set)
    return Thermodynamics.q_vap_saturation_generic.(
        Ref(thermo_params),
        T_canopy,
        ρ_canopy,
        Ref(Thermodynamics.Liquid()),
    )
end

"""
    ClimaLSM.surface_air_density(model::CanopyModel, Y, p)

A helper function which computes and returns the surface air density for the canopy
model.
"""
function ClimaLSM.surface_air_density(
    atmos::PrescribedAtmosphere,
    model::CanopyModel,
    Y,
    p,
    t,
    T_canopy,
)
    thermo_params =
        LSMP.thermodynamic_parameters(model.parameters.earth_param_set)
    ts_in = construct_atmos_ts(atmos, t, thermo_params)
    return compute_ρ_sfc.(Ref(thermo_params), Ref(ts_in), T_canopy)
end


"""

"""
function canopy_boundary_fluxes!(
    p::NamedTuple,
    canopy::CanopyModel{
        FT,
        <:AutotrophicRespirationModel,
        <:Union{BeerLambertModel, TwoStreamModel},
        <:FarquharModel,
        <:MedlynConductanceModel,
        <:PlantHydraulicsModel,
        <:PrescribedCanopyTempModel,
    },
    radiation::PrescribedRadiativeFluxes,
    atmos::PrescribedAtmosphere,
    Y::ClimaCore.Fields.FieldVector,
    t::FT,
) where {FT}

    root_water_flux = p.canopy.hydraulics.fa_roots
    fa = p.canopy.hydraulics.fa
    transpiration = p.canopy.conductance.transpiration
    shf = p.canopy.energy.shf
    lhf = p.canopy.energy.lhf
    i_end = canopy.hydraulics.n_stem + canopy.hydraulics.n_leaf

    # Compute transpiration
    (canopy_transpiration, canopy_shf, canopy_lhf) =
        canopy_turbulent_surface_fluxes(atmos, canopy, Y, p, t)
    transpiration .= canopy_transpiration
    shf .= canopy_shf
    lhf .= canopy_lhf

    # Transpiration is per unit ground area, not leaf area (mult by LAI
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

"""
function canopy_boundary_fluxes!(
    p::NamedTuple,
    canopy::CanopyModel{
        FT,
        <: AutotrophicRespirationModel,
        <:Union{BeerLambertModel, TwoStreamModel},
        <:FarquharModel,
        <:MedlynConductanceModel,
        <:PlantHydraulicsModel,
        <:BigLeafEnergyModel,
    },
    radiation::PrescribedRadiativeFluxes,
    atmos::PrescribedAtmosphere,
    Y::ClimaCore.Fields.FieldVector,
    t::FT,
) where {FT}

    # First we need to compute/update in place the boundary fluxes for
    # the component models.
    root_water_flux = p.canopy.hydraulics.fa_roots
    fa = p.canopy.hydraulics.fa
    transpiration = p.canopy.conductance.transpiration
    shf = p.canopy.energy.shf
    lhf = p.canopy.energy.lhf
    root_energy_flux = p.canopy.energy.fa_energy_roots
    i_end = canopy.hydraulics.n_stem + canopy.hydraulics.n_leaf

    # Compute transpiration
    (canopy_transpiration, canopy_shf, canopy_lhf) =
        canopy_turbulent_surface_fluxes(atmos, canopy, Y, p, t)
    transpiration .= canopy_transpiration
    shf .= canopy_shf
    lhf .= canopy_lhf

    # Transpiration is per unit ground area, not leaf area (mult by LAI
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
    # update the net SW and LW radiation absorbed by the canopy, in W/m^2
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
