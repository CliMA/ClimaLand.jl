"""
    set_eff_land_radiation_properties!(p, earth_param_set)

Sets the values of the effective land broadband albedo p.α_sfc and emissiviity p.ϵ_sfc,
and the corresponding temperature of blackbody p.T_sfc in place, using downwelling
and upwelling SW radiation, and LW upwelling.
"""
function set_eff_land_radiation_properties!(p, earth_param_set)
    # Effective (radiative) land properties
    _σ = LP.Stefan(earth_param_set)
    FT = typeof(_σ)
    @. p.α_sfc = p.SW_u / max(p.drivers.SW_d, eps(FT))
    @. p.ϵ_sfc = 1
    @. p.T_sfc = (p.LW_u / (p.ϵ_sfc * _σ))^(1 / 4)
end


"""
    Canopy.canopy_LW_energy_fluxes!(p::NamedTuple,
                                         s::PrognosticGroundConditions,
                                         canopy,
                                         radiation::PrescribedRadiativeFluxes,
                                         earth_param_set::PSE,
                                         Y::ClimaCore.Fields.FieldVector,
                                         t,
                                        ) where {PSE}

In standalone mode, this function computes and stores the net
long wave radition, in W/m^2,
absorbed by the canopy.

In integrated mode, we have already computed those quantities in
`lsm_radiant_energy_fluxes!`, so this method does nothing additional.

LW net radiation are stored in `p.canopy.radiative_transfer.LW_n`
"""
function Canopy.canopy_LW_energy_fluxes!(
    p::NamedTuple,
    s::PrognosticGroundConditions,
    canopy,
    radiation::AbstractRadiativeDrivers,
    earth_param_set::PSE,
    Y::ClimaCore.Fields.FieldVector,
    t,
) where {PSE}
    nothing
end

"""
    Canopy.ground_albedo_PAR(
        prognostic_land_components::Val{(:canopy, :soil, :soilco2)},
        ground,
        Y,
        p,
        t,
    )

A method of Canopy.Canopy.ground_albedo_PAR for a prognostic soil.
"""
function Canopy.ground_albedo_PAR(
    prognostic_land_components::Val{(:canopy, :soil, :soilco2)},
    ground,
    Y,
    p,
    t,
)
    return p.soil.PAR_albedo
end

"""
    Canopy.ground_albedo_NIR(
        prognostic_land_components::Val{(:canopy, :soil, :soilco2)},
        ground,
        Y,
        p,
        t,
    )

A method of Canopy.ground_albedo_NIR for a prognostic soil.
"""
function Canopy.ground_albedo_NIR(
    prognostic_land_components::Val{(:canopy, :soil, :soilco2)},
    ground,
    Y,
    p,
    t,
)
    return p.soil.NIR_albedo
end


"""
    Canopy.ground_albedo_PAR(
       prognostic_land_components::Union{
            Val{(:canopy, :snow, :soil, :soilco2)},
            Val{(:canopy, :snow, :soil)},
        },
        ground,
        Y,
        p,
        t,
    )

A method of Canopy.ground_albedo_PAR for a prognostic soil/snow. This function is called in
the Canopy update_aux! function.
"""
function Canopy.ground_albedo_PAR(
    prognostic_land_components::Union{
        Val{(:canopy, :snow, :soil, :soilco2)},
        Val{(:canopy, :snow, :soil)},
    },
    ground,
    Y,
    p,
    t,
)
    @. p.α_ground.PAR =
        (1 - p.snow.snow_cover_fraction) * p.soil.PAR_albedo +
        p.snow.snow_cover_fraction * p.snow.α_snow
    return p.α_ground.PAR
end

"""
    Canopy.ground_albedo_NIR(
        prognostic_land_components::Union{
            Val{(:canopy, :snow, :soil, :soilco2)},
            Val{(:canopy, :snow, :soil)},
            },
        ground,
        Y,
        p,
        t,
    )

A method of Canopy.ground_albedo_NIR for a prognostic soil/snow. This function is called in
the Canopy update_aux! function.
"""
function Canopy.ground_albedo_NIR(
    prognostic_land_components::Union{
        Val{(:canopy, :snow, :soil, :soilco2)},
        Val{(:canopy, :snow, :soil)},
    },
    ground,
    Y,
    p,
    t,
)
    @. p.α_ground.NIR =
        (1 - p.snow.snow_cover_fraction) * p.soil.NIR_albedo +
        p.snow.snow_cover_fraction * p.snow.α_snow
    return p.α_ground.NIR
end
