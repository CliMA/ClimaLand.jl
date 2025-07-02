"""
    EnergyHydrology(FT, domain, forcing;
                         prognostic_land_components = (:soil),
                         earth_param_set = LP.LandParameters(FT),
                         albedo = Soil.CLMTwoBandSoilAlbedo{FT}(; clm_soil_albedo_parameters(domain.space.surface, FT)...),
                         runoff =  ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(f_over = FT(3.28),
                                                                            R_sb = FT(1.484e-4 / 1000),
                                                                            f_max = topmodel_fmax(domain.space.surface,FT),
                                                                            ),
                         retention_parameters = soil_vangenuchten_parameters(domain.space.subsurface, FT),
                         composition_parameters = soil_composition_parameters(domain.space.subsurface, FT),
                         S_s = ClimaCore.Fields.zeros(domain.space.subsurface) .+ 1e-3,
                         )

Creates a ClimaLand.Soil.EnergyHydrology model with the given float type FT, domain, earth_param_set, forcing, and prognostic land components.

When running the soil model in standalone mode, `prognostic_land_components = (:soil,)`, while for running integrated land models,
this should be a list of the component models. This value of this argument must be the same across all components in the land model.

Default spatially varying parameters (for retention curve parameters, composition, and specific storativity) are provided but can be
changed with keyword arguments. 

The runoff and albedo parameterizations are also provided and can be changed via keyword argument.

TODO: Move scalar parameters to ClimaParams and obtain from earth_param_set, possibly use types in retention and composition arguments.
"""
function EnergyHydrology(
    FT,
    domain,
    forcing;
    prognostic_land_components = (:soil),
    earth_param_set = LP.LandParameters(FT),
    albedo::ClimaLand.Soil.AbstractSoilAlbedoParameterization = ClimaLand.Soil.CLMTwoBandSoilAlbedo{
        FT,
    }(;
        clm_soil_albedo_parameters(domain.space.surface, FT)...,
    ),
    runoff::ClimaLand.Soil.Runoff.AbstractRunoffModel = ClimaLand.Soil.Runoff.TOPMODELRunoff{
        FT,
    }(
        f_over = FT(3.28), # extract from EPS
        R_sb = FT(1.484e-4 / 1000),# extract from EPS
        f_max = topmodel_fmax(domain.space.surface, FT),
    ),
    retention_parameters = soil_vangenuchten_parameters(
        domain.space.subsurface,
        FT,
    ),
    composition_parameters = soil_composition_parameters(
        domain.space.subsurface,
        FT,
    ),
    S_s = ClimaCore.Fields.zeros(domain.space.subsurface) .+ 1e-3,
)
    top_bc = ClimaLand.Soil.AtmosDrivenFluxBC(
        forcing.atmos,
        forcing.radiation,
        runoff,
        prognostic_land_components,
    )
    bottom_bc = ClimaLand.Soil.EnergyWaterFreeDrainage()
    boundary_conditions = (; top = top_bc, bottom = bottom_bc)
    # sublimation and subsurface runoff are added automatically
    if :canopy ∈ prognostic_land_components
        sources =
            (ClimaLand.RootExtraction{FT}(), ClimaLand.Soil.PhaseChange{FT}())
    else
        sources = (ClimaLand.Soil.PhaseChange{FT}(),)
    end
    parameters = ClimaLand.Soil.EnergyHydrologyParameters(
        FT;
        retention_parameters...,
        composition_parameters...,
        albedo,
        S_s,
    )
    return ClimaLand.Soil.EnergyHydrology{FT}(;
        parameters,
        domain,
        boundary_conditions,
        sources,
    )
end
