"""
    EnergyHydrology(FT, domain, earth_param_set, forcing, prognostic_land_components;
                         runoff_model =  ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(f_over = FT(3.28), # extract from EPS
                                                                                  R_sb = FT(1.484e-4 / 1000),# extract from EPS
                                                                                  f_max = topmodel_fmax(domain.space.surface,FT),
                                                                                  ),
                         retention_parameters = soil_vangenuchten_parameters(domain.space.subsurface, FT), # Should this be vanGenuchten type?
                         composition_parameters = soil_composition_parameters(domain.space.subsurface, FT),
                         albedo_parameters = clm_soil_albedo_parameters(domain.space.surface, FT), # eventually, can be CLMSoilAlbedo type
                         S_s = ClimaCore.Fields.zeros(domain.space.subsurface) .+ 1e-3,# extract from EPS or get from file
                         )

Creates a ClimaLand.Soil.EnergyHydrology model with the given float type FT, domain, earth_param_set, forcing, and prognostic land components.

Default spatially varying parameters (for retention curve parameters, composition, albedo, and specific storativity) are provided but can be
changed with keyword arguments. 

The runoff model is also provided and can be changed with a keyword argument.
"""
function EnergyHydrology(
    FT,
    domain,
    earth_param_set,
    forcing,
    prognostic_land_components;
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
    albedo = Soil.CLMTwoBandSoilAlbedo{FT}(; clm_soil_albedo_parameters(domain.space.surface, FT)...),
    S_s = ClimaCore.Fields.zeros(domain.space.subsurface) .+ 1e-3,# extract from EPS or get from file
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
