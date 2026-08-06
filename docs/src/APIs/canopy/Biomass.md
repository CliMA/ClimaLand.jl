# Biomass

```@meta
CurrentModule = ClimaLand.Canopy
```
## Parameterizations

```@docs
ClimaLand.Canopy.PrescribedBiomassModel
ClimaLand.Canopy.PrescribedBiomassModel{FT}(
    domain,
    LAI::AbstractTimeVaryingInput,
    toml_dict::CP.ParamDict;
    SAI::FT = toml_dict["SAI"],
    RAI::FT = toml_dict["RAI"],
    rooting_depth = clm_rooting_depth(domain.space.surface),
    height = toml_dict["canopy_height"]
) where {FT <: AbstractFloat}
ClimaLand.Canopy.PrescribedBiomassModel{FT}(; LAI, SAI::FT, RAI::FT, rooting_depth, height::FT) where {FT}
ClimaLand.Canopy.PrescribedAreaIndices
ClimaLand.Canopy.PrescribedAreaIndices(
    LAI::AbstractTimeVaryingInput,
    SAI,
    RAI,
)
ClimaLand.Canopy.prescribed_lai_era5
ClimaLand.Canopy.prescribed_lai_modis
ClimaLand.Canopy.prescribed_climatological_lai_modis
ClimaLand.Canopy.update_biomass!
ClimaLand.Canopy.PrognosticCarbonModel
ClimaLand.Canopy.PrognosticCarbonModel{FT}(
    lai_model::AbstractBiomassModel{FT},
    parameters::PrognosticCarbonParameters{FT},
) where {FT}
ClimaLand.Canopy.PrognosticCarbonModel{FT}(
    lai_model::AbstractBiomassModel{FT},
    toml_dict::CP.ParamDict;
    kwargs...,
) where {FT}
ClimaLand.Canopy.PrognosticCarbonParameters
ClimaLand.Canopy.PrognosticCarbonParameters(toml_dict::CP.ParamDict)
ClimaLand.Canopy.update_carbon_fluxes!
ClimaLand.Canopy.woody_fraction
ClimaLand.Canopy.seasonality_limit
ClimaLand.Canopy.monthly_pet
ClimaLand.Canopy.tau_stem_scale
ClimaLand.Canopy.mask_biomass!(p, prognostic_land_components)
ClimaLand.Canopy.mask_biomass!(
    p,
    prognostic_land_components::Union{
        Val{(:canopy, :lake, :snow, :soil, :soilco2)},
        Val{(:canopy, :lake, :snow, :soil)},},
)
```
