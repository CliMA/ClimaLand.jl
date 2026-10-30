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
ClimaLand.Canopy.area_index_from_max_lai
ClimaLand.Canopy.AreaIndexFromMaxLAI
ClimaLand.Canopy.prescribed_lai_era5
ClimaLand.Canopy.prescribed_lai_modis
ClimaLand.Canopy.prescribed_climatological_lai_modis
ClimaLand.Canopy.update_biomass!
ClimaLand.Canopy.clip_absent_plant_area!
ClimaLand.Canopy.PrognosticCarbonModel
ClimaLand.Canopy.PrognosticCarbonModel{FT}(
    lai_model::AbstractBiomassModel{FT},
    parameters::PrognosticCarbonParameters{FT},
) where {FT}
ClimaLand.Canopy.PrognosticCarbonParameters
ClimaLand.Canopy.PrognosticCarbonParameters(toml_dict::CP.ParamDict)
ClimaLand.Canopy.update_carbon_fluxes!
ClimaLand.Canopy.leaf_allocation
ClimaLand.Canopy.stem_allocation
ClimaLand.Canopy.tau_stem_scale
ClimaLand.Canopy.equilibrium_carbon_pools
ClimaLand.Canopy.equilibrium_sugar
ClimaLand.Canopy.mask_biomass!(p, prognostic_land_components)
ClimaLand.Canopy.mask_biomass!(
    p,
    prognostic_land_components::Union{
        Val{(:canopy, :lake, :snow, :soil, :soilco2)},
        Val{(:canopy, :lake, :snow, :soil)},},
)
```
