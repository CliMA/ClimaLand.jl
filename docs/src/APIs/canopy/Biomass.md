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
ClimaLand.Canopy.PrescribedAreaIndices{FT}(
    LAI::AbstractTimeVaryingInput,
    SAI::FT,
    RAI::FT,
) where {FT <: AbstractFloat}
ClimaLand.Canopy.prescribed_lai_era5
ClimaLand.Canopy.prescribed_lai_modis
ClimaLand.Canopy.update_biomass!
```
