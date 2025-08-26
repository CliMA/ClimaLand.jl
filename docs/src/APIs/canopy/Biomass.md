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
    toml_dict::CP.AbstractTOMLDict;
    SAI::FT = FT(0),
    RAI::FT = FT(1),
    rooting_depth = clm_rooting_depth(domain.space.surface),
    ) where {FT <: AbstractFloat}
ClimaLand.Canopy.PrescribedBiomassModel{FT}(; LAI, SAI, RAI, rooting_depth) where {FT}
ClimaLand.Canopy.PrescribedAreaIndices
Climaland.Canopy.PrescribedAreaIndices{FT}(
    LAI::AbstractTimeVaryingInput,
    SAI::FT,
    RAI::FT,
) where {FT <: AbstractFloat}
```

## Biomass callback
```@docs
ClimaLand.Canopy.make_biomass_callback
```
