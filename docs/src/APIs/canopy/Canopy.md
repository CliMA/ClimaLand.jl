# Canopy

```@meta
CurrentModule = ClimaLand.Canopy
```
## Canopy Model and Parameters

```@docs
ClimaLand.Canopy.CanopyModel
ClimaLand.Canopy.CanopyModel{FT}()
ClimaLand.Canopy.CanopyModel{FT}(
    domain::Union{
        ClimaLand.Domains.Point,
        ClimaLand.Domains.Plane,
        ClimaLand.Domains.SphericalSurface,
    },
    forcing::NamedTuple,
    LAI::AbstractTimeVaryingInput,
    toml_dict::CP.ParamDict,
) where {FT}
ClimaLand.Canopy.AbstractCanopyComponent
```

## Canopy Model Boundary Fluxes

```@docs
ClimaLand.Canopy.AbstractCanopyBC
ClimaLand.Canopy.AtmosDrivenCanopyBC
ClimaLand.Canopy.MoninObukhovCanopyFluxes
```

