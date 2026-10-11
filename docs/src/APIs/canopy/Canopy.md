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
    LAI::Union{AbstractTimeVaryingInput, Nothing},
    toml_dict::CP.ParamDict,
) where {FT}
ClimaLand.Canopy.AbstractCanopyComponent
ClimaLand.Canopy.clm_canopy_height
```

## Canopy Model Boundary Fluxes

```@docs
ClimaLand.Canopy.AbstractCanopyBC
ClimaLand.Canopy.AtmosDrivenCanopyBC
ClimaLand.Canopy.canopy_boundary_fluxes!
ClimaLand.Canopy.canopy_turbulent_fluxes!
ClimaLand.Canopy.canopy_root_fluxes!
ClimaLand.Canopy.zero_canopy_fluxes_without_plants!
ClimaLand.Canopy.plants_present
ClimaLand.Canopy.MoninObukhovCanopyFluxes
ClimaLand.Canopy.subcanopy_wind
ClimaLand.Canopy.subcanopy_reference_height
ClimaLand.Canopy.subcanopy_forcing
ClimaLand.Canopy.ground_gustiness
```

