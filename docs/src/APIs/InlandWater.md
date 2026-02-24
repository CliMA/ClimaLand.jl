# InlandWater

```@meta
CurrentModule = ClimaLand.InlandWater
```
## Models and Parameters

```@docs
ClimaLand.InlandWater.AbstractInlandWaterModel
ClimaLand.InlandWater.SlabLakeParameters
ClimaLand.InlandWater.SlabLakeParameters{FT}(; earth_param_set, kwargs...) where {FT}
ClimaLand.InlandWater.SlabLakeModel
ClimaLand.InlandWater.SlabLakeModel{FT}(;
    parameters,
    domain,
    boundary_conditions,
    inland_water_mask,
    Δz_top,
) where {FT}
```

## Boundary Conditions and Fluxes

```@docs
ClimaLand.InlandWater.AtmosDrivenLakeBC
ClimaLand.InlandWater.lake_boundary_fluxes!(
    bc::ClimaLand.InlandWater.AtmosDrivenLakeBC,
    prognostic_land_components,
    model::ClimaLand.InlandWater.SlabLakeModel{FT},
    Y,
    p,
    t,
) where {FT}
```
