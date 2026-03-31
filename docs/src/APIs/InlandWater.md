# InlandWater

```@meta
CurrentModule = ClimaLand.InlandWater
```
## Models and Parameters

`SlabLakeParameters` contains the tunable lake physics parameters.

```@docs
ClimaLand.InlandWater.AbstractInlandWaterModel
ClimaLand.InlandWater.SlabLakeParameters
ClimaLand.InlandWater.SlabLakeParameters(toml_dict::CP.ParamDict; kwargs...)
ClimaLand.InlandWater.SlabLakeModel
ClimaLand.InlandWater.SlabLakeModel{FT}(;
    parameters,
    domain,
    boundary_conditions,
    inland_water_mask,
) where {FT}
ClimaLand.InlandWater.SlabLakeModel(
    FT,
    domain,
    forcing,
    toml_dict::CP.ParamDict;
    inland_water_mask,
    prognostic_land_components,
    parameters,
)
```

## Inland Water Mask

```@docs
ClimaLand.InlandWater.inland_water_mask
```

## Boundary Conditions and Fluxes

```@docs
ClimaLand.InlandWater.AtmosDrivenLakeBC
ClimaLand.InlandWater.lake_boundary_fluxes!(
    bc::ClimaLand.InlandWater.AtmosDrivenLakeBC,
    prognostic_land_components::Val{(:lake,)},
    model::ClimaLand.InlandWater.SlabLakeModel{FT},
    Y,
    p,
    t,
) where {FT}
```

## Lake Parameterizations

```@docs
ClimaLand.InlandWater.lake_energy_at_freezing
ClimaLand.InlandWater.lake_liquid_fraction
ClimaLand.InlandWater.lake_energy_from_temperature
ClimaLand.InlandWater.lake_temperature
ClimaLand.InlandWater.lake_volumetric_internal_energy
ClimaLand.InlandWater.lake_runoff_energy_flux
ClimaLand.InlandWater.lake_surface_albedo
ClimaLand.InlandWater.lake_sediment_heat_flux
ClimaLand.InlandWater.lake_specific_humidity
```
