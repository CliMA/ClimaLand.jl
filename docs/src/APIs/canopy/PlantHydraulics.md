# PlantHydraulics

```@meta
CurrentModule = ClimaLand.Canopy
```
## Models and Parameters

```@docs
ClimaLand.Canopy.AbstractPlantHydraulicsModel
ClimaLand.Canopy.PlantHydraulicsModel
ClimaLand.Canopy.PlantHydraulicsParameters
```

## Plant Hydraulics Parameterizations

```@docs
ClimaLand.Canopy.AbstractConductivityModel
ClimaLand.Canopy.Weibull
ClimaLand.Canopy.AbstractRetentionModel
ClimaLand.Canopy.LinearRetentionCurve
```

## Constructor Methods

```@docs
ClimaLand.Canopy.PlantHydraulicsModel{FT}(
    domain,
    toml_dict::CP.ParamDict;
) where {FT <: AbstractFloat}
```

## Plant Hydraulics Diagnostic Variables

```@docs
ClimaLand.Canopy.water_flux
ClimaLand.Canopy.effective_saturation
ClimaLand.Canopy.augmented_liquid_fraction
ClimaLand.Canopy.water_retention_curve
ClimaLand.Canopy.inverse_water_retention_curve
ClimaLand.Canopy.root_water_flux_per_ground_area!
ClimaLand.Canopy.hydraulic_conductivity
```
