# PlantHydraulics

```@meta
CurrentModule = ClimaLand.Canopy
```
## Models and Parameters

```@docs
ClimaLand.AbstractPlantHydraulicsModel
ClimaLand.PlantHydraulicsModel
ClimaLand.PlantHydraulicsParameters
```

## Plant Hydraulics Parameterizations

```@docs
ClimaLand.AbstractConductivityModel
ClimaLand.Weibull
ClimaLand.AbstractRetentionModel
ClimaLand.LinearRetentionCurve
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
ClimaLand.water_flux
ClimaLand.effective_saturation
ClimaLand.augmented_liquid_fraction
ClimaLand.water_retention_curve
ClimaLand.inverse_water_retention_curve
ClimaLand.root_water_flux_per_ground_area!
ClimaLand..hydraulic_conductivity
```
