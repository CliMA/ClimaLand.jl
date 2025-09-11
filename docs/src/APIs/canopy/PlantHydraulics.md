# PlantHydraulics

```@meta
CurrentModule = ClimaLand.PlantHydraulics
```
## Models and Parameters

```@docs
ClimaLand.PlantHydraulics.PlantHydraulicsModel
ClimaLand.PlantHydraulics.PlantHydraulicsParameters
```

## Plant Hydraulics Parameterizations

```@docs
ClimaLand.PlantHydraulics.Weibull
ClimaLand.PlantHydraulics.LinearRetentionCurve
ClimaLand.PlantHydraulics.PrescribedSiteAreaIndex
ClimaLand.PlantHydraulics.AbstractTranspiration
ClimaLand.PlantHydraulics.DiagnosticTranspiration
```

## Constructor Methods

```@docs
ClimaLand.Canopy.PlantHydraulicsModel{FT}(
    domain,
    LAI::AbstractTimeVaryingInput,
    toml_dict::CP.ParamDict;
) where {FT <: AbstractFloat}
```

## Plant Hydraulics Diagnostic Variables

```@docs
ClimaLand.PlantHydraulics.water_flux
ClimaLand.PlantHydraulics.effective_saturation
ClimaLand.PlantHydraulics.augmented_liquid_fraction
ClimaLand.PlantHydraulics.water_retention_curve
ClimaLand.PlantHydraulics.inverse_water_retention_curve
ClimaLand.PlantHydraulics.root_water_flux_per_ground_area!
ClimaLand.PlantHydraulics.hydraulic_conductivity
ClimaLand.PlantHydraulics.root_distribution
```
