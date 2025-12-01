# PlantHydraulics

```@meta
CurrentModule = ClimaLand.PlantHydraulics
```
## Models and Parameters

```@docs
ClimaLand.PlantHydraulics.AbstractPlantHydraulicsModel
ClimaLand.PlantHydraulics.SteadyStateModel
ClimaLand.PlantHydraulics.PlantHydraulicsModel
ClimaLand.PlantHydraulics.PlantHydraulicsParameters
```

## Prognostic Plant Hydraulics Parameterizations

```@docs
ClimaLand.PlantHydraulics.AbstractConductivityModel
ClimaLand.PlantHydraulics.Weibull
ClimaLand.PlantHydraulics.AbstractRetentionModel
ClimaLand.PlantHydraulics.LinearRetentionCurve
ClimaLand.PlantHydraulics.AbstractTranspiration
ClimaLand.PlantHydraulics.DiagnosticTranspiration
ClimaLand.PlantHydraulics.PrescribedTranspiration
```

## Constructor Methods

```@docs
ClimaLand.Canopy.PlantHydraulicsModel{FT}(
    domain,
    toml_dict::CP.ParamDict;
) where {FT <: AbstractFloat}
ClimaLand.Canopy.SteadyStateModel{FT}()
```

## Prognostic Plant Hydraulics Diagnostic Variables

```@docs
ClimaLand.PlantHydraulics.water_flux
ClimaLand.PlantHydraulics.effective_saturation
ClimaLand.PlantHydraulics.augmented_liquid_fraction
ClimaLand.PlantHydraulics.water_retention_curve
ClimaLand.PlantHydraulics.inverse_water_retention_curve
ClimaLand.PlantHydraulics.root_water_flux_per_ground_area!
ClimaLand.PlantHydraulics.hydraulic_conductivity
```
