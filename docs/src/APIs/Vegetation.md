# Vegetation

```@meta
CurrentModule = ClimaLSM.Roots
```
## Models

```@docs
ClimaLSM.Roots.AbstractVegetationModel
ClimaLSM.Roots.RootsModel
```

## Roots Diagnostic Variables

```@docs
flux,
ClimaLSM.Roots.effective_saturation
ClimaLSM.Roots.augmented_liquid_fraction
ClimaLSM.Roots.van_genuchten_volume_to_pressure
ClimaLSM.Roots.van_genuchten_pressure_to_volume
ClimaLSM.Roots.ϑ_l_to_absolute_pressure
ClimaLSM.Roots.absolute_pressure_to_ϑ_l
ClimaLSM.Roots.flux_out_roots
```

## Roots Parameters

```@docs
ClimaLSM.Roots.RootsParameters
```

## Roots Methods and Types

```@docs
ClimaLSM.Roots.flux_out_roots
ClimaLSM.Roots.PrescribedSoilPressure
ClimaLSM.Roots.PrescribedTranspiration
ClimaLSM.Roots.AbstractRootExtraction
```