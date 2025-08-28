# Autotrophic Respiration

```@meta
CurrentModule = ClimaLand.Canopy
```

## Models and Parameters

```@docs
ClimaLand.Canopy.AutotrophicRespirationModel
ClimaLand.Canopy.AutotrophicRespirationModel{FT}() where {FT <: AbstractFloat}
ClimaLand.Canopy.AutotrophicRespirationParameters
ClimaLand.Canopy.AutotrophicRespirationParameters(
    ::Type{FT};
    kwargs...,
) where {FT <: AbstractFloat}
```

## Methods

```@docs
ClimaLand.Canopy.nitrogen_content
ClimaLand.Canopy.plant_respiration_maintenance
ClimaLand.Canopy.plant_respiration_growth
```
