# Autotrophic Respiration

```@meta
CurrentModule = ClimaLand.Canopy
```

## Models and Parameters

```@docs
ClimaLand.Canopy.AutotrophicRespirationModel
ClimaLand.Canopy.AutotrophicRespirationModel{FT}(toml_dict::CP.ParamDict) where {FT <: AbstractFloat}
ClimaLand.Canopy.AutotrophicRespirationParameters
ClimaLand.Canopy.AutotrophicRespirationParameters(toml_dict::CP.ParamDict; kwargs...)
```

## Methods

```@docs
ClimaLand.Canopy.compute_autrophic_respiration
ClimaLand.Canopy.plant_respiration_growth
```
