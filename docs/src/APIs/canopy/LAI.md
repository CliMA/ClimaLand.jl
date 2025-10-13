# LAI Models

```@meta
CurrentModule = ClimaLand.Canopy
```

## Models and Parameters

```@docs
ClimaLand.Canopy.AbstractLAIModel
ClimaLand.Canopy.OptimalLAIModel
ClimaLand.Canopy.OptimalLAIModel{FT}(parameters::OptimalLAIParameters{FT}) where {FT <: AbstractFloat}
ClimaLand.Canopy.OptimalLAIModel{FT}(toml_dict::CP.ParamDict) where {FT <: AbstractFloat}
ClimaLand.Canopy.OptimalLAIParameters
ClimaLand.Canopy.OptimalLAIParameters{FT}(toml_dict::CP.ParamDict) where {FT}
```

## Methods

```@docs
ClimaLand.Canopy.initialize_LAI!
ClimaLand.Canopy.update_LAI!
ClimaLand.Canopy.make_OptimalLAI_callback
ClimaLand.Canopy.call_update_optimal_LAI
```
