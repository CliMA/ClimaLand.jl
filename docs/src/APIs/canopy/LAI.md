# Optimal LAI

```@meta
CurrentModule = ClimaLand.Canopy
```

## Models and Parameters

```@docs
ClimaLand.Canopy.ZhouOptimalLAIModel
ClimaLand.Canopy.OptimalLAIParameters
ClimaLand.Canopy.OptimalLAIParameters{FT}(toml_dict::CP.ParamDict) where {FT}
```

## Methods

```@docs
ClimaLand.Canopy.update_optimal_LAI
ClimaLand.Canopy.compute_A0_daily
ClimaLand.Canopy.compute_LAI
ClimaLand.Canopy.make_OptimalLAI_callback
ClimaLand.Canopy.call_update_optimal_LAI
ClimaLand.Canopy.compute_L_max
ClimaLand.Canopy.compute_m
ClimaLand.Canopy.compute_steady_state_LAI
ClimaLand.Canopy.lambertw0
```
