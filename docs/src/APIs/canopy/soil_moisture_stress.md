# Soil moisture stress

```@meta
CurrentModule = ClimaLand.Canopy
```

## Models and Parameters

```@docs
ClimaLand.Canopy.TuzetMoistureStressModel
ClimaLand.Canopy.TuzetMoistureStressModel{FT}(toml_dict::CP.ParamDict)
ClimaLand.Canopy.PiecewiseMoistureStressModel
ClimaLand.Canopy.PiecewiseMoistureStressModel(domain, toml_dict::CP.ParamDict)
ClimaLand.Canopy.NoMoistureStressModel
```

## Methods

```@docs
ClimaLand.Canopy.update_soil_moisture_stress!
ClimaLand.Canopy.update_piecewise_soil_moisture_stress!
```
