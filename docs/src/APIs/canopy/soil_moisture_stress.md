# Soil moisture stress

```@meta
CurrentModule = ClimaLand.Canopy
```

## Models and Parameters

```@docs
ClimaLand.Canopy.TuzetMoistureStressParameters
ClimaLand.Canopy.TuzetMoistureStressModel
ClimaLand.Canopy.TuzetMoistureStressModel{FT}()
ClimaLand.Canopy.TuzetMoistureStressModel{FT}(toml_dict; sc::FT = toml_dict["moisture_stress_sc"], pc::FT = toml_dict["moisture_stress_pc"])
ClimaLand.Canopy.PiecewiseMoistureStressParameters
ClimaLand.Canopy.PiecewiseMoistureStressModel
ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}()
Climaland.Canopy.PiecewiseMoistureStressModel{FT}(
    toml_dict::CP.AbstractTOMLDict;
    c::FT = toml_dict["moisture_stress_c"],
    porosity_residual = true,
    soil_params
) where {FT <: AbstractFloat}
ClimaLand.Canopy.NoMoistureStressModel
```

## Methods

```@docs
ClimaLand.Canopy.update_soil_moisture_stress!
ClimaLand.Canopy.update_piecewise_soil_moisture_stress!
```
