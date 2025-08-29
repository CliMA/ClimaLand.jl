# Soil Biogeochemistry

```@meta
CurrentModule = ClimaLand.Soil.Biogeochemistry
```
## Model Structure

```@docs
ClimaLand.Soil.Biogeochemistry.SoilCO2Model
```

## Parameter Structure

```@docs
ClimaLand.Soil.Biogeochemistry.SoilCO2ModelParameters
ClimaLand.Soil.SoilCO2ModelParameters(::Type{FT}) where {FT <: AbstractFloat}
```

## Model-specific Types

```@docs
ClimaLand.Soil.Biogeochemistry.MicrobeProduction
ClimaLand.Soil.Biogeochemistry.SoilCO2FluxBC
ClimaLand.Soil.Biogeochemistry.SoilCO2StateBC
ClimaLand.Soil.Biogeochemistry.AtmosCO2StateBC
ClimaLand.Soil.Biogeochemistry.AbstractSoilDriver
ClimaLand.Soil.Biogeochemistry.SoilDrivers
ClimaLand.Soil.Biogeochemistry.PrescribedMet
```

## Functions of State

```@docs
ClimaLand.Soil.Biogeochemistry.volumetric_air_content
ClimaLand.Soil.Biogeochemistry.co2_diffusivity
ClimaLand.Soil.Biogeochemistry.microbe_source
```

## Extendible Functions

```@docs
ClimaLand.Soil.Biogeochemistry.soil_moisture
ClimaLand.Soil.Biogeochemistry.soil_temperature
```
