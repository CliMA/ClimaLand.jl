# Soil Biogeochemistry

## This component model is available for use, but is still under development and is not yet fully debugged. Note that errors in this component do not propagate back to other component models.

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
ClimaLand.Soil.SoilCO2ModelParameters(toml_dict::CP.ParamDict)
```

## Model-specific Types

```@docs
ClimaLand.Soil.Biogeochemistry.MicrobeProduction
ClimaLand.Soil.Biogeochemistry.SoilCO2FluxBC
ClimaLand.Soil.Biogeochemistry.SoilCO2StateBC
ClimaLand.Soil.Biogeochemistry.AtmosCO2StateBC
ClimaLand.Soil.Biogeochemistry.AtmosO2StateBC
ClimaLand.Soil.Biogeochemistry.AbstractSoilDriver
ClimaLand.Soil.Biogeochemistry.SoilDrivers
ClimaLand.Soil.Biogeochemistry.PrescribedMet
```

## Functions of State

```@docs
ClimaLand.Soil.Biogeochemistry.volumetric_air_content
ClimaLand.Soil.Biogeochemistry.co2_diffusivity
ClimaLand.Soil.Biogeochemistry.microbe_source
ClimaLand.Soil.Biogeochemistry.o2_availability
ClimaLand.Soil.Biogeochemistry.o2_concentration
ClimaLand.Soil.Biogeochemistry.o2_fraction_from_concentration
```

## Extendible Functions

```@docs
ClimaLand.Soil.Biogeochemistry.soil_moisture
ClimaLand.Soil.Biogeochemistry.soil_temperature
```
