# Available Models and Parameterizations

## Standalone and integrated models
As described in the [ClimaLand fundamental concepts](@ref) page, ClimaLand contains a variety of standalone
and integrated models. A complete list of available models is provided here.

### ClimaLand standalone models

This table shows the abstract type (where applicable) and concrete type for each model.
Abstract types are a concept originating in object-oriented programming
that allow us to create a hierarchy of model types.
For example, we can define an abstract type `AbstractSoilModel`, and then
create specific models `RichardsModel` and `EnergyHydrology` as subtypes.
This helps organize models and makes it easier to write functions that work
with any model in the group via multiple dispatch on the abstract type.
For more information about abstract types, please see the [Julia manual](https://docs.julialang.org/en/v1/manual/types/#man-abstract-types).

| Abstract model type                | Model name                |
|------------------------------------|---------------------------|
| `AbstractSoilModel`                | `RichardsModel`           |
|                                    | `EnergyHydrology`         |
| N/A                                | `CanopyModel`             |
| `AbstractSnowModel`                | `SnowModel`               |
| `AbstractSoilBiogeochemistryModel` | `SoilCO2Model`            |
| `AbstractBucketModel`              | `BucketModel`             |
| `AbstractSurfaceWaterModel`        | `PondModel`               |

### ClimaLand integrated models
| Integrated model name      | Component standalone models |
|----------------------------|-----------------------------|
| `LandModel`                | `EnergyHydrology`           |
|                            | `CanopyModel`               |
|                            | `SoilCO2Model`              |
|                            | `SnowModel`                 |
| `SoilCanopyModel`          | `EnergyHydrology`           |
|                            | `CanopyModel`               |
|                            | `SoilCO2Model`              |
| `SoilSnowModel`            | `EnergyHydrology`           |
|                            | `SnowModel`                 |
| `LandSoilBiogeochemistry`  | `EnergyHydrologyModel`      |
|                            | `SoilCO2Model`              |
| `LandHydrology`            | `RichardsModel`             |
|                            | `PondModel`                 |

## Model parameterizations

Each model contains one or more parameterizations, often with several options
to choose from. The following tables show all available parameterizations
to clarify the different model setups.

We also note the default option for each parameterization.
This is what the model will use automatically when built with just the basic
inputs (domain, forcing data, and TOML parameter file), without customizing
the parameterizations yourself.
The defaults are only explicitly shown here for standalone models but
are also inherited by any integrated model containing a given standalone model,
unless otherwise noted.

### Soil Models

#### `EnergyHydrology`
| Parameterization type                | Available options                 |
|--------------------------------------|-----------------------------------|
| `AbstractSoilAlbedoParameterization` | `CLMTwoBandSoilAlbedo` (default)  |
|                                      | `ConstantTwoBandSoilAlbedo`       |
| `AbstractRunoffModel`                | `TOPMODELRunoff` (default)        |
|                                      | `SurfaceRunoff`                   |
|                                      | `NoRunoff`                        |
| `AbstractSoilHydrologyClosure`       | `vanGenuchten` (default)          |
|                                      | `BrooksCorey`                     |

#### `RichardsModel`
| Parameterization type                | Available options                 |
|--------------------------------------|-----------------------------------|
| `AbstractRunoffModel`                | `TOPMODELRunoff` (default)        |
|                                      | `SurfaceRunoff`                   |
|                                      | `NoRunoff`                        |
| `AbstractSoilHydrologyClosure`       | `vanGenuchten` (default)          |
|                                      | `BrooksCorey`                     |

### `CanopyModel`
| Parameterization type                 | Available options                                           |
| ------------------------------------- | ------------------------------------------------------------|
| `AbstractAutotrophicRespirationModel` | `AutotrophicRespirationModel` (default)                     |
| `AbstractRadiationModel`              | `TwoStreamModel` (default)                                  |
|                                       | `BeerLambertModel`                                          |
| `AbstractPhotosynthesisModel`         | `FarquharModel` (default)                                   |
|                                       | `PModel`                                                    |
| `AbstractStomatalConductanceModel`    | `MedlynConductanceModel` (default)                          |
|                                       | `PModelConductance`                                         |
| `AbstractPlantHydraulicsModel`        | `PlantHydraulicsModel` (default)                            |
| `AbstractSoilMoistureStressModel`     | `TuzetMoistureStressModel` (default)                        |
|                                       | `PiecewiseMoistureStressModel`                              |
| `AbstractSIFModel`                    | `Lee2015SIFModel` (default)                                 |
| `AbstractCanopyEnergyModel`           | `PrescribedCanopyTempModel` (default, in standalone canopy) |
|                                       | `BigLeafEnergyModel` (default, in integrated models)        |

### `SnowModel`
| Parameterization type            | Available options                       |
| -------------------------------- | ----------------------------------------|
| `AbstractDensityModel`           | `MinimumDensityModel` (default)         |
|                                  | `NeuralDepthModel`                      |
| `AbstractAlbedoModel`            | `ConstantAlbedoModel` (default)         |
|                                  | `ZenithAngleAlbedoModel`                |
| `AbstractSnowCoverFractionModel` | `WuWuSnowCoverFractionModel` (default)  |
