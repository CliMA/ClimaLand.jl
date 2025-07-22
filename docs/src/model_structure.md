# ClimaLand model structure

Here, we explain the internal organization of ClimaLand models.

## Standalone and integrated models
A core concept of ClimaLand's design is the presence of _standalone_ and _integrated_ models.
In ClimaLand's framework, the land system is divided into individual components: soil, vegetation (or canopy), and snow. Each of these components can be run individually (_standalone_) or combined together (_integrated_).

As an example, the soil model `EnergyHydrology` and the canopy model `CanopyModel` can each be set up and run on their own in standalone mode. Or, a user can set up the integrated `SoilCanopyModel`, which internally constructs each of these models and runs them in tandem, computing relevant interaction terms at each timestep.

## Model variables

All ClimaLand models (standalone and integrated alike) contain two types of variables:

### Prognostic variables
These variables are solved for at each timestep based on differential equations and using timestepping methods implemented in [ClimaTimeSteppers.jl](https://github.com/CliMA/ClimaTimeSteppers.jl). The specific timestepping methods used depend on the model being used, and can be chosen from the many that ClimaTimeSteppers.jl has implemented.

Prognostic variables are stored in the model state, which is represented as the variable `Y`.

For more information about timestepping, please see the [timestepping tutorial](https://clima.github.io/ClimaLand.jl/stable/generated/shared_utilities/timestepping/).

### Diagnostic variables
These variables are computed at each timestep as functions of the prognostic variables (as well as other parameters, forcing data, etc).

Diagnostic variables are stored in the model cache, which is represented as the variable `p`.

## Model internals
All standalone models in the ClimaLand ecosystem contain the following objects:
- domain: The physical space on which the simulation is being evaluated; ClimaLand domains use [ClimaCore.jl](https://github.com/CliMA/ClimaCore.jl) spaces internally.
- parameters: A set of constants that may or may not vary in space, and are required to solve the model equations.
- boundary conditions: A ClimaLand object containing information about how to handle computation at the edges of the domain.
Boundary conditions arise as a result of solving PDEs for 3D models, and we maintain this terminology for single layer models.
In general, boundary conditions specify how the component interacts with surrounding pieces of the model, whether they are
other component models or the atmosphere.

Standalone models may contain additional information, depending on the particular type of model.
For example, soil models contain a flag `lateral_flow` to control lateral flow in the soil, but this flag is not included in canopy models (where it wouldn't make sense).

Integrated models don't have this same structure, but instead contain each of the component models they consist of.
These component models, though run in integrated mode, still each contain the structure described above.
The result is a nested structure for integrated ClimaLand models.
