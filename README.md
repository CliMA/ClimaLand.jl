# ClimaLSM
an in-progress prototype interface for running land models in integrated (multi-
component) or standalone (single component) modes.

## Models
```
Model <: AbstractModel
```

All models have the same abstract supertype, and have shared functionality.
There is the the option to define new methods particular to your model, or to
fall back on defaults. The functions each new Model type can define are:
- make_rhs()
- make_update_aux()
- initialize_prognostic()
- initialize_auxiliary()
- initialize()
- prognostic_vars()
- auxiliary_vars()
    
Each model will also have some notion of a domain with coordinates, parameter sets,
and boundary conditions or other prescribed drivers.

Examples:

Component Models:
RichardsModel <: AbstractSoilModel <: AbstractModel [runnable w/o LandModel wrapper as well]

PlantHydraulicsModel <: AbstractVegetationModel <: AbstractModel  [runnable w/o LandModel wrapper as well]

PondModel <: AbstractSurfaceWaterModel  <: AbstractModel  [runnable w/o LandModel wrapper as well]

Combined Models:

SoilPlantHydrologyModel <: AbstractModel (constructs the individual ComponentModels based on arguments)

|||
|---------------------:|:----------------------------------------------|
| **Docs Build**       | [![docs build][docs-bld-img]][docs-bld-url]   |
| **Documentation**    | [![dev][docs-dev-img]][docs-dev-url]          |
| **GHA CI**           | [![gha ci][gha-ci-img]][gha-ci-url]           |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]        |
| **Bors enabled**     | [![bors][bors-img]][bors-url]                 |

[docs-bld-img]: https://github.com/CliMA/ClimaLSM.jl/actions/workflows/docs.yml/badge.svg
[docs-bld-url]: https://github.com/CliMA/ClimaLSM.jl/actions/workflows/docs.yml

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://CliMA.github.io/ClimaLSM.jl/dev/

[gha-ci-img]: https://github.com/CliMA/ClimaLSM.jl/actions/workflows/ci.yml/badge.svg
[gha-ci-url]: https://github.com/CliMA/ClimaLSM.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/CliMA/ClimaLSM.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/ClimaLSM.jl

[bors-img]: https://bors.tech/images/badge_small.svg
[bors-url]: https://app.bors.tech/repositories/40649
