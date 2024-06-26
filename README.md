<picture>
  <source media="(prefers-color-scheme: dark)" srcset="logo-white.svg">
  <source media="(prefers-color-scheme: light)" srcset="logo.svg">
  <img alt="Shows the logo of ClimaLand, with a drop and three leaves" src="logo.svg">
</picture>
<p align="center">
  <strong>Create and run  land models in integrated (multi-
component) or standalone (single component) modes.  </strong>
</p>

[![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FClimaLand&query=total_requests&suffix=%2Ftotal&label=Downloads)](http://juliapkgstats.com/pkg/ClimaLand)

Recommended Julia Version: Stable release v1.10.0. CI no longer tests earlier
versions of Julia.

Certain features, including global runs, are not currently available on
Windows due to limitations with our regridding software.

Note that ClimaLand.jl is a different model from the original CliMA Land,
which aims to utilize remote sensing data through more complex canopy RT
and plant physiology modules. For more details, please refer to
https://github.com/CliMA/Land.
- Wang, Yujie, et al. "Testing stomatal models at the stand level in deciduous angiosperm and evergreen gymnosperm forests using CliMA Land (v0. 1)." Geoscientific Model Development 14.11 (2021): 6741-6763.
- R. K. Braghiere, Y. Wang, R. Doughty, D. Souza, T. Magney, J. Widlowski, M. Longo, A. Bloom, J. Worden, P. Gentine, and C. Frankenberg. 2021. Accounting for canopy structure improves hyperspectral radiative transfer and sun-induced chlorophyll fluorescence representations in a new generation Earth System model. Remote Sensing of Environment. 261: 112497.
- Wang, Yujie, and Christian Frankenberg. "On the impact of canopy model complexity on simulated carbon, water, and solar-induced chlorophyll fluorescence fluxes." Biogeosciences 19.1 (2022): 29-45.
- Wang, Yujie, et al. "GriddingMachine, a database and software for Earth system modeling at global and regional scales." Scientific data 9.1 (2022): 258.
- Holtzman, Nataniel, et al. "Constraining plant hydraulics with microwave radiometry in a land surface model: Impacts of temporal resolution." Water Resources Research 59.11 (2023): e2023WR035481.

## Models

Component Models:
RichardsModel <: AbstractSoilModel <: AbstractImExModel <: AbstractModel [runnable w/o LandModel wrapper as well]

EnergyHydrologyModel <: AbstractSoilModel <: AbstractImExModel <: AbstractModel [runnable w/o LandModel wrapper as well]

CanopyModel <: AbstractVegetationModel <: AbstractExpModel <: AbstractModel  [runnable w/o LandModel wrapper as well]

SnowModel <: AbstractSnowModel <: AbstractExpModel <: AbstractModel [runnable w/o LandModel wrapper as well]

Combined Models:

SoilCanopyModel <: AbstractLandModel <: AbstractImExModel <: AbstractModel (constructs the individual ComponentModels based on arguments)

|||
|---------------------:|:----------------------------------------------|
| **Docs Build**       | [![docs build][docs-bld-img]][docs-bld-url]   |
| **Documentation**    | [![dev][docs-dev-img]][docs-dev-url]          |
| **GHA CI**           | [![gha ci][gha-ci-img]][gha-ci-url]           |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]        |

[docs-bld-img]: https://github.com/CliMA/ClimaLand.jl/actions/workflows/docs.yml/badge.svg
[docs-bld-url]: https://github.com/CliMA/ClimaLand.jl/actions/workflows/docs.yml

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://CliMA.github.io/ClimaLand.jl/dev/

[gha-ci-img]: https://github.com/CliMA/ClimaLand.jl/actions/workflows/ci.yml/badge.svg
[gha-ci-url]: https://github.com/CliMA/ClimaLand.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/CliMA/ClimaLand.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/ClimaLand.jl
