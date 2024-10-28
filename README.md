<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="logo-white.svg">
    <source media="(prefers-color-scheme: light)" srcset="logo.svg">
    <img alt="Shows the logo of ClimaLand, with a water drop and three leaves" src="logo.svg" width="500">
  </picture>
</h1>
<p align="center">
  <strong>Create and run  land models in integrated (multi-
component) or standalone (single component) modes.  </strong>
</p>

<div align="center">

|||
|---------------------:|:----------------------------------------------|
| **Documentation**    | [![dev][docs-stable-img]][docs-stable-url]          |
| **GHA CI**           | [![gha ci][gha-ci-img]][gha-ci-url]           |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]        |
| **Downloads**        | [![downloads][downloads-img]][downloads-url]  |

</div>

# Introduction

This is the repository of the CliMA land model code. Here are some notable features:
- ClimaLand has a modular design, models can be run as standalone (e.g., soil moisture only) or integrated (e.g., soil moisture and energy AND canopy AND snow, etc.)
- ClimaLand can simulate single columns, regional boxes, and global runs
- ClimaLand is CPU and GPU compatible
- ClimaLand welcome contributions: please feel free to reach out to us with questions about how to get started, create a branch, and extend our code. For example, a modeler might want to test a new stomatal conductance model.
- ClimaLand provides APIs and UIs at multiple levels.

## Installation

To use ClimaLand.jl, first you need to [install Julia](https://julialang.org/downloads/).
Then, you can install ClimaLand.jl by doing:

```Julia
julia> using Pkg
julia> Pkg.add(ClimaLand)
```

Which is equivalent to doing

```Julia
julia> ] # typing the ] key with go into package REPL mode
pkg> add ClimaLand
```

You are now ready to use `ClimaLand.jl`. To get started, we recommend reading the [documentation](https://clima.github.io/ClimaLand.jl/dev/).

## Models

In our code base, a "model" define a set of prognostic variables which must be timestepped. The equations which govern the time evolution likely contain parameters and are informed by parameterization and physical domain choices. Any ClimaLand model contains all of the information needed to evaluate these equations. Below are the current models we support:

<strong> Component Models: </strong>

- RichardsModel <: AbstractSoilModel <: AbstractModel (runnable only in standalone mode)

- EnergyHydrologyModel <: AbstractSoilModel <: AbstractModel (runnable in standalone mode, or as part of a land model)

- CanopyModel <: AbstractVegetationModel <: AbstractModel  (runnable in standalone mode, or as part of a land model)

- SnowModel <: AbstractSnowModel <: AbstractModel (runnable in standalone mode, or as part of a land model)

<strong> Combined Models: </strong>

- SoilCanopyModel <: AbstractLandModel <: AbstractModel (an example of a land model, made of individual component models which are solved simultaneously but taking into account interactions between the components)

## Notes

Recommended Julia Version: Stable release v1.11.1. CI tests Julia v1.10 and 1.11.

ClimaLand.jl is a different model from the original CliMA Land,
which aims to utilize remote sensing data through more complex canopy RT
and plant physiology modules. For more details, please refer to
https://github.com/CliMA/Land.
- Wang, Yujie, et al. "Testing stomatal models at the stand level in deciduous angiosperm and evergreen gymnosperm forests using CliMA Land (v0. 1)." Geoscientific Model Development 14.11 (2021): 6741-6763.
- R. K. Braghiere, Y. Wang, R. Doughty, D. Souza, T. Magney, J. Widlowski, M. Longo, A. Bloom, J. Worden, P. Gentine, and C. Frankenberg. 2021. Accounting for canopy structure improves hyperspectral radiative transfer and sun-induced chlorophyll fluorescence representations in a new generation Earth System model. Remote Sensing of Environment. 261: 112497.
- Wang, Yujie, and Christian Frankenberg. "On the impact of canopy model complexity on simulated carbon, water, and solar-induced chlorophyll fluorescence fluxes." Biogeosciences 19.1 (2022): 29-45.
- Wang, Yujie, et al. "GriddingMachine, a database and software for Earth system modeling at global and regional scales." Scientific data 9.1 (2022): 258.
- Holtzman, Nataniel, et al. "Constraining plant hydraulics with microwave radiometry in a land surface model: Impacts of temporal resolution." Water Resources Research 59.11 (2023): e2023WR035481.
- Braghiere, R. K., Wang, Y., Gagné-Landmann, A., Brodrick, P. G., Bloom, A. A., Norton, A. J., Ma, S., Levine, P., Longo, M., Deck, K., Gentine, P., Worden, J. R., Frankenberg, C., & Schneider, T. (2023). The Importance of Hyperspectral Soil Albedo Information for Improving Earth System Model Projections. AGU Advances, 4(4), e2023AV000910. [link](https://doi.org/10.1029/2023AV000910)
- Wang, Y., Braghiere, R. K., Yin, Y., Yao, Y., Hao, D., & Frankenberg, C. (2024). Beyond the visible: Accounting for ultraviolet and far-red radiation in vegetation productivity and surface energy budgets. Global Change Biology, 30(5), e17346. [link](https://doi.org/10.1111/GCB.17346)

[docs-bld-img]: https://github.com/CliMA/ClimaLand.jl/actions/workflows/docs.yml/badge.svg
[docs-bld-url]: https://github.com/CliMA/ClimaLand.jl/actions/workflows/docs.yml

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://CliMA.github.io/ClimaLand.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://CliMA.github.io/ClimaLand.jl/stable/

[gha-ci-img]: https://github.com/CliMA/ClimaLand.jl/actions/workflows/ci.yml/badge.svg
[gha-ci-url]: https://github.com/CliMA/ClimaLand.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/CliMA/ClimaLand.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/ClimaLand.jl

[downloads-img]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FClimaLand&query=total_requests&suffix=%2Ftotal&label=Downloads
[downloads-url]: http://juliapkgstats.com/pkg/ClimaLand
