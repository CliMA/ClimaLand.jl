<div align="center">
  <img src="docs/src/assets/logo.svg" alt="ClimaParams.jl Logo" width="128" height="128">
</div>

# ClimaParams.jl

A centralized parameter management system for climate modeling. ClimaParams.jl supports physical constants, planetary properties, and tunable parameters designed for calibration with data assimilation and machine learning tools. 

|                           |                                                                          |
|--------------------------:|:-------------------------------------------------------------------------|
| **Stable Release**        | [![stable][stable-img]][stable-url] [![docs-stable][docs-stable-img]][docs-stable-url] |
| **Latest Documentation**  | [![dev][docs-latest-img]][docs-latest-url]                                |
| **Unit Tests**            | [![unit tests][unit-tests-img]][unit-tests-url] [![codecov][codecov-img]][codecov-url] |
| **Downloads**             | [![Downloads][dlt-img]][dlt-url]                                          |

[docs-latest-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-latest-url]: https://CliMA.github.io/ClimaParams.jl/dev/

[stable-img]: https://img.shields.io/github/v/release/CliMA/ClimaParams.jl?label=stable
[stable-url]: https://github.com/CliMA/ClimaParams.jl/releases/latest

[docs-stable-img]: https://img.shields.io/badge/docs-stable-green.svg
[docs-stable-url]: https://CliMA.github.io/ClimaParams.jl/stable/

[unit-tests-img]: https://github.com/CliMA/ClimaParams.jl/actions/workflows/ci.yml/badge.svg
[unit-tests-url]: https://github.com/CliMA/ClimaParams.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/CliMA/ClimaParams.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/ClimaParams.jl

[dlt-img]: https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FClimaParams&query=total_requests&label=Downloads
[dlt-url]: https://juliapkgstats.com/pkg/ClimaParams

## Overview

ClimaParams.jl provides a single source of truth for the parameters used in the [Climate Modeling Alliance (CliMA)](https://github.com/CliMA) ecosystem. By centralizing parameters across all model components (atmosphere, ocean, land, etc.), it enables joint calibration of interconnected climate processes through data assimilation and machine learning pipelines. This unified approach ensures that parameters shared between components remain consistent and can be optimized together, leading to more physically coherent model calibration.

The package manages two fundamental types of values:
- Physical and planetary *constants* (e.g., speed of light or planet radius)
- Tunable model *parameters* that can be calibrated individually or jointly across components

## Features

- **Centralized Management**: A single, authoritative source for all model constants and parameters, enabling joint calibration of interconnected processes.
- **Type-Safe Retrieval**: Guaranteed type-safety with automatic validation and conversion for floating-point, integer, string, and boolean types.
- **Override System**: Support for parameter overrides with precedence handling for flexible experimentation.
- **Parameter Tagging**: Logically group values by model component (e.g., atmosphere, land) for easy filtering and retrieval.
- **Machine Learning Integration**: Seamless integration with data assimilation and ML calibration workflows, including joint parameter optimization across components.
- **TOML Configuration**: Human-readable TOML files define parameters and their metadata.
- **Reproducibility**: Automatic logging of parameter sets used in model runs to ensure scientific reproducibility.
- **Multi-Planet Support**: An extensible framework for simulations of Earth and other planetary bodies.

## Quick Example

```julia
using ClimaParams

# Floating-point type for parameters
FT = Float64

# Create a dictionary containing all default parameters, cast to the chosen float type
param_dict = create_toml_dict(FT)

# Retrieve a struct of physical constants by name
constants = get_parameter_values(
    param_dict, 
    ["gravitational_acceleration", "planet_radius", "light_speed"],
)
# Access constants, e.g.: 
constants.gravitational_acceleration

# Retrieve parameters and assign them custom names for convenience
params = get_parameter_values(
    param_dict,
    Dict("universal_gas_constant" => "R", "gravitational_acceleration" => "g"),
)
# Access parameters, e.g. 
params.R

# Get all parameters associated with a specific tag
atmospheric_params = get_tagged_parameter_values(param_dict, "atmosphere")
```

## Parameter Categories

ClimaParams manages two main categories of parameters:

- **Constants**: Immutable physical values. This includes universal constants (e.g., speed of light) and planet-specific properties (e.g., gravitational acceleration or planetary radius)
- **Model Parameters**: Tunable values that are subject to calibration via data assimilation or machine learning. These often correspond to specific physical processes, including:
  - Atmospheric physics 
  - Land surface 
  - Turbulence closures
  - Biogeochemistry

## Integration with CliMA Models

By providing a centralized source of parameters, ClimaParams.jl is essential for the interoperability and reliability of the CliMA ecosystem. It is designed to:
- Enable joint calibration of parameters across model components through unified data assimilation and machine learning pipelines
- Ensure consistent parameter usage across all models (atmosphere, ocean, land, etc.)
- Facilitate parameter sensitivity analysis and uncertainty quantification studies, jointly across model components
- Enable reproducible experiments by explicitly logging the exact parameter sets used for each simulation