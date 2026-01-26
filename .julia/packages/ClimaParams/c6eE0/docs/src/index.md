# ClimaParams.jl

A centralized parameter management system for climate modeling, ClimaParams.jl supports physical constants, planetary properties, and tunable parameters designed for calibration with data assimilation and machine learning tools.

## Overview

ClimaParams.jl provides a single source of truth for the parameters used in the [Climate Modeling Alliance (CliMA)](https://github.com/CliMA) ecosystem. By centralizing parameters across all model components (atmosphere, ocean, land, etc.), it enables joint calibration of interconnected climate processes through data assimilation and machine learning pipelines. This unified approach ensures that parameters shared between components remain consistent and can be optimized together, leading to more physically coherent model calibration.

The package manages two fundamental types of values:
- Physical and planetary *constants* (e.g., speed of light or planet radius)
- Tunable model *parameters* that can be calibrated individually or jointly across components

## What parameters are good candidates for ClimaParams?

A parameter is a good candidate for ClimaParams if it has _all_ of the following attributes:

 - The parameter does not vary in space
 - The parameter does not vary in time (per climate simulation)
 - The parameter is a function of only constants and other ClimaParams

## Getting Started

The basic flow is as follows:
1. Create the parameter dictionary with your desired floating point type
2. Retrieve parameters

```@example howto
import ClimaParams as CP

# Create parameter dictionary with default values
param_dict = CP.create_toml_dict(Float64)

# Retrieve physical constants
constants = CP.get_parameter_values(
    param_dict, 
    ["gravitational_acceleration", "planet_radius", "light_speed"]
)

# Retrieve parameters with custom names
params = CP.get_parameter_values(
    param_dict,
    Dict("universal_gas_constant" => "R", "gravitational_acceleration" => "g")
)
```
Or retrieve all parameters associated with a specific tag, for example, "SurfaceFluxes"

```@example howto
sf_params = CP.get_tagged_parameter_values(param_dict, "surfacefluxes") 
```

## Best Practices

1. **Use descriptive parameter names** in your TOML files
2. **Log component usage** by specifying component names
3. **Use name maps** for shorter, more intuitive variable names
4. **Tag related parameters** for easy filtering
5. **Create parameter structs** to encapsulate related parameters and their relationships
6. **Validate parameter values** before using them in simulations

## Documentation Overview

### Core Usage

- **[Parameter retrieval](param_retrieval.md)**: Learn how to retrieve parameters from dictionaries, use name maps, and work with tagged parameters. Includes practical examples from the CliMA ecosystem.

### Configuration and File Format

- **[TOML file interface](toml.md)**: Understand how to define parameters in TOML files, including value types, descriptions, tags, and advanced features for calibration workflows.

### API Reference

- **[API](API.md)**: Complete reference documentation for all functions and types in ClimaParams.jl, organized by functionality and typical usage patterns.

### Advanced Topics

- **[Parameter Structs](param_retrieval.md#Parameter-Structs)**: Learn how to create custom parameter structs for your models.
- **[Component Logging](param_retrieval.md#Component-Logging)**: Understand how to track parameter usage across model components for logging and validation.
- **[Override Files](toml.md#Override-Files)**: See how to customize parameters for specific experiments.
- **[Tagged Parameters](param_retrieval.md#Tagged-Parameters)**: Discover how to organize and retrieve related parameters.

For detailed usage examples and integration into your code, start with the Parameter retrieval guide, then explore the TOML file interface for configuration details.
