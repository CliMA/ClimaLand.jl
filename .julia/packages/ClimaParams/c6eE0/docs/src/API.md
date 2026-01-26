```@meta
CurrentModule = ClimaParams
```

# API

This page documents the functions and types available in ClimaParams.jl. The API is organized to follow a typical user workflow.

## 1. Core Data Structures

These are the main types for holding and interacting with parameters.

```@docs
ParamDict
Base.getindex
float_type
```

## 2. Creating a Parameter Dictionary

The primary entry point is [`create_toml_dict`](@ref), which can be customized by merging multiple files.

```@docs
create_toml_dict
merge_toml_files
```

## 3. Accessing Parameter Values

Once a [`ParamDict`](@ref) is created, you can retrieve parameter values in several ways. The most common method is [`get_parameter_values`](@ref).

```@docs
get_parameter_values
```

### Tag-Based Retrieval

Parameters can be organized in the TOML file with `tag` entries. These functions allow you to retrieve all parameters associated with one or more tags.

```@docs
get_tagged_parameter_values
get_tagged_parameter_names
fuzzy_match
```

## 4. Utilities for Integration and Reproducibility

These functions support logging, validation, and integration with user-defined parameter structs.

```@docs
log_parameter_information
write_log_file
check_override_parameter_usage
log_component!
create_parameter_struct
```