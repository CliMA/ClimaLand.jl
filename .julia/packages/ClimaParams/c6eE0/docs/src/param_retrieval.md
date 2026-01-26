# Basic Parameter Retrieval

ClimaParams.jl provides a centralized system for managing climate model parameters. The core workflow involves creating a parameter dictionary and then retrieving specific parameters from it.

## Core Functions

There are three key functions for parameter retrieval:
- [`create_toml_dict`](@ref) constructs a parameter dictionary from TOML files
- [`get_parameter_values`](@ref) retrieves parameters from the dictionary

## Creating Parameter Dictionaries

To construct a parameter dictionary, pass in the desired floating point type. 
This will source parameter values from the global default list stored in `src/parameters.toml`

```@example howto
import ClimaParams as CP
toml_dict = CP.create_toml_dict(Float64)
nothing # hide
```

You can also specify custom override and default files:

```julia
# With custom files
toml_dict = CP.create_toml_dict(
    Float64,
    override_file = "my_parameters.toml",
    default_file = "default_parameters.toml"
)
```

## Retrieving Parameters

To retrieve parameters, pass in the TOML dictionary and the parameter names that match those in the TOML file.

```@example howto
params = CP.get_parameter_values(toml_dict, ["universal_gas_constant", "gravitational_acceleration"])
params.universal_gas_constant
params.gravitational_acceleration
nothing #hide
```

You can also use direct indexing to obtain values from the parameter dictionary:
```@example howto
toml_dict["gravitational_acceleration"]
```

## Name Maps

Name maps allow you to map global parameter names to local variable names for
convenience. This is especially useful when you want shorter, more intuitive
variable names in your code.

### Using NamedTuples

```@example howto
name_map = (;
    :gravitational_acceleration => :g,
    :angular_velocity_planet_rotation => :omega
)
params = CP.get_parameter_values(toml_dict, name_map)
params.g  # gives value field of gravitational_acceleration
params.omega
```

### Using Dictionaries

```@example howto
name_map = Dict("gravitational_acceleration" => "g", "angular_velocity_planet_rotation" => "omega")
params = CP.get_parameter_values(toml_dict, name_map)
nothing # hide
```

### Using Varargs

```@example howto
params = CP.get_parameter_values(toml_dict, 
    :gravitational_acceleration => :g,
    :angular_velocity_planet_rotation => :omega
)
nothing # hide
```

## Component Logging

You can specify a component name when retrieving parameters. This logs which parameters are used by which model component, which is useful for reproducibility:

```@example howto
params = CP.get_parameter_values(toml_dict, ["gravitational_acceleration"], "Ocean")
nothing # hide
```

## Tagged Parameters

ClimaParams supports parameter tagging for easy filtering. You can retrieve all parameters with a specific tag:

```@example howto
# Get all atmospheric parameters
atmospheric_params = CP.get_tagged_parameter_values(toml_dict, "atmosphere")

# Get parameters with multiple tags
physics_params = CP.get_tagged_parameter_values(toml_dict, ["atmosphere", "turbulence"])
nothing # hide
```

## Example Usage

### Simple Parameter Retrieval

Here's a basic example showing how to retrieve parameters for use in a simulation:

```julia
import ClimaParams as CP

# Create parameter dictionary
toml_dict = CP.create_toml_dict(Float64)

# Retrieve specific parameters
params = CP.get_parameter_values(toml_dict, [
    "gravitational_acceleration",
    "universal_gas_constant", 
    "planet_radius"
])

# Use parameters in your code
g = params.gravitational_acceleration
R = params.universal_gas_constant

# Alternatively, you can index directly into the parameter dict
toml_dict["gravitational_acceleration"]

```

## Parameter Structs

For more complex applications, you can build parameter structs that encapsulate
related parameters. Here's a complete example from the CliMA ecosystem:

### Building Parameter Structs

```julia
Base.@kwdef struct ThermodynamicsParameters{FT}
    universal_gas_constant::FT
    molmass_dryair::FT
    # derived parameters
    R_d::FT = universal_gas_constant / molmass_dryair
end

# Float-type constructor
ThermodynamicsParameters(::Type{FT}) = ThermodynamicsParameters(CP.create_toml_dict(FT))

# TOML dictionary constructor
function ThermodynamicsParameters(toml_dict::ParamDict{FT}) where {FT}
    return ThermodynamicsParameters{FT}(;
        temperature_triple_point = toml_dict["T_triple"],
        adiabatic_exponent_dry_air = toml_dict["kappa_d"],
        pressure_triple_point = toml_dict["press_triple"],
        thermodynamics_temperature_reference = toml_dict["T_0"],
        temperature_water_freeze = toml_dict["T_freeze"],
        isobaric_specific_heat_ice = toml_dict["cp_i"],
    )
end
nothing # hide
```

### Hierarchical Parameter Sets

For complex models with multiple components, you can build hierarchical
parameter sets that maintain parameter relationships:

```julia
# Build individual component parameter sets
thermodynamics_params = ThermodynamicsParameters(toml_dict)
params_0M = CloudMicrophysics.Microphysics_0M_Parameters(toml_dict)

# Combine into a hierarchical parameter set
parameter_set = CloudMicrophysics.CloudMicrophysicsParameters(
    toml_dict,
    params_0M,
    thermodynamics_params
    thermodynamics_params,
    microphysics_params
)
```

### Parameters-as-functions

Parameters can be accessed as functions for added flexibility:

```julia
K_therm(param_set) = param_set.K_therm
```

This can be useful for derived parameters:

```julia
derived_param(param_set) = param_set.param1 * param_set.param2
```

Or to forward parameters from nested parameter structs:

```julia
forwarded_param(ps::ParamSet) = ps.nested_params.forwarded_param
```

Functions can be autogenerated using `@eval`:

```julia
for fn in fieldnames(ParamSet)
    @eval $(fn)(ps::ParamSet) = ps.$(fn)
end
```

## Parameter Types

ClimaParams supports several parameter types:

- **float**: Numeric values (default)
- **integer**: Whole numbers
- **string**: Text values
- **bool**: Boolean values
- **datetime**: DateTime - an RFC 3339 formatted date-time with the offset omitted or an offset of `z`

The type is specified in the TOML file:

```toml
[gravitational_acceleration]
value = 9.81
type = "float"
description = "Gravitational acceleration on the planet (m s⁻²)."

[epoch_time]
value = 1970-01-01T00:00:00.0
type = "datetime"
description = "Unix epoch"
```
