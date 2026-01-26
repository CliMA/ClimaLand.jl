# The TOML Parameter File Interface

Parameters for CliMA models are defined in `.toml` files. ClimaParams.jl is designed to work with two main sources of parameters, which are merged together:

1.  **A default parameter file**: This file is bundled with ClimaParams.jl and contains a comprehensive list of default values for the entire CliMA ecosystem.
2.  **A user-defined override file**: This file is provided by the user for a specific experiment. It only needs to contain the parameters that deviate from the defaults.

## Parameter Format

Each parameter is defined by its unique name as a TOML table header (e.g., `[my_parameter_name]`). It can have the following attributes:

- `value`: (Required) The value of the parameter. Can be a scalar or an array.
- `type`: (Required) The data type. See [supported types](param_retrieval.md#Parameter-Types).
- `description`: (Recommended) A string explaining the parameter's purpose and its physical units.
- `tag`: An optional array of strings used to group related parameters.

Additional attributes, for example, used by [EnsembleKalmanProcesses.jl](https://github.com/CliMA/EnsembleKalmanProcesses.jl), may include:
- `prior`: An optional string describing a prior distribution, for use in calibration and data assimilation workflows.
- `transformation`: An optional string describing a transformation for the parameter, used in calibration.

!!! warn "On Array Types"
    Array values use the same `type` declaration as their scalar counterparts. For example, a vector of floats is specified with `type = "float"`.

### Basic Parameter Definition

At a minimum, a parameter requires a `value` and a `type`.

```TOML
[molar_mass_dry_air]
value = 0.03
type = "float"
```

It is highly recommended to include a `description` with units (CliMA generally uses SI units), as found in the default parameter files.

```TOML
[molar_mass_dry_air]
value = 0.02897
type = "float"
description = "Molecular weight of dry air (kg/mol)"
```

### Tagging Parameters

Tags provide a way to group related parameters. They do not create namespaces, and all parameter names must remain globally unique. To add tags, provide a list of strings to the `tag` field.

A recommended convention is to tag parameters with the model component(s) where they are used.

```TOML
[prandtl_number_0_grachev]
value = 0.98
type = "float"
description = "The turbulent Prandtl number in neutral conditions ($Pr_0$) for the Grachev universal functions (unitless). Source: Grachev et al. (2007), DOI: 10.1007/s10546-007-9177-6."
tag = ["SurfaceFluxes"]
```

Parameters with a specific tag can then be retrieved easily in Julia. Tag matching is case-insensitive and ignores punctuation. For more information, see the API for [`fuzzy_match`](@ref).

```julia
# This will retrieve all parameters tagged with "SurfaceFluxes"
sf_params = get_tagged_parameter_values(toml_dict, "surfacefluxes")
```

### Advanced: ClimaParams.jl for Calibration

For calibration workflows, parameters can include additional metadata to guide the calibration process:

```TOML
[entrainment_parameter]
value = 0.2
type = "float"
description = "Entrainment rate parameter for convective plumes"
tag = ["Convection", "Turbulence"]
prior = "LogNormal(-1.6, 0.4)"
transformation = "log"
```

The `prior` and `transformation` fields help guide the calibration process in [EnsembleKalmanProcesses.jl](https://github.com/CliMA/EnsembleKalmanProcesses.jl).

## Override Files

When an override file is provided, its values for any given parameter **take precedence** over the default values. Other attributes from the default file (like `description` or `tag`) are merged if they are not present in the override file.

### Override Mechanism

For example, if the user's override file contains:
```TOML
[molar_mass_dry_air]
value = 0.03
type = "float"
```
The final, merged parameter used in the simulation will be:
```TOML
[molar_mass_dry_air]
value = 0.03  # <-- Overwritten by the user's value
type = "float"
description = "Molar mass of dry air (kg mol⁻¹)." # <-- Merged from the default file
```

## Interacting with Parameters in Julia

ClimaParams.jl provides a clear workflow for using parameters in your code.

### 1. Loading Parameters

The main entry point is [`create_toml_dict`](@ref), which loads, merges, and types the parameters.

```julia
create_toml_dict(FT; override_file=nothing, default_file=...)
```

The first argument, `FT`, must be a float type (e.g., `Float64` or `Float32`) and determines the precision of all floating-point parameters.

A typical use case involves providing the path to a local override file:
```julia
import ClimaParams

FT = Float64
local_experiment_file = joinpath(@__DIR__, "local_exp_parameters.toml")
toml_dict = ClimaParams.create_toml_dict(FT; override_file = local_experiment_file)
```

If `override_file` is omitted, only the default parameters are loaded. You can also pass Julia `Dict`s directly instead of file paths. To combine more than two files, see the API for `merge_toml_files`.

### 2. Using and Logging Parameters

The returned `toml_dict` is then used to construct parameter structs for different model components.

```julia
# Retrieve values and construct the component-specific parameter struct
thermo_params = Thermodynamics.ThermodynamicsParameters(toml_dict)

# ... build the rest of the model components ...

# After all components are built, log the used parameters before running
log_file = joinpath(@__DIR__, "parameter_log.toml")
ClimaParams.log_parameter_information(toml_dict, log_file)

# ... run_model(...) ...
```

The function [`log_parameter_information`](@ref) performs two key tasks:
1.  **Writes a log file**: It saves a complete record of every parameter *actually used* in the simulation to `log_file`.
2.  **Performs sanity checks**: It verifies that all parameters in your override file were used.

The log file includes a `used_in` field, which lists every component that requested the parameter. Continuing the example, the log file would contain:

```TOML
[molar_mass_dry_air]
value = 0.03
type = "float"
description = "Molar mass of dry air (kg mol⁻¹)."
used_in = ["Thermodynamics"]
```

!!! note "Reproducibility"
    The generated log file is a valid TOML parameter file and can be used as an `override_file` to exactly reproduce an experiment.

!!! warn "Unused Parameter Checks"
    By default, [`log_parameter_information`](@ref) will issue a warning if any parameter in your override file was not requested by any component. To treat this as a fatal error, set its argument `strict=true`.