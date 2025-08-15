# Parameters

## How can I check that a parameter is actually being used?

`ClimaParams` provides the function [`ClimaParams.log_parameter_information`](https://clima.github.io/ClimaParams.jl/dev/API/#ClimaParams.log_parameter_information)
which logs the parameters used in a simulation. Before calling
`ClimaLand.Simulations.solve!`, log the parameters and inspect the generated
TOML file to verify that the particular parameter you're interested in is being
used.

!!! note "Keyword arguments"
    Note that any parameter overwritten by a keyword argument will not be
    logged.

```julia
# -- Set up simulation --
filepath_to_save_params = "parameters.log"
ClimaParams.log_parameter_information(toml_dict, filepath_to_save_params)
ClimaLand.Simulations.solve!(simulation)
```

!!! note "Logging parameters with versions above ClimaParams 0.12"
    With how `ClimaParams` is used with `ClimaLand`, the
    `ClimaParams.log_parameter_information` function requires ClimaParams
    version above 0.12 to behave correctly.

## Changing parameter values

`ClimaLand` provides two ways to change parameter values: via keyword arguments
of the functions that create the model or parameter struct, or via TOML files.

!!! warning "Keyword arguments over TOML files"
    Any keyword arguments take *precedence* over the parameters in the TOML
    files.

Parameters stored in TOML files are handled by
[`ClimaParams`](https://github.com/CliMA/ClimaParams.jl) which also stores
parameters in a
[single TOML file](https://github.com/CliMA/ClimaParams.jl/blob/main/src/parameters.toml).
This allows for a single set of default parameters that can be used universally
across multiple experiments and setups.

!!! note "How should I change parameters?"
    It depends on your use case! For example, if you are running calibration
    with `EnsembleKalmanProcesses.jl`, TOML dicts would be a natural choice, as
    they would be for any simulation using the default model. On the other hand,
    we use keyword arguments in many of our tutorials to avoid needing a
    different TOML file for each simulation. However, we recommend that you do
    not modify the TOML files and keyword arguments at the same time. This can
    easily lead to confusion about whether the parameter in the TOML file or the
    keyword argument value is being used in the simulation.

## How do I use the TOML files?

Most of the constructors for the land model take in a `toml_dict`, which is a
`ClimaParams.ParamDict`. `ClimaLand` provides the function `create_toml_dict`
for creating a `ClimaParams.ParamDict` from a TOML file. See the example below
of constructing one.

```julia
import ClimaLand.Parameters as LP

# Store parameters for running long runs
default_params_filepath =
    joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml")
toml_dict = LP.create_toml_dict(FT, default_params_filepath)

# Store parameters for running bucket model
bucket_params_filepath =
    joinpath(pkgdir(ClimaLand), "toml", "bucket_parameters.toml")
toml_dict = LP.create_toml_dict(FT, bucket_params_filepath)
```

# FAQ

## I added a new parameterization with new parameters. How do I add a parameter to a TOML file?

To add a parameter to the TOML file, a parameter needs a name, value, and type.
For example, it could look like this:

```
[example_land_parameter]
value = 42.0
type = "float"
description = "Represents an example land parameter in a new parameterization. Units are m/s."
tag = "NewParameterizationModel"
```

The possible types are `bool`, `float`, `integer`, or `string`. All parameter
names must be unique. The tag attribute is optional, but we recommend adding
a tag to make it easier for other users to know where a particular parameter is
used.

## What parameters are available?

All parameters are stored in:
- The `toml` [directory](https://github.com/CliMA/ClimaLand.jl/tree/main/toml)
  of `ClimaLand`
- The [TOML file](https://github.com/CliMA/ClimaParams.jl/blob/main/src/parameters.toml)
  in `ClimaParams`

## I changed a parameter in the TOML file(s) and it did not change the simulation. How can I check that a parameter is actually being used?

`ClimaParams` provides the function [`ClimaParams.log_parameter_information`](https://clima.github.io/ClimaParams.jl/dev/API/#ClimaParams.log_parameter_information)
which logs the parameters used in a simulation. See the [first section](#How-can-I-check-that-a-parameter-is-actually-being-used?) of this
page for more information.

!!! note "Keyword arguments"
    Note that any parameter overwritten by a keyword argument will not be
    logged.

# API

```@docs
ClimaLand.Parameters.LandParameters
ClimaLand.Parameters.create_toml_dict
```
