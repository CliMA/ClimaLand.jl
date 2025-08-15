# Parameters

## Changing parameter values

ClimaLand provides two ways to change the parameter values which are via the
keyword arguments of the function that create the model or parameter struct or
via the TOML files.

!!! warning "Keyword arguments over TOML files"
    Any keyword arguments take *precedence* over the parameters in the TOML
    files.

Parameters stored in TOML files are handled by `ClimaParams`[add link] which
also stores parameters in a single TOML file. (add link to this). This allows
for a single set of default parameters that can be used universally across
multiple experiments and setups.

!!! note "How should I change parameters?"
    We recommend that you override parameters with keyword arguments instead of
    changing the parameter values in the TOML files. In addition, we also
    recommend that you do not modify the TOML file and keyword arguments at the
    same time. This can easily lead to confusion, where you are not sure if the
    parameter in the TOML file or the value of the keyword argument is used in
    the simulation.

## List of TOML files

TODO: Make a list or table of the toml files and what each one of them is used for

# FAQ

## I added a new parameterization with new parameters. How do I add a parameter to a TOML file?

At the minimum, a parameter needs a parameter name, value, and type. For
example, it could look like this:

```
[example_land_parameter]
value = 42.0
type = "float"
description = "Represent an example land parameter in a new parameterization. Units are m/s."
tag = "NewParameterizationModel"
```

The possible types are `bool`, `float`, `integer`, or `string`. All parameter
names must be unique. The tag attribute which is  is optional, but it is
recommended to add one to make it easier for other users to know why a
particular parameter was added.

## What parameters are available?

All the parameters are stored in `toml` directory of `ClimaLand` (include link) and
in ClimaParams (include link).

## I changed a parameter in the TOML file(s) and it did not change the simulation. How can I check that a parameter is actually used?

`ClimaParams` provides the function `CP.log_parameter_information` which logs
all the parameters used in a simulation. Before calling
`ClimaLand.Simulations.solve!`, log the parameters, and inspect the TOML file
that is generated and see that the particular parameter you are interested in is
being used.
