# Getting Started

## Installation of Julia and ClimaLand

ClimaLand is provided as a Julia package, so it requires having Julia installed. Information about Julia packages is available on the [Julia website](https://julialang.org/packages/).

First, download and install Julia by following the instructions at [https://julialang.org/downloads/](https://julialang.org/downloads/).
Then, you can launch a [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) and install the
ClimaLand.jl package by running:

```julia
julia> using Pkg
julia> Pkg.add(ClimaLand)
julia> using ClimaLand
```

A typical land simulation employs several different parameterizations to model the various land-surface processes. Let's start our journet into ClimaLand by looking at one of those.

### Parameterization

Let's start with a basic example: compute canopy gross photosynthesis (GPP).

```@repl
using ClimaLand
@doc ClimaLand.Canopy.compute_GPP
```

As you can see, our parameterization for GPP is located in the [`Canopy`](@ref Photosynthesis) Module, and requires four arguments.
For example, with An = 5 µmol m⁻² s⁻¹, K = 0.5, LAI = 3 m² m⁻², Ω = 0.7, you can compute GPP like below:

```@repl
import ClimaLand.Canopy as canopy
canopy.compute_GPP(5.0, 0.5, 3.0, 0.7)
```

Et voilà!

Note that our package [ParamViz](https://github.com/CliMA/ParamViz.jl) allows interactive visualisation of
our parameterizations. See examples in the standalone models pages.

### ClimaLand structure

ClimaLand contains multiple modules. They are listed below:

```@repl
using MethodAnalysis, ClimaLand
child_modules(ClimaLand)
```

To explore what modules, functions and types are exported in a particular module, you can use [About.jl](https://github.com/tecosaur/About.jl):

```@repl
using ClimaLand
using About
about(ClimaLand.Soil.Biogeochemistry)
```

Where modules are shown in red, functions are shown in blue, and types are shown in yellow.

To see the documentation about a particular module, function or type, you can use ? to go in help mode
in the REPL, or `@doc` as in [Parameterization above](#Parameterization).
