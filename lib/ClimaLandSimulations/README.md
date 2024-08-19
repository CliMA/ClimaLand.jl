# ClimaLandSimulations

A library of methods to run ClimaLand:
- at single sites 
- globally (WIP)

## Installation

```jl
julia> ] # to go into pkg mode
pkg> add https://github.com/CliMA/ClimaLand.jl:lib/ClimaLandSimulations # because unregistered for now
pkg> *backspace* # press backspace button to go back into julia mode
julia> using ClimaLandSimulations
```

For examples, see the [experiments](https://github.com/CliMA/ClimaLand.jl/lib/ClimaLandSimulations/experiments) folder

### Development of the `ClimaLandSimulations` subpackage

    cd ClimaLand/lib/ClimaLandSimulations

    # Add ClimaLand to subpackage environment
    julia --project -e 'using Pkg; Pkg.develop(path="../../")'

    # Instantiate ClimaLandSimulations project environment
    julia --project -e 'using Pkg; Pkg.instantiate()'
    julia --project -e 'using Pkg; Pkg.test()'
