# SurfaceFluxes.jl

```@docs
SurfaceFluxes
```

## Core input types

```@docs
SurfaceFluxes.StateValues
```

## Dispatch types

```@docs
SurfaceFluxes.Fluxes
SurfaceFluxes.FluxesAndFrictionVelocity
SurfaceFluxes.Coefficients
SurfaceFluxes.ValuesOnly
```

## User-facing methods

```@docs
SurfaceFluxes.surface_conditions
SurfaceFluxes.recover_profile
```

# Parameters
Convenience constructors are provided for the `SurfaceFluxesParameters` and the various `UniversalFunctions` parameter structs.
To use them, you must first import ClimaParams:
```julia
import ClimaParams as CP
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF

FT = Float64

# SurfaceFluxesParameters requires a float type and a UniversalFunctionsParameters type
SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

# Or a TOML dict instead of a float type
toml_dict = CP.create_toml_dict(Float64)
SFP.SurfaceFluxesParameters(toml_dict, UF.GrachevParams)

# UniversalFunctionsParameters only require a float type or a TOML dict.
UF.BusingerParams(FT)
UF.GryanikParams(FT)
UF.GrachevParams(FT)
```

## Universal Functions

```@docs
SurfaceFluxes.UniversalFunctions
```

```@docs
SurfaceFluxes.UniversalFunctions.Gryanik
SurfaceFluxes.UniversalFunctions.Grachev
SurfaceFluxes.UniversalFunctions.Businger
```

```@docs
SurfaceFluxes.UniversalFunctions.phi
SurfaceFluxes.UniversalFunctions.psi
SurfaceFluxes.UniversalFunctions.Psi
```
