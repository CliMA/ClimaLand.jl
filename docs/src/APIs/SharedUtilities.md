# Shared Utilities

```@meta
CurrentModule = ClimaLSM
```
## Domains

```@docs
ClimaLSM.Domains.AbstractDomain
ClimaLSM.Domains.RootDomain
ClimaLSM.Domains.AbstractVegetationDomain
ClimaLSM.Domains.HybridBox
ClimaLSM.Domains.Column
ClimaLSM.Domains.Plane
ClimaLSM.Domains.Point
ClimaLSM.Domains.LSMSingleColumnDomain
ClimaLSM.Domains.coordinates
```

## Models

```@docs
ClimaLSM.AbstractModel
ClimaLSM.make_ode_function
ClimaLSM.make_rhs
ClimaLSM.make_update_aux
ClimaLSM.prognostic_vars
ClimaLSM.auxiliary_vars
ClimaLSM.initialize_prognostic
ClimaLSM.initialize_auxiliary
ClimaLSM.initialize
ClimaLSM.name
```