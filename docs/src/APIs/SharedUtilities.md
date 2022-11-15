# Shared Utilities

```@meta
CurrentModule = ClimaLSM
```
## Domains

```@docs
ClimaLSM.Domains.AbstractDomain
ClimaLSM.Domains.AbstractLSMDomain
ClimaLSM.Domains.PlantHydraulicsDomain
ClimaLSM.Domains.AbstractVegetationDomain
ClimaLSM.Domains.SphericalShell
ClimaLSM.Domains.SphericalSurface
ClimaLSM.Domains.HybridBox
ClimaLSM.Domains.Column
ClimaLSM.Domains.Plane
ClimaLSM.Domains.Point
ClimaLSM.Domains.LSMSingleColumnDomain
ClimaLSM.Domains.LSMMultiColumnDomain
ClimaLSM.Domains.LSMSphericalShellDomain
ClimaLSM.Domains.coordinates
ClimaLSM.Domains.obtain_face_space
ClimaLSM.Domains.obtain_surface_space
```

## Models

```@docs
ClimaLSM.AbstractModel
ClimaLSM.make_ode_function
ClimaLSM.make_rhs
ClimaLSM.make_update_aux
ClimaLSM.prognostic_vars
ClimaLSM.prognostic_types
ClimaLSM.auxiliary_vars
ClimaLSM.auxiliary_types
ClimaLSM.initialize_prognostic
ClimaLSM.initialize_auxiliary
ClimaLSM.initialize
ClimaLSM.name
ClimaLSM.AbstractBC
ClimaLSM.AbstractSource
ClimaLSM.source!
ClimaLSM.AbstractBoundary
ClimaLSM.TopBoundary
ClimaLSM.BottomBoundary
ClimaLSM.boundary_flux
ClimaLSM.diffusive_flux
ClimaLSM.get_Δz
```