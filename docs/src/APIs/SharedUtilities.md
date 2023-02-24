# Shared Utilities

```@meta
CurrentModule = ClimaLSM
```
## Domains

```@docs
ClimaLSM.Domains.AbstractDomain
ClimaLSM.Domains.AbstractLSMDomain
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
ClimaLSM.make_set_initial_aux_state
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

## Drivers
```@docs
ClimaLSM.Drivers.PrescribedAtmosphere
ClimaLSM.Drivers.PrescribedRadiativeFluxes
ClimaLSM.Drivers.AbstractAtmosphericDrivers
ClimaLSM.Drivers.AbstractRadiativeDrivers
ClimaLSM.Drivers.surface_fluxes_at_a_point
ClimaLSM.Drivers.radiative_fluxes_at_a_point
ClimaLSM.Drivers.construct_atmos_ts
ClimaLSM.Drivers.surface_air_density
ClimaLSM.Drivers.liquid_precipitation
ClimaLSM.Drivers.snow_precipitation
