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
ClimaLSM.Domains.top_center_to_surface
```

## Models

```@docs
ClimaLSM.AbstractModel
ClimaLSM.AbstractImExModel
ClimaLSM.AbstractExpModel
ClimaLSM.make_exp_tendency
ClimaLSM.make_imp_tendency
ClimaLSM.make_compute_exp_tendency
ClimaLSM.make_compute_imp_tendency
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
ClimaLSM.make_tendency_jacobian
ClimaLSM.make_update_jacobian
ClimaLSM.∂tendencyBC∂Y
ClimaLSM.AbstractTridiagonalW
```

## Drivers
```@docs
ClimaLSM.PrescribedAtmosphere
ClimaLSM.PrescribedRadiativeFluxes
ClimaLSM.AbstractAtmosphericDrivers
ClimaLSM.AbstractRadiativeDrivers
ClimaLSM.surface_fluxes_at_a_point
ClimaLSM.radiative_fluxes_at_a_point
ClimaLSM.construct_atmos_ts
ClimaLSM.surface_air_density
ClimaLSM.liquid_precipitation
ClimaLSM.snow_precipitation
```
