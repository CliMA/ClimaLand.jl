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
ClimaLSM.Domains.coordinates
ClimaLSM.Domains.obtain_face_space
ClimaLSM.Domains.obtain_surface_space
ClimaLSM.Domains.obtain_surface_domain
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
ClimaLSM.make_update_boundary_fluxes
ClimaLSM.make_set_initial_cache
ClimaLSM.make_update_drivers
ClimaLSM.prognostic_vars
ClimaLSM.prognostic_types
ClimaLSM.prognostic_domain_names
ClimaLSM.auxiliary_vars
ClimaLSM.auxiliary_types
ClimaLSM.auxiliary_domain_names
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
ClimaLSM.get_drivers
```

## Drivers
```@docs
ClimaLSM.PrescribedAtmosphere
ClimaLSM.PrescribedRadiativeFluxes
ClimaLSM.AbstractAtmosphericDrivers
ClimaLSM.AbstractRadiativeDrivers
ClimaLSM.turbulent_fluxes
ClimaLSM.turbulent_fluxes_at_a_point
ClimaLSM.radiative_fluxes_at_a_point
ClimaLSM.construct_atmos_ts
ClimaLSM.surface_air_density
ClimaLSM.liquid_precipitation
ClimaLSM.snow_precipitation
ClimaLSM.surface_temperature
ClimaLSM.surface_resistance
ClimaLSM.surface_specific_humidity
ClimaLSM.make_update_drivers
```
