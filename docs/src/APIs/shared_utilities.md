<!-- # Shared Utilities

```@meta
CurrentModule = ClimaLand
```
## Domains

```@docs
ClimaLand.Domains.AbstractDomain
ClimaLand.Domains.SphericalShell
ClimaLand.Domains.SphericalSurface
ClimaLand.Domains.HybridBox
ClimaLand.Domains.Column
ClimaLand.Domains.Plane
ClimaLand.Domains.Point
ClimaLand.Domains.coordinates
ClimaLand.Domains.obtain_face_space
ClimaLand.Domains.obtain_surface_space
ClimaLand.Domains.obtain_surface_domain
ClimaLand.Domains.top_center_to_surface
ClimaLand.Domains.top_face_to_surface
ClimaLand.Domains.linear_interpolation_to_surface!
ClimaLand.Domains.get_Δz

```

## Models

```@docs
ClimaLand.AbstractModel
ClimaLand.AbstractImExModel
ClimaLand.AbstractExpModel
ClimaLand.make_exp_tendency
ClimaLand.make_imp_tendency
ClimaLand.make_compute_exp_tendency
ClimaLand.make_compute_imp_tendency
ClimaLand.make_update_aux
ClimaLand.make_update_boundary_fluxes
ClimaLand.make_set_initial_cache
ClimaLand.make_update_drivers
ClimaLand.prognostic_vars
ClimaLand.prognostic_types
ClimaLand.prognostic_domain_names
ClimaLand.auxiliary_vars
ClimaLand.auxiliary_types
ClimaLand.auxiliary_domain_names
ClimaLand.initialize_prognostic
ClimaLand.initialize_auxiliary
ClimaLand.initialize
ClimaLand.name
ClimaLand.AbstractBC
ClimaLand.AbstractSource
ClimaLand.source!
ClimaLand.AbstractBoundary
ClimaLand.TopBoundary
ClimaLand.BottomBoundary
ClimaLand.boundary_flux!
ClimaLand.diffusive_flux
ClimaLand.boundary_vars
ClimaLand.boundary_var_domain_names
ClimaLand.boundary_var_types
ClimaLand.make_jacobian
ClimaLand.make_compute_jacobian
ClimaLand.set_dfluxBCdY!
ClimaLand.get_drivers
```

## Drivers
```@docs
ClimaLand.PrescribedAtmosphere
ClimaLand.PrescribedPrecipitation
ClimaLand.PrescribedRadiativeFluxes
ClimaLand.PrescribedSoilOrganicCarbon
ClimaLand.CoupledAtmosphere
ClimaLand.CoupledRadiativeFluxes
ClimaLand.AbstractAtmosphericDrivers
ClimaLand.AbstractRadiativeDrivers
ClimaLand.turbulent_fluxes!
ClimaLand.turbulent_fluxes_at_a_point
ClimaLand.set_atmos_ts!
ClimaLand.surface_air_density
ClimaLand.surface_temperature
ClimaLand.surface_resistance
ClimaLand.surface_specific_humidity
``` -->
