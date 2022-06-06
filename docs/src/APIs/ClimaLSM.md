# ClimaLSM

```@meta
CurrentModule = ClimaLSM
```
## LSM Model Types and methods

```@docs
ClimaLSM.RootSoilModel
ClimaLSM.LandHydrology
ClimaLSM.make_interactions_update_aux
ClimaLSM.initialize_interactions
ClimaLSM.land_components
ClimaLSM.interaction_vars
ClimaLSM.interaction_types
ClimaLSM.interaction_domains
ClimaLSM.domain
```

## Land Hydrology

```@docs
ClimaLSM.infiltration_capacity
ClimaLSM.infiltration_at_point
ClimaLSM.PrognosticRunoff
ClimaLSM.RunoffBC
```

## RootSoilModel

```@docs
ClimaLSM.PrognosticSoilPressure
ClimaLSM.RootExtraction
```

## BucketModel

```@docs
ClimaLSM.BucketModelParameters
ClimaLSM.PrescribedAtmosphere
ClimaLSM.PrescribedRadiation
ClimaLSM.AbstractAtmosphericDrivers
ClimaLSM.AbstractRadiativeDrivers
ClimaLSM.BulkAlbedo
ClimaLSM.BucketModel
ClimaLSM.surface_fluxes
ClimaLSM.surface_air_density
ClimaLSM.liquid_precipitation
ClimaLSM.surface_albedo