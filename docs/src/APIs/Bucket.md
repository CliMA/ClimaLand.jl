# Bucket

```@meta
CurrentModule = ClimaLSM.Bucket
```
## Types

```@docs
ClimaLSM.Bucket.BucketModelParameters
ClimaLSM.Bucket.PrescribedAtmosphere
ClimaLSM.Bucket.PrescribedRadiativeFluxes
ClimaLSM.Bucket.AbstractAtmosphericDrivers
ClimaLSM.Bucket.AbstractRadiativeDrivers
ClimaLSM.Bucket.BulkAlbedoMap
ClimaLSM.Bucket.BulkAlbedoFunction
ClimaLSM.Bucket.BucketModel
```

## Extendible Methods

```@docs
ClimaLSM.Bucket.surface_fluxes
ClimaLSM.Bucket.net_radiation
ClimaLSM.Bucket.surface_air_density
ClimaLSM.Bucket.liquid_precipitation
ClimaLSM.Bucket.snow_precipitation
ClimaLSM.Bucket.surface_albedo
ClimaLSM.Bucket.beta_factor
```

## Methods for handling parameters read in from NetCDF files

```@docs
ClimaLSM.Bucket.set_initial_parameter_field!
ClimaLSM.Bucket.update_soil_albedo
```

## Artifact Path Functions

```@docs
ClimaLSM.Bucket.cesm2_land_albedo_dataset_path
```
