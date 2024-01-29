Snow Model

```@meta
CurrentModule = ClimaLand.Snow
```
## Snow Parameters

```@docs
ClimaLand.Snow.SnowParameters
```

## Snow Functions of State

```@docs
ClimaLand.Snow.specific_heat_capacity
ClimaLand.Snow.snow_surface_temperature
ClimaLand.Snow.snow_depth
ClimaLand.Snow.snow_thermal_conductivity
ClimaLand.Snow.snow_bulk_temperature
ClimaLand.Snow.snow_liquid_mass_fraction
ClimaLand.Snow.maximum_liquid_mass_fraction
ClimaLand.Snow.runoff_timescale
```

## Tools for AI-Enabled Snow Forecasting

```@docs
ClimaLand.Snow.DataTools.sitedata_daily
ClimaLand.Snow.DataTools.sitedata_hourly
ClimaLand.Snow.DataTools.apply_bounds
ClimaLand.Snow.DataTools.hourly2daily
ClimaLand.Snow.DataTools.rectify_daily_hourly
ClimaLand.Snow.DataTools.scale_cols
ClimaLand.Snow.DataTools.makediffs
ClimaLand.Snow.DataTools.rolldata
ClimaLand.Snow.DataTools.snowsplit
ClimaLand.Snow.DataTools.filter_phys!
ClimaLand.Snow.DataTools.prep_data
ClimaLand.Snow.DataTools.make_data
ClimaLand.Snow.DataTools.snotel_metadata
ClimaLand.Snow.ModelTools.make_model
ClimaLand.Snow.ModelTools.get_model_ps
ClimaLand.Snow.ModelTools.settimescale!
ClimaLand.Snow.ModelTools.setoutscale!
ClimaLand.Snow.ModelTools.LinearModel
ClimaLand.Snow.ModelTools.make_timeseries
ClimaLand.Snow.ModelTools.trainmodel!
```