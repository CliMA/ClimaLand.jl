Snow Model

```@meta
CurrentModule = ClimaLSM.Snow
```
## Snow Parameters

```@docs
ClimaLSM.Snow.SnowParameters
```

## Snow Functions of State

```@docs
ClimaLSM.Snow.specific_heat_capacity
ClimaLSM.Snow.snow_surface_temperature
ClimaLSM.Snow.snow_depth
ClimaLSM.Snow.snow_thermal_conductivity
ClimaLSM.Snow.snow_bulk_temperature
ClimaLSM.Snow.snow_liquid_mass_fraction
ClimaLSM.Snow.maximum_liquid_mass_fraction
ClimaLSM.Snow.runoff_timescale
```

## Tools for AI-Enabled Snow Forecasting

```@docs
ClimaLSM.Snow.DataTools.sitedata_daily
ClimaLSM.Snow.DataTools.sitedata_hourly
ClimaLSM.Snow.DataTools.apply_bounds
ClimaLSM.Snow.DataTools.hourly2daily
ClimaLSM.Snow.DataTools.rectify_daily_hourly
ClimaLSM.Snow.DataTools.scale_cols
ClimaLSM.Snow.DataTools.makediffs
ClimaLSM.Snow.DataTools.rolldata
ClimaLSM.Snow.DataTools.snowsplit
ClimaLSM.Snow.DataTools.filter_phys!
ClimaLSM.Snow.DataTools.prep_data
ClimaLSM.Snow.DataTools.make_data
ClimaLSM.Snow.DataTools.snotel_metadata
ClimaLSM.Snow.ModelTools.make_model
ClimaLSM.Snow.ModelTools.get_model_ps
ClimaLSM.Snow.ModelTools.settimescale!
ClimaLSM.Snow.ModelTools.setoutscale!
ClimaLSM.Snow.ModelTools.LinearModel
ClimaLSM.Snow.ModelTools.make_timeseries
ClimaLSM.Snow.ModelTools.trainmodel!
```