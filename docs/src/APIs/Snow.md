Snow Model

```@meta
CurrentModule = ClimaLand.Snow
```
## Snow Model and Parameters

```@docs
ClimaLand.Snow.SnowModel
ClimaLand.Snow.SnowModel()
ClimaLand.Snow.SnowModel(
    FT,
    domain,
    forcing,
    toml_dict::CP.ParamDict,
    Δt;
)
ClimaLand.Snow.SnowParameters
ClimaLand.Snow.SnowParameters(::Type{FT}, Δt) where {FT <: AbstractFloat}
```

## Snow Functions of State

```@docs
ClimaLand.Snow.specific_heat_capacity
ClimaLand.Snow.snow_surface_temperature
ClimaLand.Snow.snow_thermal_conductivity
ClimaLand.Snow.snow_bulk_temperature
ClimaLand.Snow.maximum_liquid_mass_fraction
ClimaLand.Snow.runoff_timescale
ClimaLand.Snow.compute_water_runoff
ClimaLand.Snow.energy_from_q_l_and_swe
ClimaLand.Snow.energy_from_T_and_swe
ClimaLand.Snow.energy_flux_falling_rain
ClimaLand.Snow.energy_flux_falling_snow
```

## Computing fluxes for snow

```@docs
ClimaLand.Snow.snow_boundary_fluxes!
ClimaLand.Snow.phase_change_flux
ClimaLand.Snow.AtmosDrivenSnowBC
```

## Snow parameterizations

```@docs
ClimaLand.Snow.WuWuSnowCoverFractionModel
ClimaLand.Snow.ZenithAngleAlbedoModel
ClimaLand.Snow.ConstantAlbedoModel
```
