Snow Model

```@meta
CurrentModule = ClimaLand.Snow
```
## Snow Model and Parameters

```@docs
ClimaLand.Snow.SnowModel
ClimaLand.Snow.SnowModel(;
    parameters::SnowParameters{FT, DM, PSE},
    domain::ClimaLand.Domains.AbstractDomain,
    boundary_conditions::BC,
) where {FT, DM, PSE, BC}
ClimaLand.Snow.SnowModel(
    FT,
    domain,
    forcing,
    toml_dict::CP.AbstractTOMLDict,
    Δt;
    prognostic_land_components = (:snow,),
    z_0m = toml_dict["snow_momentum_roughness_length"],
    z_0b = toml_dict["snow_scalar_roughness_length"],
    ϵ_snow = toml_dict["snow_emissivity"],
    α_snow = ConstantAlbedoModel(toml_dict["snow_albedo"]),
    density = MinimumDensityModel(toml_dict["snow_density"]),
    scf = WuWuSnowCoverFractionModel(toml_dict, FT(1)),
    θ_r = toml_dict["holding_capacity_of_water_in_snow"],
    Ksat = toml_dict["wet_snow_hydraulic_conductivity"],
    ΔS = toml_dict["delta_S"],
)
ClimaLand.Snow.SnowParameters
ClimaLand.Snow.SnowParameters(::Type{FT}, Δt; kwargs...) where {FT <: AbstractFloat}
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
