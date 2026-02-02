# Soil Models

```@meta
CurrentModule = ClimaLand.Soil
```

## Soil Models

```@docs
ClimaLand.Soil.AbstractSoilModel
ClimaLand.Soil.RichardsModel
ClimaLand.Soil.RichardsModel{FT}(;
    parameters::RichardsParameters,
    domain::D,
    boundary_conditions::NamedTuple,
    sources::Tuple,
    lateral_flow::Bool = false,
) where {FT, D}
ClimaLand.Soil.RichardsModel{FT}(
    domain,
    forcing;
    runoff::Runoff.AbstractRunoffModel = Runoff.TOPMODELRunoff{FT}(
        f_over = FT(3.28), # extract from EPS
        R_sb = FT(1.484e-4 / 1000),# extract from EPS
        f_max = topmodel_fmax(domain.space.surface, FT),
    ),
    retention_parameters = soil_vangenuchten_parameters(
        domain.space.subsurface,
        FT,
    ),
    S_s = ClimaCore.Fields.zeros(domain.space.subsurface) .+ 1e-3,
) where {FT}
ClimaLand.Soil.EnergyHydrology
ClimaLand.Soil.EnergyHydrology{FT}(;
    parameters::EnergyHydrologyParameters{FT, PSE},
    domain::D,
    boundary_conditions::NamedTuple,
    sources::Tuple,
    lateral_flow::Bool = false,
) where {FT, D, PSE}
ClimaLand.Soil.EnergyHydrology{FT}(
    domain,
    forcing,
    toml_dict::CP.ParamDict
) where {FT <: AbstractFloat}
```

## Soil Parameter Structs

```@docs
ClimaLand.Soil.RichardsParameters
ClimaLand.Soil.RichardsParameters()
ClimaLand.Soil.EnergyHydrologyParameters
ClimaLand.Soil.EnergyHydrologyParameters(
    ::Type{FT};
    ν,
    ν_ss_om,
    ν_ss_quartz,
    ν_ss_gravel,
    hydrology_cm,
    K_sat,
    S_s,
    θ_r,
    albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(),
    kwargs...,
) where {FT <: AbstractFloat}
```

## Soil Hydrology Parameterizations

```@docs
ClimaLand.Soil.volumetric_liquid_fraction
ClimaLand.Soil.pressure_head
ClimaLand.Soil.hydraulic_conductivity
ClimaLand.Soil.impedance_factor
ClimaLand.Soil.viscosity_factor
ClimaLand.Soil.effective_saturation
ClimaLand.Soil.matric_potential
ClimaLand.Soil.dψdϑ
ClimaLand.Soil.inverse_matric_potential
ClimaLand.Soil.AbstractSoilHydrologyClosure
ClimaLand.Soil.vanGenuchten
ClimaLand.Soil.BrooksCorey
```

## Soil Heat Parameterizations

```@docs
ClimaLand.Soil.volumetric_heat_capacity
ClimaLand.Soil.κ_solid
ClimaLand.Soil.κ_sat_frozen
ClimaLand.Soil.κ_sat_unfrozen
ClimaLand.Soil.κ_sat
ClimaLand.Soil.κ_dry
ClimaLand.Soil.kersten_number
ClimaLand.Soil.relative_saturation
ClimaLand.Soil.volumetric_internal_energy
ClimaLand.Soil.volumetric_internal_energy_liq
ClimaLand.Soil.temperature_from_ρe_int
ClimaLand.Soil.thermal_conductivity
ClimaLand.Soil.phase_change_source
ClimaLand.Soil.thermal_time
```

## Soil Runoff Types and Methods

```@docs
ClimaLand.Soil.Runoff.AbstractRunoffModel
ClimaLand.Soil.NoRunoff
ClimaLand.Soil.SurfaceRunoff
ClimaLand.Soil.TOPMODELRunoff
ClimaLand.Soil.TOPMODELSubsurfaceRunoff
ClimaLand.Soil.subsurface_runoff_source
ClimaLand.Soil.update_infiltration_water_flux!
```

## Soil Albedo Types and Methods

```@docs
ClimaLand.Soil.CLMTwoBandSoilAlbedo
ClimaLand.Soil.ConstantTwoBandSoilAlbedo
ClimaLand.Soil.update_albedo!
```

## Soil BC Methods and Types

```@docs
ClimaLand.Soil.MoistureStateBC
ClimaLand.Soil.HeatFluxBC
ClimaLand.Soil.WaterFluxBC
ClimaLand.Soil.TemperatureStateBC
ClimaLand.Soil.FreeDrainage
ClimaLand.Soil.EnergyWaterFreeDrainage
ClimaLand.Soil.RichardsAtmosDrivenFluxBC
ClimaLand.Soil.AtmosDrivenFluxBC
ClimaLand.Soil.WaterHeatBC
ClimaLand.Soil.soil_boundary_fluxes!
```

## Soil Source Types

```@docs
ClimaLand.Soil.AbstractSoilSource
ClimaLand.Soil.PhaseChange
```
