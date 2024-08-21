# Soil Models

```@meta
CurrentModule = ClimaLand.Soil
```
## Soil Models

```@docs
ClimaLand.Soil.AbstractSoilModel
ClimaLand.Soil.RichardsModel
ClimaLand.Soil.EnergyHydrology
```
## Soil Parameter Structs

```@docs
ClimaLand.Soil.RichardsParameters
ClimaLand.Soil.EnergyHydrologyParameters
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
ClimaLand.Soil.NoRunoff
ClimaLand.Soil.SurfaceRunoff
ClimaLand.Soil.TOPMODELRunoff
ClimaLand.Soil.TOPMODELSubsurfaceRunoff
ClimaLand.Soil.subsurface_runoff_source
ClimaLand.Soil.update_runoff!
```

## Soil BC Methods and Types

```@docs
ClimaLand.Soil.MoistureStateBC
ClimaLand.Soil.HeatFluxBC
ClimaLand.Soil.WaterFluxBC
ClimaLand.Soil.TemperatureStateBC
ClimaLand.Soil.FreeDrainage
ClimaLand.Soil.RichardsAtmosDrivenFluxBC
ClimaLand.Soil.AtmosDrivenFluxBC
ClimaLand.Soil.WaterHeatBC
ClimaLand.Soil.soil_boundary_fluxes!
```

## Soil Source Types

```@docs
ClimaLand.Soil.AbstractSoilSource
ClimaLand.Soil.PhaseChange
ClimaLand.Soil.RootExtraction
```

## Soil Jacobian Structures

```@docs
ClimaLand.Soil.ImplicitEquationJacobian
```
