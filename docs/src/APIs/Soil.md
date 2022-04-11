# Soil Models

```@meta
CurrentModule = ClimaLSM.Soil
```
## Soil Models

```@docs
ClimaLSM.Soil.AbstractSoilModel
ClimaLSM.Soil.RichardsModel
ClimaLSM.Soil.EnergyHydrology
```

## Soil Functions of State

```@docs
ClimaLSM.Soil.volumetric_liquid_fraction
ClimaLSM.Soil.pressure_head
ClimaLSM.Soil.hydraulic_conductivity
ClimaLSM.Soil.impedance_factor
ClimaLSM.Soil.viscosity_factor
ClimaLSM.Soil.effective_saturation
ClimaLSM.Soil.matric_potential
ClimaLSM.Soil.volumetric_heat_capacity
ClimaLSM.Soil.κ_solid
ClimaLSM.Soil.κ_sat_frozen
ClimaLSM.Soil.κ_sat_unfrozen
ClimaLSM.Soil.κ_sat
ClimaLSM.Soil.κ_dry
ClimaLSM.Soil.kersten_number
ClimaLSM.Soil.relative_saturation
ClimaLSM.Soil.volumetric_internal_energy
ClimaLSM.Soil.volumetric_internal_energy_liq
ClimaLSM.Soil.temperature_from_ρe_int
ClimaLSM.Soil.thermal_conductivity
```

## Soil Parameters

```@docs
ClimaLSM.Soil.RichardsParameters
ClimaLSM.Soil.EnergyHydrologyParameters
```

## Soil Methods and Types

```@docs
ClimaLSM.Soil.boundary_fluxes
ClimaLSM.Soil.FluxBC
ClimaLSM.Soil.AbstractSoilSource
ClimaLSM.Soil.source!
```