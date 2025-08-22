# Canopy Energy Model

```@meta
CurrentModule = ClimaLand.Canopy
```

## Models and Parameters

```@docs
ClimaLand.Canopy.PrescribedCanopyTempModel
ClimaLand.Canopy.BigLeafEnergyModel
ClimaLand.Canopy.BigLeafEnergyModel{FT}(;
    ac_canopy::FT = FT(2e3),
) where {FT <: AbstractFloat}
ClimaLand.Canopy.BigLeafEnergyParameters
ClimaLand.Canopy.BigLeafEnergyParameters(
    toml_dict::CP.AbstractTOMLDict;
    ac_canopy = toml_dict["ac_canopy"],
)
ClimaLand.Canopy.AbstractCanopyEnergyModel
```

## Methods

```@docs
ClimaLand.Canopy.canopy_temperature
ClimaLand.Canopy.root_energy_flux_per_ground_area!
ClimaLand.Canopy.total_energy_per_area!
```
