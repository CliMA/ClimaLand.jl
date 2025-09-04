# Stomatal Conductance

```@meta
CurrentModule = ClimaLand.Canopy
```

## Models and Parameters

```@docs
ClimaLand.Canopy.MedlynConductanceModel
ClimaLand.Canopy.MedlynConductanceModel{FT}(
    domain,
    toml_dict::CP.ParamDict;
) where {FT <: AbstractFloat}
ClimaLand.Canopy.MedlynConductanceParameters
ClimaLand.Canopy.MedlynConductanceParameters(toml_dict::CP.ParamDict)
ClimaLand.Canopy.PModelConductance
ClimaLand.Canopy.PModelConductance{FT}(
    toml_dict::CP.ParamDict
) where {FT <: AbstractFloat}
ClimaLand.Canopy.PModelConductanceParameters
```

## Methods

```@docs
ClimaLand.Canopy.medlyn_term
ClimaLand.Canopy.medlyn_conductance
ClimaLand.Canopy.penman_monteith
```
