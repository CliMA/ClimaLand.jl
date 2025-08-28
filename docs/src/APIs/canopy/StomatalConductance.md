# Stomatal Conductance

```@meta
CurrentModule = ClimaLand.Canopy
```

## Models and Parameters

```@docs
ClimaLand.Canopy.MedlynConductanceModel
ClimaLand.Canopy.MedlynConductanceModel{FT}(
    domain;
    g1 = clm_medlyn_g1(domain.space.surface),
    g0::FT = LP.get_default_parameter(FT, :min_stomatal_conductance),
) where {FT <: AbstractFloat}
ClimaLand.Canopy.MedlynConductanceParameters
ClimaLand.Canopy.MedlynConductanceParameters(::Type{FT}; kwargs...) where {FT <: AbstractFloat}
ClimaLand.Canopy.PModelConductance
ClimaLand.Canopy.PModelConductance{FT}(; Drel = FT(1.6)) where {FT <: AbstractFloat}
ClimaLand.Canopy.PModelConductanceParameters
```

## Methods

```@docs
ClimaLand.Canopy.medlyn_term
ClimaLand.Canopy.medlyn_conductance
ClimaLand.Canopy.upscale_leaf_conductance
ClimaLand.Canopy.penman_monteith
```
