# Radiative Transfer

```@meta
CurrentModule = ClimaLand.Canopy
```

## Models and Parameters

```@docs
ClimaLand.Canopy.TwoStreamModel
ClimaLand.Canopy.TwoStreamModel{FT}(
    domain;
    radiation_parameters = clm_canopy_radiation_parameters(
        domain.space.surface,
    ),
    ϵ_canopy::FT = LP.get_default_parameter(FT, :canopy_emissivity),
    n_layers::Int = 20,
) where {FT <: AbstractFloat}
ClimaLand.Canopy.TwoStreamParameters
ClimaLand.Canopy.BeerLambertModel
ClimaLand.Canopy.BeerLambertModel{FT}(
    domain;
    radiation_parameters = clm_canopy_radiation_parameters(
        domain.space.surface,
    ),
    ϵ_canopy::FT = LP.get_default_parameter(FT, :canopy_emissivity),
) where {FT <: AbstractFloat}
ClimaLand.Canopy.BeerLambertParameters
ClimaLand.Canopy.BeerLambertParameters(::Type{FT}; kwargs...) where {FT <: AbstractFloat}
ClimaLand.Canopy.TwoStreamParameters(::Type{FT}; kwargs...) where {FT <: AbstractFloat}
```

## Radiative Transfer Parameterizations
```@docs
ClimaLand.Canopy.ConstantGFunction
ClimaLand.Canopy.CLMGFunction
```

## Methods

```@docs
ClimaLand.Canopy.canopy_radiant_energy_fluxes!
ClimaLand.Canopy.ground_albedo_PAR
ClimaLand.Canopy.ground_albedo_NIR
ClimaLand.Canopy.compute_fractional_absorbances!
ClimaLand.Canopy.canopy_sw_rt_beer_lambert
ClimaLand.Canopy.canopy_sw_rt_two_stream
ClimaLand.Canopy.extinction_coeff
ClimaLand.Canopy.compute_G
```
