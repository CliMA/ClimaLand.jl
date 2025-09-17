# Radiative Transfer

```@meta
CurrentModule = ClimaLand.Canopy
```

## Models and Parameters

```@docs
ClimaLand.Canopy.TwoStreamModel
ClimaLand.Canopy.TwoStreamModel{FT}(
    domain,
    toml_dict::CP.ParamDict
) where {FT <: AbstractFloat}
ClimaLand.Canopy.BeerLambertModel
ClimaLand.Canopy.BeerLambertModel{FT}(
    domain,
    toml_dict::CP.ParamDict;
) where {FT <: AbstractFloat}
ClimaLand.Canopy.BeerLambertParameters
ClimaLand.Canopy.BeerLambertParameters(toml_dict::CP.ParamDict)
ClimaLand.Canopy.TwoStreamParameters
ClimaLand.Canopy.TwoStreamParameters(toml_dict::CP.ParamDict)
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
