# Insolation.jl

```@meta
CurrentModule = Insolation
```

`Insolation.jl` is a library that calculates the zenith angle and insolation 
at a given point point in space-time on an arbitrary planet.
The spatial location is defined by a (lat, lon) and the date specified with a 
DateTime object.
The planet is defined by the orbital parameters 
(obliquity, eccentricity, and longitude of perihelion).
The library is split between two files, `ZenithAngleCalc.jl` 
which calculates the zenith angle, azimuth angle, and planet-star distance, 
and `InsolationCalc.jl` which calculates the insolation given a zenith angle.

The zenith angle and insolation can both be calculated either as instantaneous 
values or as daily averaged values. The functions in `ZenithAngleCalc.jl` are 
overwritten to accept a variety of inputs, either calculating the orbital parameters 
given a DateTime object, or prescribing the orbital parameters as input.

The equations used for calculating orbital parameters are from Tapio Schneider's textbook draft. 
See [Zenith Angle Equations](@ref) for more details.

This package has been designed to be flexible for an arbitrary planet, but designed 
to calculate the insolation on Earth as the input for a climate model.

## Authors
`Insolation.jl` is being developed by [the Climate Modeling Alliance](https://clima.caltech.edu).
Specifically it has been developed to be used for the input to [`RRTMGP.jl`](https://github.com/CliMA/RRTMGP.jl) from [`ClimaAtmos.jl`](https://github.com/CliMA/ClimaAtmos.jl).