## ITime

`ITime`, or _integer time_, is a time type used by CliMA simulations to keep
track of simulation time. For more information, refer to the
[TimeManager section](https://clima.github.io/ClimaUtilities.jl/dev/timemanager/)
in ClimaUtilities and the 
[ITime section](https://clima.github.io/ClimaAtmos.jl/dev/itime/) in ClimaAtmos.

### How do I use ITime?

We do not support the functionality of running simulations using either `ITime` or
floating point. See the simulations below that automatically use `ITime`.

1. Snowy land long run
2. Soil-Canopy long run
3. California regional simulation long run
4. Soil model long run
5. Global bucket simulation long run
6. Bucket ERA5 model experiment
7. Richards runoff standalone experiment

If there is a particular simulation that is not compatible with `ITime` and you
would like to use `ITime` with it, then please open an issue for it!
