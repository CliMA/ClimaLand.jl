# Frequently Asked Questions

## Is it possible to preprocess the data in `TimeVaryingInput` or `SpaceVaryingInput`, for instance, to remove NaNs or change units?

Yes, [`TimeVaryingInput`](@ref timevaryinginput) and [`SpaceVaryingInput`](@ref spacevaryinginput) support this
feature. [`TimeVaryingInput`](@ref timevaryinginput) and [`SpaceVaryingInput`](@ref spacevaryinginput) that read
NetCDF files use [`NCFileReader`](@ref ncfilereaders) under the hood. `NCFileReader`s can be
constructed with an optional keyword argument `preprocess_func`, a pointwise
function that transforms the data read into something else. `Input`s can be
constructed to pass down this keyword argument. Let us have a look at an
example. Suppose the file `distances.nc` contains a time-varying variable
`distance` in centimeters, but we want it in meters.

Schematically,
```julia
import ClimaUtilities: TimeVaryingInputs
import ClimaCore
import Interpolations
# Loading ClimaCore, NCDatasets, Interpolations loads the extensions we need

# target_space is defined somewhere and is our computational grid
# start_date is the simulation date

distance_tv = TimeVaryingInputs.TimeVaryingInput("distances.nc",
                                                  "distance",
                                                  target_space;
                                                  start_date,
                                                  file_reader_kwargs = (; preprocess_func = x -> 10x))
```
