export soil_surface_infiltration

"""
    AbstractRunoffModel

The abstract type for soil runoff models to be used with
the following boundary condition types:
- `ClimaLSM.Soil.AtmosDrivenFluxBC`
- `ClimaLSM.Soil.RichardsAtmosDrivenFluxBC`,

and for these functions:
-`ClimaLSM.Soil.soil_surface_infiltration`
- `ClimaLSM.Soil.subsurface_runoff_source`
- `ClimaLSM.source!`.

Please see the documentation for these for more details.
The model should specify the subsurface runoff sink term as well
as the surface runoff implementation.
"""
abstract type AbstractRunoffModel end

"""
    NoRunoff <: AbstractRunoffModel

A concrete type of soil runoff model; the 
default choice which does not include the 
effects of runoff.
"""
struct NoRunoff <: AbstractRunoffModel end

"""
    soil_surface_infiltration(::NoRunoff, net_water_flux, _...)

A function which computes the infiltration into the soil
 for the default of `NoRunoff`.

If `net_water_flux = P+E`, where `P` is the precipitation and `E`
is the evaporation (both negative if towards the soil), 
this returns `P+E` as the water boundary flux for the soil.
"""
soil_surface_infiltration(::NoRunoff, net_water_flux, _...) = net_water_flux

"""
    subsurface_runoff_source(runoff::AbstractRunoffModel)::Union{Nothing, AbstractSoilSource} 

A function which returns the soil source for the runoff model 
`runoff`; the default returns nothing in which case no source is added.

"""
subsurface_runoff_source(
    runoff::AbstractRunoffModel,
)::Union{Nothing, AbstractSoilSource} = nothing

"""
    append_source(src::AbstractSoilSource, srcs::Tuple)::Tuple

Appends `src` to the tuple of sources `srcs` if `src` is of type `AbstractSoilSource`.
"""
append_source(src::AbstractSoilSource, srcs::Tuple)::Tuple = (srcs..., src)

"""
    append_source(src::Nothing , srcs::Tuple)::Tuple

Appends `src` to the tuple of sources `srcs` if `src` is of type `AbstractSoilSource`.
"""
append_source(src::Nothing, srcs::Tuple)::Tuple = srcs
