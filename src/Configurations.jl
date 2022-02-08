module Configurations

export AbstractConfiguration, RootSoilConfiguration
"""

    AbstractComponentExchange{FT <: AbstractFloat}

A general abstract type for any type of exchange between components of the 
Land Surface Model, either in standalone or integrated modes.

This encompasses true boundary conditions (e.g. Dirichlet or Neumann, 
at the boundary of the soil domain), source/sink terms (e.g. transpiration 
of leaves within the canopy airspace), and more general flows (e.g. a prescribed
soil pressure leading to root extraction of water from soil).

This is to be used both for standalone component runs, in which case exchanges 
are computed using only the component state and the prescribed quantities driving
the system, or for Land Surface Models, in which case exchanges are computed using
the entire land surface state and atmospheric quantities, which drive the system.

The user must define a concrete type for dispatching on
when computing necessary exchange/flux/flow quantities.

"""
abstract type AbstractConfiguration{FT <: AbstractFloat} end

"""
    RootSoilConfiguration{FT} <: AbstractConfiguration{FT}

Root-soil configuration. 
"""
Base.@kwdef struct RootSoilConfiguration{FT} <: AbstractConfiguration{FT}
    "Time dependent transpiration, given in moles/sec"
    T::Function = (t) -> FT(0.0)
    "Time dependent precipitation, given in m/s"
    P::Function = (t) -> FT(0.0)

end
end
