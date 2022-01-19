module ComponentExchanges

export AbstractComponentExchange
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
abstract type AbstractComponentExchange{FT <: AbstractFloat} end

end
