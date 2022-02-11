module Configurations

"""

    AbstractConfiguration{FT <: AbstractFloat}

An abstract type for land surface and component simulation configurations.

Configuration types for each component hold the necessary information for
setting up the right hand side, including boundary conditions, source terms,
and forcing terms. 

Each component has the freedom to define the attributes
of the configuration for that component as it makes sense for that 
component. In general, these attributes will be of a type which can be used
for multiple dispatch, such that the RHS functions call different methods
as appropriate.

For example, a `apply_boundary_conditions(bc)` function, appearing in
the right hand side of a discretized PDE, will execute different methods if
the boundary condition `bc` is of type `Dirichlet` or `Neumann`. A 
`compute_source(source)` function will execute different methods
depending on the type of the `source` argument.

"""
abstract type AbstractConfiguration{FT <: AbstractFloat} end

end
