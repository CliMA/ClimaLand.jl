module Configurations

"""

    AbstractConfiguration{FT <: AbstractFloat}

An abstract type for land surface and component simulation configurations.

Configuration types for each component hold the necessary information for
setting up the right hand side, including boundary conditions, source terms,
and forcing terms. 

This is currently not used, but may be a convenient shorthand way of
signifying a certain setup for different land models.
"""
abstract type AbstractConfiguration{FT <: AbstractFloat} end

end
