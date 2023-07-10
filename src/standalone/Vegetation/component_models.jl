

"""
    AbstractCanopyComponent{FT <: AbstractFloat}

An abstract type for canopy component parameterizations.

Canopy component parameterizations do not run in standalone
mode, but only as part of a `CanopyModel`. As such, they
do not require all of the functionality of
`AbstractModel`s,  and they are not `AbstractModel`s
themselves. The `CanopyModel` is an `AbstractModel`.

However, some of the same functionality is nice to have
for canopy components, especially when defining the variables,
which is why we introduce the
`AbstractCanopyComponent` type and extend many of the
methods for `ClimaLSM.AbstractModel`s for
the canopy component parameterizations.
"""
abstract type AbstractCanopyComponent{FT <: AbstractFloat} end

"""
    ClimaLSM.auxiliary_types(::AbstractCanopyComponent, canopy)

Returns the auxiliary types of the canopy component
passed in as an argument.
"""
ClimaLSM.auxiliary_types(::AbstractCanopyComponent, canopy) = ()

"""
    ClimaLSM.prognostic_types(::AbstractCanopyComponent, canopy)

Returns the prognostic types of the canopy component
passed in as an argument.
"""
ClimaLSM.prognostic_types(::AbstractCanopyComponent, canopy) = ()

"""
    ClimaLSM.prognostic_vars(::AbstractCanopyComponent)

Returns the prognostic vars of the canopy component
passed in as an argument.
"""
ClimaLSM.prognostic_vars(::AbstractCanopyComponent) = ()

"""
    ClimaLSM.auxiliary_vars(::AbstractCanopyComponent)

Returns the auxiliary types of the canopy component
passed in as an argument.
"""
ClimaLSM.auxiliary_vars(::AbstractCanopyComponent) = ()

"""
     ClimaLSM.make_compute_exp_tendency(component::AbstractCanopyComponent, canopy)

Creates the compute_exp_tendency!(dY,Y,p,t) function for the canopy `component`.

Since component models are not standalone models, other information
may be needed and passed in (via the `canopy` model itself).
The right hand side for the entire canopy model can make use of
these functions for the individual components.
"""
function ClimaLSM.make_compute_exp_tendency(::AbstractCanopyComponent, canopy)
    function compute_exp_tendency!(dY, Y, p, t) end
    return compute_exp_tendency!
end
