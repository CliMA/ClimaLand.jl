

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
    ClimaLSM.auxiliary_types(::AbstractCanopyComponent)

Returns the auxiliary types of the canopy component
passed in as an argument.
"""
ClimaLSM.auxiliary_types(::AbstractCanopyComponent) = ()

"""
    ClimaLSM.prognostic_types(::AbstractCanopyComponent)

Returns the prognostic types of the canopy component
passed in as an argument.
"""
ClimaLSM.prognostic_types(::AbstractCanopyComponent) = ()

"""
    ClimaLSM.prognostic_vars(::AbstractCanopyComponent)

Returns the prognostic vars of the canopy component
passed in as an argument.
"""
ClimaLSM.prognostic_vars(::AbstractCanopyComponent) = ()

"""
   prognostic_domain_names(m::AbstractCanopyComponent)

Returns the domain names for the prognostic variables in the form of a tuple.
"""
ClimaLSM.prognostic_domain_names(::AbstractCanopyComponent) = ()

"""
    ClimaLSM.auxiliary_vars(::AbstractCanopyComponent)

Returns the auxiliary types of the canopy component
passed in as an argument.
"""
ClimaLSM.auxiliary_vars(::AbstractCanopyComponent) = ()

"""
   auxiliary_domain_names(m::AbstractCanopyComponent)

Returns the domain names for the auxiliary variables in the form of a tuple.
"""
ClimaLSM.auxiliary_domain_names(::AbstractCanopyComponent) = ()

"""
    initialize_prognostic(
        component::AbstractCanopyComponent,
        state,
    )

Creates and returns a ClimaCore.Fields.FieldVector
with the prognostic variables of the canopy component
 `component`, stored using the name of the component.

The input `state` is usually a ClimaCore Field object.
"""
function initialize_prognostic(component::AbstractCanopyComponent, state)
    ClimaLSM.initialize_vars(
        prognostic_vars(component),
        prognostic_types(component),
        prognostic_domain_names(component),
        state,
        name(component),
    )
end

"""
    initialize_auxiliary(
        component::AbstractCanopyComponent,
        state,
    )

Creates and returns a ClimaCore.Fields.FieldVector
with the auxiliary variables of the canopy component
 `component`, stored using the name of the component.

The input `state` is usually a ClimaCore Field object.
"""
function initialize_auxiliary(
    component::AbstractCanopyComponent{FT},
    state,
) where {FT}
    ClimaLSM.initialize_vars(
        auxiliary_vars(component),
        auxiliary_types(component),
        auxiliary_domain_names(component),
        state,
        name(component),
    )
end

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


"""
     set_canopy_prescribed_field!(component::AbstractCanopyComponent,
                                  p,
                                  t0,
                                 ) end

Sets the spatially and temporally varying prescribed fields of the `component`
with their initial values.

These fields are stored in the aux-state and *should not* depend on the prognostic
state `Y` or other diagnostic variables stored in `p`; this allows them
to be updated first, prior to updating the rest of the aux state and prognostic state.

However, there is no
guarantee on the order of operations in terms of when diagnostic auxiliary
variables are updated vs. prescribed field auxiliary variables; any required
order of operations must be enforced by the developer who writes the update_aux
function.
"""
function set_canopy_prescribed_field!(component::AbstractCanopyComponent, p, t0) end

"""
     update_canopy_prescribed_field!(component::AbstractCanopyComponent,
                                     p,
                                     t,
                                     ) end
 
Updates the spatially and temporally varying prescribed fields of the `component`
with their values at time `t`.

These fields are stored in the aux-state and *should not* depend on the prognostic
state `Y` or other diagnostic variables stored in `p`; this allows them
to be updated first, prior to updating the rest of the aux state and prognostic state.

However, there is no
guarantee on the order of operations in terms of when diagnostic auxiliary
variables are updated vs. prescribed field auxiliary variables; any required
order of operations must be enforced by the developer who writes the update_aux
function.
"""
function update_canopy_prescribed_field!(
    component::AbstractCanopyComponent,
    p,
    t,
) end
