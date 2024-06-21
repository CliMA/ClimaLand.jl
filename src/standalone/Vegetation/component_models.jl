

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
methods for `ClimaLand.AbstractModel`s for
the canopy component parameterizations.
"""
abstract type AbstractCanopyComponent{FT <: AbstractFloat} end

"""
    ClimaLand.auxiliary_types(::AbstractCanopyComponent)

Returns the auxiliary types of the canopy component
passed in as an argument.
"""
ClimaLand.auxiliary_types(::AbstractCanopyComponent) = ()

"""
    ClimaLand.prognostic_types(::AbstractCanopyComponent)

Returns the prognostic types of the canopy component
passed in as an argument.
"""
ClimaLand.prognostic_types(::AbstractCanopyComponent) = ()

"""
    ClimaLand.prognostic_vars(::AbstractCanopyComponent)

Returns the prognostic vars of the canopy component
passed in as an argument.
"""
ClimaLand.prognostic_vars(::AbstractCanopyComponent) = ()

"""
   prognostic_domain_names(m::AbstractCanopyComponent)

Returns the domain names for the prognostic variables in the form of a tuple.
"""
ClimaLand.prognostic_domain_names(::AbstractCanopyComponent) = ()

"""
    ClimaLand.auxiliary_vars(::AbstractCanopyComponent)

Returns the auxiliary types of the canopy component
passed in as an argument.
"""
ClimaLand.auxiliary_vars(::AbstractCanopyComponent) = ()

"""
   auxiliary_domain_names(m::AbstractCanopyComponent)

Returns the domain names for the auxiliary variables in the form of a tuple.
"""
ClimaLand.auxiliary_domain_names(::AbstractCanopyComponent) = ()

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
    ClimaLand.initialize_vars(
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
    ClimaLand.initialize_vars(
        auxiliary_vars(component),
        auxiliary_types(component),
        auxiliary_domain_names(component),
        state,
        name(component),
    )
end

"""
     ClimaLand.make_compute_exp_tendency(component::AbstractCanopyComponent, canopy)

Creates the compute_exp_tendency!(dY,Y,p,t) function for the canopy `component`.

Since component models are not standalone models, other information
may be needed and passed in (via the `canopy` model itself).
The right hand side for the entire canopy model can make use of
these functions for the individual components.
"""
function ClimaLand.make_compute_exp_tendency(
    component::AbstractCanopyComponent,
    canopy,
)
    function compute_exp_tendency!(dY, Y, p, t)
        vars = prognostic_vars(component)
        dY_canopy = getproperty(dY, name(canopy))
        if !isempty(vars)
            getproperty(dY_canopy, name(component)) .= 0
        end
    end
    return compute_exp_tendency!
end

"""
     ClimaLand.make_compute_imp_tendency(component::AbstractCanopyComponent, canopy)

Creates the compute_imp_tendency!(dY,Y,p,t) function for the canopy `component`.

Since component models are not standalone models, other information
may be needed and passed in (via the `canopy` model itself).
The right hand side for the entire canopy model can make use of
these functions for the individual components.
"""
function ClimaLand.make_compute_imp_tendency(
    component::AbstractCanopyComponent,
    canopy,
)
    function compute_imp_tendency!(dY, Y, p, t)
        vars = prognostic_vars(component)
        dY_canopy = getproperty(dY, name(canopy))
        if !isempty(vars)
            getproperty(dY_canopy, name(component)) .= 0
        end

    end
    return compute_imp_tendency!
end


function make_compute_jacobian(component::AbstractCanopyComponent, canopy)
    function compute_jacobian!(W, Y, p, dtγ, t) end
    return compute_jacobian!
end

"""
     set_canopy_prescribed_field!(component::AbstractCanopyComponent,
                                  p,
                                  t,
                                 ) end

Updates the spatio-temporally varying prescribed fields of the `component`
with their values at time `t`.

These fields are stored in the aux-state, and currently are updated at the
beginning of the `update_aux!` function. Any required
order of operations must be enforced by the developer who writes the `update_aux!`
function.
"""
function set_canopy_prescribed_field!(component::AbstractCanopyComponent, p, t) end
