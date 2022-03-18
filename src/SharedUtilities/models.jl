export AbstractModel,
    make_rhs,
    make_ode_function,
    make_update_aux,
    initialize_prognostic,
    initialize_auxiliary,
    initialize,
    prognostic_vars,
    auxiliary_vars,
    name

import .Domains: coordinates
## Default methods for all models - to be in a seperate module at some point.

abstract type AbstractModel{FT <: AbstractFloat} end

"""
    name(model::AbstractModel)

Returns a symbol of the model component name, e.g. :soil or :vegetation.
"""
name(model::AbstractModel) =
    error("`name` not implemented for $(Base.typename(typeof(model)).wrapper)")

"""
   prognostic_vars(m::AbstractModel)

Returns the prognostic variable symbols for the model in the form of a tuple.
"""
prognostic_vars(m::AbstractModel) = ()

"""
   auxiliary_vars(m::AbstractModel)

Returns the auxiliary variable symbols for the model in the form of a tuple.
"""
auxiliary_vars(m::AbstractModel) = ()

"""
    make_rhs(model::AbstractModel)

Return a `rhs!` function that updates state variables.

`rhs!` should be compatible with OrdinaryDiffEq.jl solvers.
"""
function make_rhs(model::AbstractModel)
    function rhs!(dY, Y, p, t) end
    return rhs!
end

"""
    make_update_aux(model::AbstractModel)

Return an `update_aux!` function that updates auxiliary parameters `p`.
"""
function make_update_aux(model::AbstractModel)
    function update_aux!(p, Y, t) end
    return update_aux!
end

"""
    make_ode_function(model::AbstractModel)

Returns an `ode_function` that updates auxiliary variables and
updates the prognostic state.

`ode_function!` should be compatible with OrdinaryDiffEq.jl solvers.
"""
function make_ode_function(model::AbstractModel)
    rhs! = make_rhs(model)
    update_aux! = make_update_aux(model)
    function ode_function!(dY, Y, p, t)
        update_aux!(p, Y, t)
        rhs!(dY, Y, p, t)
    end
    return ode_function!
end

"""
    initialize_prognostic(model::AbstractModel, state::Union{ClimaCore.Fields.Field, Vector{FT}, FT})

Returns a FieldVector of prognostic variables for `model` with the required
structure, with values equal to `similar(state)`. This assumes that all 
prognostic variables are defined over the entire domain,
and that all prognostic variables have the same dimension and type.

If a model has no prognostic variables,
the returned FieldVector contains only an empty array.

Adjustments to this - for example because different prognostic variables
have different dimensions - require defining a new method.
"""
function initialize_prognostic(
    model::AbstractModel{FT},
    state::Union{ClimaCore.Fields.Field, Vector{FT}, FT},
) where {FT}
    model_name = name(model)
    keys = prognostic_vars(model)
    if length(keys) == 0
        return ClimaCore.Fields.FieldVector(; model_name => FT[])
    else
        values = map((x) -> similar(state), keys)
        return ClimaCore.Fields.FieldVector(;
            model_name => (; zip(keys, values)...),
        )
    end

end

"""
    initialize_auxiliary(model::AbstractModel,state::Union{ClimaCore.Fields.Field, Vector{FT}, FT})

Returns a FieldVector of auxiliary variables for `model` with the required
structure, with values equal to `similar(state)`. This assumes that all
 auxiliary variables are defined over the entire domain,
and that all auxiliary variables have the same dimension and type.

If a model has no auxiliary variables,
the returned FieldVector contains only an empty array.

Adjustments to this - for example because different auxiliary variables
have different dimensions - require defining a new method.
"""
function initialize_auxiliary(
    model::AbstractModel{FT},
    state::Union{ClimaCore.Fields.Field, Vector{FT}, FT},
) where {FT}
    model_name = name(model)
    keys = auxiliary_vars(model)
    if length(keys) == 0
        return ClimaCore.Fields.FieldVector(; model_name => FT[])
    else
        values = map((x) -> similar(state), keys)
        return ClimaCore.Fields.FieldVector(;
            model_name => (; zip(keys, values)...),
        )
    end
end

Domains.coordinates(model::AbstractModel) = Domains.coordinates(model.domain)

"""
    initialize(model::AbstractModel)

Creates the prognostic and auxiliary states structures, but with unset 
values; constructs and returns the coordinates for the `model` domain.
We may need to consider this default more as we add diverse components and 
`Simulations`.
"""
function initialize(model::AbstractModel{FT}) where {FT}
    coords = Domains.coordinates(model)
    zero_state = map(_ -> zero(FT), coords)
    Y = initialize_prognostic(model, zero_state)
    p = initialize_auxiliary(model, zero_state)
    return Y, p, coords
end
