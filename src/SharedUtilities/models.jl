export AbstractModel,
    make_rhs,
    make_ode_function,
    make_update_aux,
    initialize_prognostic,
    initialize_auxiliary,
    initialize_vars,
    initialize,
    prognostic_vars,
    auxiliary_vars,
    prognostic_types,
    auxiliary_types,
    make_set_initial_aux_state,
    name

import .Domains: coordinates
## Default methods for all models - to be in a seperate module at some point.
"""
    abstract type AbstractModel{FT <: AbstractFloat}

An abstract type for all models.
"""
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
   prognostic_types(m::AbstractModel{FT}) where {FT}

Returns the prognostic variable types for the model in the form of a tuple.

Types provided must have `zero(T::DataType)` defined. Common examples
 include
- Float64, Float32 for scalar variables (a scalar value at each 
coordinate point)
- SVector{k,Float64} for a mutable but statically sized array of
 length `k` at each coordinate point.

Here, the coordinate points are those returned by coordinates(model).
"""
prognostic_types(m::AbstractModel) = ()
"""
   auxiliary_vars(m::AbstractModel)

Returns the auxiliary variable symbols for the model in the form of a tuple.
"""
auxiliary_vars(m::AbstractModel) = ()


"""
   auxiliary_types(m::AbstractModel{FT}) where {FT}

Returns the auxiliary variable types for the model in the form of a tuple.

Types provided must have `zero(T::DataType)` defined. Common examples
 include
- Float64, Float32 for scalar variables (a scalar value at each 
coordinate point)
- SVector{k,Float64} for a mutable but statically sized array of
 length `k` at each coordinate point.
- Note that Arrays, MVectors are not isbits and cannot be used.

Here, the coordinate points are those returned by coordinates(model).
"""
auxiliary_types(m::AbstractModel) = ()

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
    make_set_initial_aux_state(model::AbstractModel)

Returns the set_initial_aux_state! function, which updates the auxiliary
state `p` in place with the initial values corresponding to Y(t=t0) = Y0.

In principle, this function is not needed, because in the very first evaluation
of the `ode_function`, at t=t0, the auxiliary state is updated using the initial
conditions for Y=Y0. However,
without setting the initial `p` state prior to running the simulation,
 the value of `p` in the saved output at t=t0 will be unset.

Furthermore, specific methods of this function may be useful for models
which store time indepedent spatially varying parameter fields in 
the auxiliary state. In this case, `update_aux!` does not need to do
anything, but they do need to be set with the initial (constant) values
before the simulation can be carried out.
"""
function make_set_initial_aux_state(model::AbstractModel)
    update_aux! = make_update_aux(model)
    function set_initial_aux_state!(p, Y0, t0)
        update_aux!(p, Y0, t0)
    end
    return set_initial_aux_state!
end


"""
    initialize_prognostic(model::AbstractModel, state::Union{ClimaCore.Fields.Field, Vector{FT}})

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
    state::Union{ClimaCore.Fields.Field, Vector{FT}},
) where {FT}
    initialize_vars(
        prognostic_vars(model),
        prognostic_types(model),
        state,
        name(model),
    )
end

"""
    initialize_auxiliary(model::AbstractModel,state::Union{ClimaCore.Fields.Field, Vector{FT}})

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
    state::Union{ClimaCore.Fields.Field, Vector{FT}},
) where {FT}
    initialize_vars(
        auxiliary_vars(model),
        auxiliary_types(model),
        state,
        name(model),
    )
end

function initialize_vars(keys, types, state, model_name)
    FT = eltype(state)
    if length(keys) == 0
        return ClimaCore.Fields.FieldVector(; model_name => FT[])
    else
        zero_states = map(types) do (T)
            zero_instance = zero(T)
            map(_ -> zero_instance, state)
        end
        return ClimaCore.Fields.FieldVector(;
            model_name => (; zip(keys, zero_states)...),
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

    # Q: do we need separate initialize_prognostic and initialize_auxiliary methods?
    # A: yes - code is the same other than call to auxiliary_vars or prognostic_vars
    # however we do need to build separate FieldVectors => shared initialize_vars() method?
    # auxiliary_types, auxiliary_spaces, prognostic_types, prognostic_spaces
    Y = initialize_prognostic(model, coords)
    p = initialize_auxiliary(model, coords)
    return Y, p, coords
end
