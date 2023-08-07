export AbstractModel,
    AbstractImExModel,
    AbstractExpModel,
    make_imp_tendency,
    make_exp_tendency,
    make_compute_imp_tendency,
    make_compute_exp_tendency,
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
    name,
    domain_name

import .Domains: coordinates
## Default methods for all models - to be in a seperate module at some point.
"""
    abstract type AbstractModel{FT <: AbstractFloat}

An abstract type for all models.
"""
abstract type AbstractModel{FT <: AbstractFloat} end

"""
    AbstractImExModel{FT} <: AbstractModel{FT}

An abstract type for models which must be treated implicitly (and which may
also have tendency terms that can be treated explicitly).
This inherits all the default function definitions from AbstractModel, as well
as `make_imp_tendency` and `make_compute_imp_tendency` defaults.
"""
abstract type AbstractImExModel{FT} <: AbstractModel{FT} end

"""
    AbstractExpModel{FT} <: AbstractModel{FT}

An abstract type for models which must be treated explicitly.
This inherits all the default function definitions from AbstractModel, as well
as a `make_imp_tendency` default.
"""
abstract type AbstractExpModel{FT} <: AbstractModel{FT} end

"""
    name(model::AbstractModel)

Returns a symbol of the model component name, e.g. :soil or :vegetation.
"""
name(model::AbstractModel) =
    error("`name` not implemented for $(Base.typename(typeof(model)).wrapper)")

"""
    domain_name(model::AbstractModel)

Returns a symbol indicating the model's domain name, e.g. :surface or :subsurface. Only required for models that will be used as part of an LSM.
"""
domain_name(model::AbstractModel) = error(
    "`domain` not implemented for $(Base.typename(typeof(model)).wrapper)",
)



"""
   prognostic_vars(m::AbstractModel)

Returns the prognostic variable symbols for the model in the form of a tuple.
"""
prognostic_vars(m::AbstractModel) = ()

"""
   prognostic_types(m::AbstractModel{FT}) where {FT}

Returns the prognostic variable types for the model in the form of a tuple.

Types provided must have `ClimaCore.RecursiveApply.rzero(T::DataType)`
 defined. Common examples
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

Types provided must have `ClimaCore.RecursiveApply.rzero(T::DataType)`
defined. Common examples
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
    make_update_aux(model::AbstractModel)

Return an `update_aux!` function that updates auxiliary parameters `p`.
"""
function make_update_aux(model::AbstractModel)
    function update_aux!(p, Y, t) end
    return update_aux!
end

"""
    make_imp_tendency(model::AbstractImExModel)

Returns an `imp_tendency` that updates auxiliary variables and
updates the prognostic state of variables that are stepped implicitly.

`compute_imp_tendency!` should be compatible with OrdinaryDiffEq.jl solvers.
"""
function make_imp_tendency(model::AbstractImExModel)
    compute_imp_tendency! = make_compute_imp_tendency(model)
    update_aux! = make_update_aux(model)
    function imp_tendency!(dY, Y, p, t)
        update_aux!(p, Y, t)
        compute_imp_tendency!(dY, Y, p, t)
    end
    return imp_tendency!
end

"""
    make_imp_tendency(model::AbstractModel)

Returns an `imp_tendency` that does nothing. This model type is not
stepped explicity.
"""
function make_imp_tendency(model::AbstractModel)
    function imp_tendency!(dY, Y, p, t) end
end

"""
    make_exp_tendency(model::AbstractModel)

Returns an `exp_tendency` that updates auxiliary variables and
updates the prognostic state of variables that are stepped explicitly.

`compute_exp_tendency!` should be compatible with OrdinaryDiffEq.jl solvers.
"""
function make_exp_tendency(model::AbstractModel)
    compute_exp_tendency! = make_compute_exp_tendency(model)
    update_aux! = make_update_aux(model)
    function exp_tendency!(dY, Y, p, t)
        update_aux!(p, Y, t)
        compute_exp_tendency!(dY, Y, p, t)
    end
    return exp_tendency!
end

"""
    make_compute_imp_tendency(model::AbstractModel)

Return a `compute_imp_tendency!` function that updates state variables
that we will be stepped implicitly.

`compute_imp_tendency!` should be compatible with OrdinaryDiffEq.jl solvers.
"""
function make_compute_imp_tendency(model::AbstractModel)
    function compute_imp_tendency!(dY, Y, p, t) end
    return compute_imp_tendency!
end

"""
    make_compute_exp_tendency(model::AbstractModel)

Return a `compute_exp_tendency!` function that updates state variables
that we will be stepped explicitly.

`compute_exp_tendency!` should be compatible with OrdinaryDiffEq.jl solvers.
"""
function make_compute_exp_tendency(model::AbstractModel)
    function compute_exp_tendency!(dY, Y, p, t) end
    return compute_exp_tendency!
end

"""
    make_set_initial_aux_state(model::AbstractModel)

Returns the set_initial_aux_state! function, which updates the auxiliary
state `p` in place with the initial values corresponding to Y(t=t0) = Y0.

In principle, this function is not needed, because in the very first evaluation
of either `explicit_tendency` or `implicit_tendency`, at t=t0, the auxiliary
state is updated using the initial conditions for Y=Y0. However,
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

The input `state` is an array-like object, usually a
ClimaCore Field or a Vector{FT}.

Adjustments to this - for example because different prognostic variables
have different dimensions - require defining a new method.
"""
function initialize_prognostic(model::AbstractModel{FT}, state) where {FT}
    initialize_vars(
        prognostic_vars(model),
        prognostic_types(model),
        state,
        name(model),
    )
end

"""
    initialize_auxiliary(model::AbstractModel, state::Union{ClimaCore.Fields.Field, Vector{FT}})

Returns a FieldVector of auxiliary variables for `model` with the required
structure, with values equal to `similar(state)`. This assumes that all
 auxiliary variables are defined over the entire domain,
and that all auxiliary variables have the same dimension and type.

If a model has no auxiliary variables,
the returned FieldVector contains only an empty array.

The input `state` is an array-like object, usually a
ClimaCore Field or a Vector{FT}.

Adjustments to this - for example because different auxiliary variables
have different dimensions - require defining a new method.
"""
function initialize_auxiliary(model::AbstractModel{FT}, state) where {FT}
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
            zero_instance = ClimaCore.RecursiveApply.rzero(T)
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
    Y = initialize_prognostic(model, coords)
    p = initialize_auxiliary(model, coords)
    return Y, p, coords
end
