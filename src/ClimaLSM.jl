module ClimaLSM

using UnPack
using ClimaCore
import ClimaCore: Fields
include("Domains.jl")
using .Domains
import .Domains: coordinates
include("ComponentExchanges.jl")
using .ComponentExchanges

export AbstractModel,
    make_rhs,
    make_ode_function,
    make_update_aux,
    initialize_prognostic,
    initialize_auxiliary,
    initialize,
    LandModel,
    prognostic_vars,
    auxiliary_vars
## Default methods for all models - to be in a seperate module at some point.

abstract type AbstractModel{FT <: AbstractFloat} end

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
structure, with values equal to `similar(state)`. This assumes that all prognostic
variables are defined over the entire domain. 

The default is a model with no prognostic variables, in which case 
the returned FieldVector contains only an empty array.

Adjustments to this - for example because different prognostic variables
have different dimensions - require defining a new method.
"""
function initialize_prognostic(
    model::AbstractModel{FT},
    state::Union{ClimaCore.Fields.Field, Vector{FT}, FT},
) where {FT}
    keys = prognostic_vars(model)
    if length(keys) == 0
        return ClimaCore.Fields.FieldVector(; model.model_name => FT[])
    else
        values = map((x) -> similar(state), keys)
        return ClimaCore.Fields.FieldVector(;
            model.model_name => (; zip(keys, values)...),
        )
    end

end

"""
    initialize_auxiliary(model::AbstractModel,state::Union{ClimaCore.Fields.Field, Vector{FT}, FT})

Returns a FieldVector of auxiliary variables for `model` with the required
structure, with values equal to `similar(state)`. This assumes that all auxiliary
variables are defined over the entire domain. 

The default is a model with no auxiliary variables, in which case 
the returned FieldVector contains only an empty array.

Adjustments to this - for example because different auxiliary variables
have different dimensions - require defining a new method.
"""
function initialize_auxiliary(
    model::AbstractModel{FT},
    state::Union{ClimaCore.Fields.Field, Vector{FT}, FT},
) where {FT}
    keys = auxiliary_vars(model)
    if length(keys) == 0
        return ClimaCore.Fields.FieldVector(; model.model_name => FT[])
    else
        values = map((x) -> similar(state), keys)
        return ClimaCore.Fields.FieldVector(;
            model.model_name => (; zip(keys, values)...),
        )
    end
end

"""
    initialize(model::AbstractModel)

Creates the prognostic and auxiliary states structures, but with unset values
not set; constructs and returns the coordinates for the `model` domain.
We may need to consider this default more as we add diverse components and 
`Simulations`.
"""
function initialize(model::AbstractModel)
    coords = coordinates(model)
    Y = initialize_prognostic(model, coords)
    p = initialize_auxiliary(model, coords)
    return Y, p, coords
end

Domains.coordinates(model::AbstractModel) = Domains.coordinates(model.domain)

#### LandModel Specific

"""
    struct LandModel{FT, SM <: AbstractModel{FT}, RM <: AbstractModel{FT}} <: AbstractModel{FT}

A concrete type of `AbstractModel` for use in land surface modeling. Each component model of the
`LandModel` is itself an `AbstractModel`.
"""
struct LandModel{FT, SM <: AbstractModel{FT}, RM <: AbstractModel{FT}} <:
       AbstractModel{FT}
    soil::SM
    roots::RM
end

function initialize(land::LandModel)
    Y_soil, p_soil, coords_soil = initialize(land.soil)
    Y_roots, p_roots, coords_roots = initialize(land.roots)
    Y = ClimaCore.Fields.FieldVector(;
        soil = Y_soil.soil,
        roots = Y_roots.roots,
    )
    p = ClimaCore.Fields.FieldVector(;
        soil = p_soil.soil,
        roots = p_roots.roots,
    )
    coords =
        ClimaCore.Fields.FieldVector(; soil = coords_soil, roots = coords_roots)
    return Y, p, coords
end

function make_update_aux(land::LandModel)
    soil_update_aux! = make_update_aux(land.soil)
    roots_update_aux! = make_update_aux(land.roots)
    function update_aux!(p, Y, t)
        soil_update_aux!(p, Y, t)
        roots_update_aux!(p, Y, t)
    end
    return update_aux!
end

function make_ode_function(land::LandModel)
    rhs_soil! = make_rhs(land.soil)
    rhs_roots! = make_rhs(land.roots)
    update_aux! = make_update_aux(land)
    function ode_function!(dY, Y, p, t)
        update_aux!(p, Y, t)
        rhs_soil!(dY.soil, Y.soil, p, t)
        rhs_roots!(dY.roots, Y.roots, p, t)
    end
    return ode_function!
end

include("Soil.jl")
include("Roots.jl")

end # module
