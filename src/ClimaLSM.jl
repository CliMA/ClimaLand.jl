module ClimaLSM
using UnPack
using DocStringExtensions
using ClimaCore
import ClimaCore: Fields
include("SharedUtilities/Domains.jl")
using .Domains
include("SharedUtilities/models.jl")
include("Bucket/Bucket.jl")

"""
     AbstractLandModel{FT} <: AbstractModel{FT} 

An abstract type for all land model types, which are used
to simulated multiple land surface components as
a single system. Standalone component runs do not require
this interface and it should not be used for that purpose.

Many methods taking an argument of type `AbstractLandModel` are
extensions of functions defined for `AbstractModel`s. 
There are default methods that apply for all `AbstractLandModel`s,
including `make_update_aux`, `make_ode_function, `make_rhs`, 
`initialize_prognostic`, `initialize_auxiliary`, `initialize`,
and `coordinates`.

Methods which dispatch on a specific type of AbstractLandModel
include any function involving interactions between components,
as these interactions depend on the components in the land model 
and the versions of these component models being used.
"""
abstract type AbstractLandModel{FT} <: AbstractModel{FT} end

ClimaLSM.name(::AbstractLandModel) = :land

function Domains.coordinates(model::AbstractLandModel)
    components = land_components(model)
    coords_list = map(components) do (component)
        Domains.coordinates(getproperty(model, component))
    end
    coords =
        ClimaCore.Fields.FieldVector(; NamedTuple{components}(coords_list)...)
    return coords
end

function initialize_prognostic(model::AbstractLandModel{FT}, coords) where {FT}
    components = propertynames(coords)
    Y_state_list = map(components) do (component)
        zero_state = map(_ -> zero(FT), getproperty(coords, component))
        getproperty(
            initialize_prognostic(getproperty(model, component), zero_state),
            component,
        )
    end
    Y = ClimaCore.Fields.FieldVector(; NamedTuple{components}(Y_state_list)...)
    return Y
end

function initialize_auxiliary(model::AbstractLandModel{FT}, coords) where {FT}
    components = propertynames(coords)
    p_state_list = map(components) do (component)
        zero_state = map(_ -> zero(FT), getproperty(coords, component))
        getproperty(
            initialize_auxiliary(getproperty(model, component), zero_state),
            component,
        )
    end
    p_interactions = initialize_interactions(model, coords)
    p = ClimaCore.Fields.FieldVector(;
        p_interactions...,
        NamedTuple{components}(p_state_list)...,
    )
    return p
end

function make_update_aux(land::AbstractLandModel)
    interactions_update_aux! = make_interactions_update_aux(land)
    components = land_components(land)
    update_aux_function_list =
        map(x -> make_update_aux(getproperty(land, x)), components)
    function update_aux!(p, Y, t)
        for f! in update_aux_function_list
            f!(p, Y, t)
        end
        interactions_update_aux!(p, Y, t) # this has to come last if it uses p.component.value!!
    end


    return update_aux!
end

function make_ode_function(land::AbstractLandModel)
    components = land_components(land)
    rhs_function_list = map(x -> make_rhs(getproperty(land, x)), components)
    update_aux! = make_update_aux(land)
    function ode_function!(dY, Y, p, t)
        update_aux!(p, Y, t)
        for f! in rhs_function_list
            f!(dY, Y, p, t)
        end
    end
    return ode_function!
end

"""
    land_components(land::AbstractLandModel)

Returns the component names of the `land` model, by calling
`propertynames(land)`.
"""
land_components(land::AbstractLandModel) = propertynames(land)

function prognostic_vars(land::AbstractLandModel)
    components = land_components(land)
    prognostic_list = map(components) do model
        prognostic_vars(getproperty(land, model))
    end
    return NamedTuple{components}(prognostic_list)
end

function prognostic_types(land::AbstractLandModel)
    components = land_components(land)
    prognostic_list = map(components) do model
        prognostic_types(getproperty(land, model))
    end
    return NamedTuple{components}(prognostic_list)
end

function auxiliary_vars(land::AbstractLandModel)
    components = land_components(land)
    auxiliary_list = map(components) do model
        auxiliary_vars(getproperty(land, model))
    end
    return NamedTuple{(components..., :interactions)}((
        auxiliary_list...,
        interaction_vars(land),
    ))
end

function auxiliary_types(land::AbstractLandModel)
    components = land_components(land)
    auxiliary_list = map(components) do model
        auxiliary_types(getproperty(land, model))
    end
    return NamedTuple{(components..., :interactions)}((
        auxiliary_list...,
        interaction_types(land),
    ))
end

"""
   interaction_vars(m::AbstractModel)

Returns the interaction variable symbols for the model in the form of a tuple.
"""
interaction_vars(m::AbstractLandModel) = ()

"""
   interaction_types(m::AbstractModel)

Returns the shared interaction variable types for the model in the form of a tuple.
"""
interaction_types(m::AbstractLandModel) = ()

"""
   interaction_domains(m::AbstractModel)

Returns the interaction domain symbols in the form of a tuple.

For example, for a boundary interaction term, one would specify a model with a
domain consisting of the boundary space. Most commonly, these are surface models.
This is only required for variables shared between land submodels, and only needed
for multi-component models, not standalone components. Component-specific variables
should be listed as prognostic or auxiliary variables.
"""
interaction_domains(m::AbstractLandModel) = ()

"""
    initialize_interactions(land::AbstractLandModel) end

Initializes interaction variables, which are a type of auxiliary
variable, to empty objects of the correct type for the model. 

Interaction variables are specified by `interaction_vars`, their types
by `interaction_types`, and their domains by `interaction_domains`. 
This function should be called during `initialize_auxiliary` step.
"""
function initialize_interactions(land::AbstractLandModel, land_coords)
    vars = interaction_vars(land)
    types = interaction_types(land)
    domains = interaction_domains(land)
    interactions = map(zip(types, domains)) do (T, domain)
        zero_instance = zero(T)
        map(_ -> zero_instance, getproperty(land_coords, domain))
    end

    return NamedTuple{vars}(interactions)
end
"""
    make_interactions_update_aux(land::AbstractLandModel) end

Makes and returns a function which updates the interaction variables, 
which are a type of auxiliary variable.

The `update_aux!` function returned is evaluted during the right hand
side evaluation.

This is a stub which concrete types of LSMs extend.
"""
function make_interactions_update_aux(land::AbstractLandModel) end

# Methods extended by the LSM models we support
include("SurfaceWater/Pond.jl")
using .Pond
import .Pond: surface_runoff
include("Soil/Soil.jl")
using .Soil
import .Soil: source!, boundary_fluxes
include("Vegetation/Roots.jl")
using .Roots
import .Roots: flow_out_roots

### Concrete types of AbstractLandModels
### and associated methods
include("./root_soil_model.jl")
include("./pond_soil_model.jl")
end
