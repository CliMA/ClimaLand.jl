module ClimaLSM

using UnPack
using ClimaCore
import ClimaCore: Fields
include("Domains.jl")
using .Domains
import .Domains: coordinates
include("Configurations.jl")
using .Configurations

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

struct NotIncluded{FT}<: AbstractModel{FT}
    model_name::Symbol
end



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
initialize(model::NotIncluded{FT}) where {FT} = (;model.model_name => FT[]), (;model.model_name => FT[]), (;model.model_name => FT[])

Domains.coordinates(model::AbstractModel) = Domains.coordinates(model.domain)
#### LandModel Specific
include("Soil.jl")
using .Soil
include("Roots.jl")
using .Roots
"""
    struct LandModel{FT, SM <: AbstractModel{FT}, RM <: AbstractModel{FT}} <: AbstractModel{FT}

A concrete type of `AbstractModel` for use in land surface modeling. Each component model of the
`LandModel` is itself an `AbstractModel`.

If a user wants to run in standalone, would they use this interface?
No, but should it work ?
"""
struct LandModel{
    FT,
    SM <: AbstractModel{FT},
    VM <: AbstractModel{FT},
    RM <: AbstractModel{FT},
    SNM <: AbstractModel{FT},
    CM <: AbstractModel{FT},
} <: AbstractModel{FT}
    soil::SM
    vegetation::VM
    river::RM
    snow::SNM
    carbon::CM
end

function LandModel{FT}(;
    soil_args::NamedTuple = (;),
    vegetation_args::NamedTuple = (;),
    river_args::NamedTuple = (;),
    snow_args::NamedTuple = (;),
    carbon_args::NamedTuple = (;),
) where {FT}

    configuration = LSMConfiguration{FT}()
    soil =
        isempty(soil_args) ? NotIncluded{FT}(:soil) :
        Soil.RichardsModel{FT}(; configuration = configuration, soil_args...)
    vegetation =
        isempty(vegetation_args) ? NotIncluded{FT}(:vegetation) :
        Roots.RootsModel{FT}(;
            configuration = configuration,
            vegetation_args...,
        )
    river = NotIncluded{FT}(:river)
    snow = NotIncluded{FT}(:snow)
    carbon = NotIncluded{FT}(:carbon)
    args = (soil, vegetation, river, snow, carbon)
    return LandModel{FT, typeof.(args)...}(args...)
end

auxiliary_vars(
    land::LandModel{FT, RichardsModel{FT}, RootsModel{FT}},
) where {FT} = (:root_extraction,)
#=
function initialize_interactions(
    land::LandModel{{FT, RichardsModel{FT}, RootsModel{FT}}
                    Y::ClimaCore.Fields.FieldVector
                    ) where {FT}
p.exchanges.variable?
    return 
=#
function initialize(land::LandModel)
    Y_soil, p_soil, coords_soil = initialize(land.soil)
    Y_vegetation, p_vegetation, coords_vegetation = initialize(land.vegetation)
    Y_river, p_river, coords_river = initialize(land.river)
    Y_snow, p_snow, coords_snow = initialize(land.snow)
    Y_carbon, p_carbon, coords_carbon = initialize(land.carbon)
    #p_interactions  = initialize_interactions(land{})
    # Do we just hardwire all of the interactions? This is so hardcoded, how would
    # we do this for other combos of components?
    # Each interaction will have a different domain as well
    # do these need coordinates??
    Y = ClimaCore.Fields.FieldVector(;
                                     soil = Y_soil.soil,
                                     vegetation = Y_vegetation.vegetation,
                                     river = Y_river.river,
                                     snow = Y_snow.snow,
                                     carbon = Y_carbon.carbon
    )
    p = ClimaCore.Fields.FieldVector(;
                                     #root_extraction = similar(Y_soil.soil.ϑ_l),
                                     soil = p_soil.soil,
                                     vegetation = p_vegetation.vegetation,
                                     river = p_river.river,
                                     snow = p_snow.snow,
                                     carbon = p_carbon.carbon,
    )
    coords =
        ClimaCore.Fields.FieldVector(; soil = coords_soil,
                                     vegetation = coords_vegetation,
                                     river = coords_river,
                                     snow = coords_snow,
                                     carbon = coords_carbon)
    return Y, p, coords
end

function make_update_aux(land::LandModel)
    interactions_update_aux! = make_interactions_update_aux(land)
    soil_update_aux! = make_update_aux(land.soil)
    vegetation_update_aux! = make_update_aux(land.vegetation)
    river_update_aux! = make_update_aux(land.river)
    snow_update_aux! = make_update_aux(land.snow)
    carbon_update_aux! = make_update_aux(land.carbon)
    function update_aux!(p, Y, t)
        interactions_update_aux!(p, Y, t)
        soil_update_aux!(p, Y, t)
        vegetation_update_aux!(p, Y, t)
        river_update_aux!(p, Y, t)
        snow_update_aux!(p, Y, t)
        carbon_update_aux!(p, Y, t)
    end
    return update_aux!
end

function make_interactions_update_aux(land::LandModel{FT}) where {FT}
    function update_aux!(p, Y, t)

       # @. p.root_extraction = FT(0.0)
        #=
        root_params = land.vegetation.param_set
        z = coordinates(land.soil)
        z_up = land.vegetation.domain.compartment_heights[1]
        rhog_MPa = FT(0.0098)
        @unpack a_root,
        b_root,
        K_max_root_moles,
        size_reservoir_stem_moles =  land.vegetation.param_set
        ρm = FT(1e6/18) # moles/m^3
        p_soil = p.soil.ψ .*rhog_MPa
        p_stem = theta_to_p(Y.vegetation.rwc[1] / size_reservoir_stem_moles)
        # computing root extraction as if there was a root in each layer
        # multiply by mask (P(root in that layer)?)
        # shouldnt do each step
        mask = zeros(length(parent(z)))
        map(x -> mask[argmin(parent(abs.(z .- x)))] = 1.0, z_root_depths)
        @. p.root_extraction = compute_flow(z, z_up, p_soil, p_stem, a_root, b_root, K_max_root_moles) / ρm # m^3/s need to convert to θ̇...

        ## we should: flow/ρm/Δz* (bio count of vegetation per area = M_r/M_R) OR
        ## flow/ρm *n(z) -> number density of vegetation per unit volume (z)
        =#
    end
    return update_aux!
end


function make_ode_function(land::LandModel)
    rhs_soil! = make_rhs(land.soil)
    rhs_vegetation! = make_rhs(land.vegetation)
    rhs_river! = make_rhs(land.river)
    rhs_carbon! = make_rhs(land.carbon)
    rhs_snow! = make_rhs(land.snow)
    
    update_aux! = make_update_aux(land)
    function ode_function!(dY, Y, p, t)
        update_aux!(p, Y, t)
        rhs_soil!(dY, Y, p, t)
        rhs_vegetation!(dY, Y, p, t)
        rhs_river!(dY, Y, p, t)
        rhs_snow!(dY, Y, p, t)
        rhs_carbon!(dY, Y, p, t)
    end
    return ode_function!
end

end # module
