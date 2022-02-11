module ClimaLSM
using UnPack
using ClimaCore
import ClimaCore: Fields
include("Domains.jl")
using .Domains
include("models.jl")
include("Soil.jl")
using .Soil
import .Soil: source
include("Roots.jl")
using .Roots
import .Roots: flow_out_roots

export RootSoilModel

abstract type AbstractLandModel{FT} <: AbstractModel{FT} end

"""
    struct RootSoilModel{FT, SM <: AbstractModel{FT}, RM <: AbstractModel{FT}} <: AbstractLandModel{FT}

A concrete type of `AbstractModel` for use in land surface modeling. Each component model of the
`LandModel` is itself an `AbstractModel`.

If a user wants to run in standalone, would they use this interface?
No, but should it work ?
"""
struct RootSoilModel{
    FT,
    SM <: AbstractSoilModel{FT},
    VM <: AbstractVegetationModel{FT},
} <: AbstractLandModel{FT}
    soil::SM
    vegetation::VM
end

function RootSoilModel{FT}(;
    soil_model_type::Type{SM},
    soil_args::NamedTuple = (;),
    vegetation_model_type::Type{VM},
    vegetation_args::NamedTuple = (;),
) where {FT, SM <: AbstractSoilModel{FT}, VM <: AbstractVegetationModel{FT}}
    sources = (RootExtraction{FT}(),)
    boundary_fluxes = FluxBC{FT}(FT(0.0), FT(0.0))
    transpiration = PrescribedTranspiration{FT}((t::FT) -> FT(0.0))
    root_extraction = PrognosticSoilPressure{FT}()
    soil = soil_model_type(;
        boundary_conditions = boundary_fluxes,
        sources = sources,
        soil_args...,
    )
    vegetation = vegetation_model_type(;
        root_extraction = root_extraction,
        transpiration = transpiration,
        vegetation_args...,
    )
    args = (soil, vegetation)
    return RootSoilModel{FT, typeof.(args)...}(args...)
end

function initialize_interactions(
    land::RootSoilModel{FT, SM, RM},
) where {FT, SM <: Soil.RichardsModel{FT}, RM <: Roots.RootsModel{FT}}

    soil_coords = land.soil.coordinates
    return (root_extraction = similar(soil_coords),)
end

land_components(land::RootSoilModel) = (:soil, :vegetation)

function initialize(land::AbstractLandModel)
    components = land_components(land)
    coords_list = map(components) do (component)
        Domains.coordinates(getproperty(land, component))
    end
    Y_state_list = map(zip(components, coords_list)) do (component, coords)
        getproperty(
            initialize_prognostic(getproperty(land, component), coords),
            component,
        )
    end
    p_state_list = map(zip(components, coords_list)) do (component, coords)
        getproperty(
            initialize_auxiliary(getproperty(land, component), coords),
            component,
        )
    end
    p_interactions = initialize_interactions(land)

    Y = ClimaCore.Fields.FieldVector(;
        NamedTuple(zip(components, Y_state_list))...,
    )
    p = ClimaCore.Fields.FieldVector(;
        p_interactions...,
        NamedTuple(zip(components, p_state_list))...,
    )
    coords = ClimaCore.Fields.FieldVector(;
        NamedTuple(zip(components, coords_list))...,
    )
    return Y, p, coords
end

function make_update_aux(land::AbstractLandModel)
    interactions_update_aux! = make_interactions_update_aux(land)
    components = land_components(land)
    update_aux_function_list =
        map(x -> make_update_aux(getproperty(land, x)), components)
    function update_aux!(p, Y, t)
        interactions_update_aux!(p, Y, t)
        for f! in update_aux_function_list
            f!(p, Y, t)
        end
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

function make_interactions_update_aux(
    land::RootSoilModel{FT, SM, RM},
) where {FT, SM <: Soil.RichardsModel{FT}, RM <: Roots.RootsModel{FT}}
    function update_aux!(p, Y, t)
        @. p.root_extraction = FT(0.0)
        ##Science goes here
    end
    return update_aux!
end


## Extending methods of the Roots and Soil Model
struct PrognosticSoilPressure{FT} <: Roots.AbstractRootExtraction{FT} end

function Roots.flow_out_roots(
    re::PrognosticSoilPressure{FT},
    model::Roots.RootsModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
)::FT where {FT}
    return sum(p.root_extraction)
end


struct RootExtraction{FT} <: Soil.AbstractSoilSource{FT} end

function Soil.source(src::RootExtraction{FT}, Y, p) where {FT}
    V_layer = FT(1.0 * 0.15) # hardcoded
    ρm = FT(1e6 / 18) # moles/m^3
    return p.root_extraction ./ ρm ./ V_layer # Field
end

end # module
