export RootSoilModel
"""
    struct RootSoilModel{
        FT,
        SM <: Soil.AbstractSoilModel{FT},
        VM <: Roots.AbstractVegetationModel{FT},
    } <: AbstractLandModel{FT}
        soil::SM
        vegetation::VM
    end

A concrete type of land model used for simulating systems with a 
vegetation and a soil component.
$(DocStringExtensions.FIELDS)
"""
struct RootSoilModel{
    FT,
    SM <: Soil.AbstractSoilModel{FT},
    VM <: Roots.AbstractVegetationModel{FT},
} <: AbstractLandModel{FT}
    "The soil model to be used"
    soil::SM
    "The vegetation model to be used"
    vegetation::VM
end


"""
    RootSoilModel{FT}(;
                      soil_model_type::Type{SM},
                      soil_args::NamedTuple = (;),
                      vegetation_model_type::Type{VM},
                      vegetation_args::NamedTuple = (;),
                      ) where {FT,
                               SM <: Soil.AbstractSoilModel{FT},
                               VM <: Roots.AbstractVegetationModel{FT}}

A constructor for the `RootSoilModel`, which takes in the concrete model
type and required arguments for each component, constructs those models,
and constructs the `RootSoilModel` from them.

Each component model is constructed with everything it needs to be stepped
forward in time, including boundary conditions, source terms, and interaction
terms.
"""
function RootSoilModel{FT}(;
    soil_model_type::Type{SM},
    soil_args::NamedTuple = (;),
    vegetation_model_type::Type{VM},
    vegetation_args::NamedTuple = (;),
) where {
    FT,
    SM <: Soil.AbstractSoilModel{FT},
    VM <: Roots.AbstractVegetationModel{FT},
}

    #These may be passed in, or set, depending on use scenario
    boundary_fluxes = FluxBC{FT}(FT(0.0), FT(0.0))
    transpiration = PrescribedTranspiration{FT}((t::FT) -> FT(0.0))

    ##These should always be set by the constructor.
    sources = (RootExtraction{FT}(),)
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

interaction_vars(m::RootSoilModel) = (:root_extraction,)

interaction_types(m::RootSoilModel{FT}) where {FT} = (FT,)

interaction_domains(m::RootSoilModel) = (:soil,)

"""
    make_interactions_update_aux(
        land::RootSoilModel{FT, SM, RM},
    ) where {FT, SM <: Soil.RichardsModel{FT}, RM <: Roots.RootsModel{FT}}

A method which makes a function; the returned function 
updates the auxiliary variable `p.root_extraction`, which
is needed for both the boundary condition for the soil model and the source
term (runoff) for the surface water model.

This function is called each ode function evaluation.
"""
function make_interactions_update_aux(#Do we want defaults, for land::AbstractLandModel?
    land::RootSoilModel{FT, SM, RM},
) where {FT, SM <: Soil.RichardsModel{FT}, RM <: Roots.RootsModel{FT}}
    function update_aux!(p, Y, t)
        @. p.root_extraction = FT(0.0)
        ##Science goes here
    end
    return update_aux!
end


## Extending methods of the Roots and Soil Model
## TBD if these should be linked more explicitly.
"""
   PrognosticSoilPressure{FT} <: Roots.AbstractRootExtraction{FT}

Concrete type of Roots.AbstractRootExtraction, used for dispatch 
in an LSM with both soil and plant hydraulic components.

This is paired with the source term `Soil.RootExtraction`:both 
are used at the same time,
ensuring that the water flow into the roots is extracted correctly
from the soil.
"""
struct PrognosticSoilPressure{FT} <: Roots.AbstractRootExtraction{FT} end

"""
    Roots.flow_out_roots(
        re::PrognosticSoilPressure{FT},
        model::Roots.RootsModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    )::FT where {FT}

An extension of the `Roots.flow_out_roots` function,
 which returns the
net flow of water between the
roots and the soil, when both soil and plant
hydraulics are modeled prognostically. This is for use in an LSM.

It is computed by summing the flow of water between
roots and soil at each soil layer.
"""
function Roots.flow_out_roots(
    re::PrognosticSoilPressure{FT},
    model::Roots.RootsModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
)::FT where {FT}
    return sum(p.root_extraction)
end

"""
    RootExtraction{FT} <: Soil.AbstractSoilSource{FT}

Concrete type of Soil.AbstractSoilSource, used for dispatch 
in an LSM with both soil and plant hydraulic components.

This is paired with the source term `Roots.PrognosticSoilPressure`:both 
are used at the same time,
ensuring that the water flow into the roots is extracted correctly
from the soil.
"""
struct RootExtraction{FT} <: Soil.AbstractSoilSource{FT} end

"""
    Soil.source!(dY::ClimaCore.Fields.FieldVector,
                          src::RootExtraction{FT},
                          Y::ClimaCore.Fields.FieldVector,
                          p::ClimaCore.Fields.FieldVector)::ClimaCore.Fields.Field  where {FT}

An extension of the `Soil.source!` function,
 which computes source terms for the 
soil model; this method returns the water loss or gain due
to roots when a plant hydraulic prognostic model is included.
"""
function Soil.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::RootExtraction{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
)::ClimaCore.Fields.Field where {FT}
    return dY.soil.ϑ_l .+= p.root_extraction
end
