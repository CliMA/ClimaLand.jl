export SoilPlantHydrologyModel
"""
    struct SoilPlantHydrologyModel{
        FT,
        SM <: Soil.AbstractSoilModel{FT},
        VM <: Canopy.CanopyModel{FT},
    } <: AbstractLandModel{FT}
        soil::SM
        canopy::VM
    end

A concrete type of land model used for simulating systems with a 
canopy and a soil component.
$(DocStringExtensions.FIELDS)
"""
struct SoilPlantHydrologyModel{
    FT,
    SM <: Soil.AbstractSoilModel{FT},
    VM <: Canopy.CanopyModel{FT},
} <: AbstractLandModel{FT}
    "The soil model to be used"
    soil::SM
    "The canopy model to be used"
    canopy::VM
end


"""
    SoilPlantHydrologyModel{FT}(;
                                 land_args::NamedTuple = (;),
                                 soil_model_type::Type{SM},
                                 soil_args::NamedTuple = (;),
                                 canopy_component_types::NamedTuple = (;),
                                 canopy_component_args::NamedTuple = (;),
                                 canopy_model_args::NamedTuple = (;),
                                 ) where {FT, SM <: Soil.AbstractSoilModel{FT}}
A constructor for the `SoilPlantHydrologyModel`, which takes in the concrete model
type and required arguments for each component, constructs those models,
and constructs the `SoilPlantHydrologyModel` from them.

Each component model is constructed with everything it needs to be stepped
forward in time, including boundary conditions, source terms, and interaction
terms.
"""
function SoilPlantHydrologyModel{FT}(;
    land_args::NamedTuple = (;),
    soil_model_type::Type{SM},
    soil_args::NamedTuple = (;),
    canopy_component_types::NamedTuple = (;),
    canopy_component_args::NamedTuple = (;),
    canopy_model_args::NamedTuple = (;),
) where {FT, SM <: Soil.AbstractSoilModel{FT}}

    # These may be passed in, or set, depending on use scenario.
    (; atmos, radiation) = land_args
    precipitation = atmos.liquid_precip
    # These should always be set by the constructor.
    sources = (RootExtraction{FT}(),)
    root_extraction = PrognosticSoilPressure{FT}()

    boundary_conditions = (;
        top = (water = FluxBC((p, t) -> eltype(t)(precipitation(t))),),
        bottom = (water = Soil.FreeDrainage(),),
    )

    soil = soil_model_type(;
        boundary_conditions = boundary_conditions,
        sources = sources,
        soil_args...,
    )

    if :transpiration in propertynames(land_args)
        transpiration = land_args.transpiration
    else
        transpiration = Canopy.PlantHydraulics.DiagnosticTranspiration{FT}()
    end

    canopy = Canopy.CanopyModel{FT}(;
        radiative_transfer = canopy_component_types.radiative_transfer(
            canopy_component_args.radiative_transfer...,
        ),
        photosynthesis = canopy_component_types.photosynthesis(
            canopy_component_args.photosynthesis...,
        ),
        conductance = canopy_component_types.conductance(
            canopy_component_args.conductance...,
        ),
        hydraulics = canopy_component_types.hydraulics(;
            root_extraction = root_extraction,
            transpiration = transpiration,
            canopy_component_args.hydraulics...,
        ),
        atmos = atmos,
        radiation = radiation,
        canopy_model_args...,
    )

    return SoilPlantHydrologyModel{FT, typeof(soil), typeof(canopy)}(
        soil,
        canopy,
    )
end

interaction_vars(m::SoilPlantHydrologyModel) = (:root_extraction,)

interaction_types(m::SoilPlantHydrologyModel{FT}) where {FT} = (FT,)

interaction_domains(m::SoilPlantHydrologyModel) = (:subsurface,)

"""
    make_interactions_update_aux(
        land::SoilPlantHydrologyModel{FT, SM, RM},
    ) where {FT, SM <: Soil.RichardsModel{FT}, RM <: Canopy.CanopyModel{FT}}

A method which makes a function; the returned function 
updates the auxiliary variable `p.root_extraction`, which
is needed for both the boundary condition for the soil model and the source
term (runoff) for the surface water model.

This function is called each ode function evaluation.

Root extraction is in units of 1/s and is equivalent to:
RAI * flux per cross section of roots * root distribution (z).
"""
function make_interactions_update_aux(
    land::SoilPlantHydrologyModel{FT, SM, RM},
) where {FT, SM <: Soil.RichardsModel{FT}, RM <: Canopy.CanopyModel{FT}}
    function update_aux!(p, Y, t)
        z = ClimaCore.Fields.coordinate_field(land.soil.domain.space).z
        (; area_index, conductivity_model) = land.canopy.hydraulics.parameters
        @. p.root_extraction =
            (
                area_index[:root] +
                area_index[land.canopy.hydraulics.compartment_labels[1]]
            ) / 2 *
            PlantHydraulics.flux(
                z,
                land.canopy.hydraulics.compartment_midpoints[1],
                p.soil.ψ,
                p.canopy.hydraulics.ψ[1],
                PlantHydraulics.hydraulic_conductivity(
                    conductivity_model,
                    p.soil.ψ,
                ),
                PlantHydraulics.hydraulic_conductivity(
                    conductivity_model,
                    p.canopy.hydraulics.ψ[1],
                ),
            ) *
            (land.canopy.hydraulics.parameters.root_distribution(z))
    end
    return update_aux!
end


## Extending methods of the Plant Hydraulics and Soil Model
## TBD if these should be linked more explicitly.
"""
   PrognosticSoilPressure{FT} <: PlantHydraulics.AbstractRootExtraction{FT}

Concrete type of PlantHydraulics.AbstractRootExtraction, used for dispatch 
in an LSM with both soil and plant hydraulic components.

This is paired with the source term `Soil.RootExtraction`:both 
are used at the same time,
ensuring that the water flux into the roots is extracted correctly
from the soil.
"""
struct PrognosticSoilPressure{FT} <:
       Canopy.PlantHydraulics.AbstractRootExtraction{FT} end

"""
    PlantHydraulics.root_flux_per_ground_area!(
        fa::ClimaCore.Fields.Field,
        re::PrognosticSoilPressure{FT},
        model::Canopy.PlantHydraulics.PlantHydraulicsModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    ) where {FT}

An extension of the `PlantHydraulics.root_flux_per_ground_area!` function,
 which returns the
net flux of water between the
roots and the soil, per unti ground area, 
when both soil and plant
hydraulics are modeled prognostically. This is for use in an LSM.

It is computed by summing the flux of water per ground area between
roots and soil at each soil layer.
"""
function PlantHydraulics.root_flux_per_ground_area!(
    fa::ClimaCore.Fields.Field,
    re::PrognosticSoilPressure{FT},
    model::Canopy.PlantHydraulics.PlantHydraulicsModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
) where {FT}
    fa .= sum(p.root_extraction)
end

"""
    RootExtraction{FT} <: Soil.AbstractSoilSource{FT}

Concrete type of Soil.AbstractSoilSource, used for dispatch 
in an LSM with both soil and plant hydraulic components.

This is paired with the source term `Canopy.PlantHydraulics.PrognosticSoilPressure`:both 
are used at the same time,
ensuring that the water flux into the roots is extracted correctly
from the soil.
"""
struct RootExtraction{FT} <: Soil.AbstractSoilSource{FT} end

"""
    ClimaLSM.source!(dY::ClimaCore.Fields.FieldVector,
                     src::RootExtraction,
                     Y::ClimaCore.Fields.FieldVector,
                     p::ClimaCore.Fields.FieldVector
                     model::RichardsModel)

An extension of the `ClimaLSM.source!` function,
 which computes source terms for the 
soil model; this method returns the water loss or gain due
to roots when a plant hydraulic prognostic model is included.
"""
function ClimaLSM.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::RootExtraction,
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    model::RichardsModel,
)
    @. dY.soil.ϑ_l += -1 * p.root_extraction
    # if flow is negative, towards soil -> soil water increases, add in sign here.
end
