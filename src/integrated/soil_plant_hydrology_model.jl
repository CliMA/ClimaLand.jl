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
canopy and a soil hydrology component.

This is primarily for testing purposes.
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

In this model, the soil_driver, of type PrognosticSoil, must be defined by the 
user and passed in within the land_args tuple. The soil_driver is used to
provide the canopy model with the soil albedo, and for the full soil + canopy 
model, the soil model contains the albedo, and the constructor for the model is 
able to define and pass the soil_driver, but here for the soil model, the soil 
albedo is not defined, and therefore the user must manually define the
soil driver and provide it in land_args.


Each component model is constructed with everything it needs to be stepped
forward in time, including boundary conditions, source terms, and interaction
terms.

Runoff is not currently supported for this configuration, as this model is for testing only;
 please see the documentation for the
`SoilCanopyModel`, `EnergyHydrology`, or `RichardsModel` for models that support runoff.
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
    (; soil_driver, atmos, radiation) = land_args
    precipitation = atmos.liquid_precip
    # These should always be set by the constructor.
    sources = (RootExtraction{FT}(),)

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
            transpiration = transpiration,
            canopy_component_args.hydraulics...,
        ),
        soil_driver = soil_driver,
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
is needed for the soil model and 
for the canopy model.

This function is called each ode function evaluation.

Root extraction is in units of 1/s and is equivalent to:
RAI * flux per cross section of roots * root distribution (z).
"""
function make_interactions_update_aux(
    land::SoilPlantHydrologyModel{FT, SM, RM},
) where {FT, SM <: Soil.RichardsModel{FT}, RM <: Canopy.CanopyModel{FT}}
    function update_aux!(p, Y, t)
        z = ClimaCore.Fields.coordinate_field(land.soil.domain.space).z
        (; conductivity_model) = land.canopy.hydraulics.parameters
        area_index = p.canopy.hydraulics.area_index
        @. p.root_extraction =
            (
                area_index.root + getproperty(
                    area_index,
                    land.canopy.hydraulics.compartment_labels[1],
                )
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

"""
    PlantHydraulics.root_flux_per_ground_area!(
        fa::ClimaCore.Fields.Field,
        s::PrognosticSoil,
        model::Canopy.PlantHydraulics.PlantHydraulicsModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    ) where {FT}

An extension of the `PlantHydraulics.root_flux_per_ground_area!` function,
 which returns the
net flux of water between the
roots and the soil, per unit ground area, 
when both soil and plant
hydraulics are modeled prognostically. This is for use in an LSM.

It is computed by summing the flux of water per ground area between
roots and soil at each soil layer.
"""
function PlantHydraulics.root_flux_per_ground_area!(
    fa::ClimaCore.Fields.Field,
    s::PrognosticSoil,
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
