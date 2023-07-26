export AbstractSoilDriver, PrescribedSoil, PrognosticSoil

"""
An abstract type of soil drivers of the canopy model.
"""
abstract type AbstractSoilDriver{FT <: AbstractFloat} end

"""
     PrescribedSoil{FT} <: AbstractSoilDriver{FT}

A container for holding prescribed soil parameters needed by the canopy model
when running the canopy in standalone mode, including the soil pressure and 
albedo.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedSoil{FT} <: AbstractSoilDriver{FT}
    "The depth of the root tips, in meters"
    root_depths::Vector{FT}
    "Prescribed soil potential (m) as a function of time"
    ψ_soil::Function
    "Soil albedo"
    soil_α::FT
end

"""
     function PrescribedSoil{FT}(;
         root_depths::Vector{FT},
         ψ_soil::Function,
         soil_α::FT
     ) where {FT}

An outer constructor for the PrescribedSoil soil driver allowing the user to 
specify the soil parameters by keyword arguments.
"""
function PrescribedSoil(;
    root_depths::Vector{FT} = FT.(-Array(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0),
    ψ_soil::Function = t -> eltype(t)(FT(0.0)),
    soil_α::FT = FT(0.2),
) where {FT}
    return PrescribedSoil{FT}(root_depths, ψ_soil, soil_α)
end

"""
     PrognosticSoil{FT} <: AbstractSoilDriver{FT}

Concrete type of AbstractSoilDriver used for dispatch in cases where both 
a canopy model and soil model are run. Contains the soil albedo to be shared 
between the canopy and soil models.

When running the SoilCanopyModel, the soil model specifies the albedo. In this 
case, the constructor for the model reads the albedo from the soil model and the
user only specifies it once, for the soil model. When running the 
SoilPlantHydrologyModel (which is intended for internal usage/testing primarily)
the user must specify the albedo because the soil model in that case does not
have an albedo.
$(DocStringExtensions.FIELDS)
"""
struct PrognosticSoil{FT} <: AbstractSoilDriver{FT}
    soil_α::FT
end

"""
    function PrognosticSoil{FT}(; soil_α::FT) where {FT}

An outer constructor for the PrognosticSoil soil driver allowing the user to
specify the soil albedo by keyword argument.
"""
function PrognosticSoil(; soil_α::FT = FT(0.2)) where {FT}
    return PrognosticSoil{FT}(soil_α)
end
