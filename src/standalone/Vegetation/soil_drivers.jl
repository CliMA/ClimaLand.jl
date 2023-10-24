export AbstractSoilDriver, PrescribedSoil

"""
An abstract type of soil drivers of the canopy model.
"""
abstract type AbstractSoilDriver{FT <: AbstractFloat} end

"""
     PrescribedSoil{FT} <: AbstractSoilDriver{FT}

A container for holding prescribed soil parameters needed by the canopy model
when running the canopy in standalone mode.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedSoil{FT} <: AbstractSoilDriver{FT}
    "The depth of the root tips, in meters"
    root_depths::Vector{FT}
    "Prescribed soil potential (m) as a function of time"
    ψ_soil::Function
    "Prescribed soil temperature (K) as a function of time"
    T_soil::Function
    "Soil albedo for PAR"
    soil_α_PAR::FT
    "Soil albedo for NIR"
    soil_α_NIR::FT
    "Soil emissivity"
    ϵ_soil::FT
end

"""
     function PrescribedSoil{FT}(;
                                 root_depths::Vector{FT} = FT.(-Array(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0),
                                 ψ_soil::Function = t -> eltype(t)(FT(0.0)),
                                 T_soil::Function = t -> eltype(t)(FT(0.0)),
                                 soil_α_PAR::FT = FT(0.2),
                                 soil_α_NIR::FT = FT(0.4),
                                 ϵ_soil::FT = FT(0.98),
                                 ) where {FT}

An outer constructor for the PrescribedSoil soil driver allowing the user to 
specify the soil parameters by keyword arguments.
"""
function PrescribedSoil(;
    root_depths::Vector{FT} = FT.(-Array(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0),
    ψ_soil::Function = t -> eltype(t)(FT(0.0)),
    T_soil::Function = t -> eltype(t)(FT(0.0)),
    soil_α_PAR::FT = FT(0.2),
    soil_α_NIR::FT = FT(0.4),
    ϵ_soil::FT = FT(0.98),
) where {FT}
    return PrescribedSoil{FT}(
        root_depths,
        ψ_soil,
        T_soil,
        soil_α_PAR,
        soil_α_NIR,
        ϵ_soil,
    )
end
