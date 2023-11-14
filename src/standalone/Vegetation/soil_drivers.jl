export AbstractSoilDriver, PrescribedSoil, ground_albedo_NIR, ground_albedo_PAR

"""
An abstract type of soil drivers of the canopy model.
"""
abstract type AbstractSoilDriver{FT <: AbstractFloat} end

"""
     PrescribedSoil{FT} <: AbstractSoilDriver{FT}

A container for holding prescribed soil parameters needed by the canopy model
when running the canopy in standalone mode, including the soil pressure, surface
temperature, and albedo.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedSoil{FT} <: AbstractSoilDriver{FT}
    "The depth of the root tips, in meters"
    root_depths::Vector{FT}
    "Prescribed soil potential (m) in the root zone  a function of time"
    ψ::Function
    "Prescribed soil surface temperature (K) as a function of time"
    T::Function
    "Soil albedo for PAR"
    α_PAR::FT
    "Soil albedo for NIR"
    α_NIR::FT
    "Soil emissivity"
    ϵ::FT
end

"""
     function PrescribedSoil{FT}(;
         root_depths::Vector{FT},
         ψ::Function,
         T::Function,
         α_PAR::FT,
         α_NIR::FT,
         ϵ::FT
     ) where {FT}

An outer constructor for the PrescribedSoil soil driver allowing the user to
specify the soil parameters by keyword arguments.
"""
function PrescribedSoil{FT}(;
    root_depths::Vector{FT} = FT.(-Array(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0),
    ψ::Function = t -> 0.0,
    T::Function = t -> 298.0,
    α_PAR::FT = FT(0.2),
    α_NIR::FT = FT(0.4),
    ϵ::FT = FT(0.99),
) where {FT}
    return PrescribedSoil{FT}(root_depths, ψ, T, α_PAR, α_NIR, ϵ)
end


"""
    ground_albedo_PAR(soil_driver::PrescribedSoil, _...)

Returns the soil albedo in the PAR for a PrescribedSoil driver.
"""
function ground_albedo_PAR(soil_driver::PrescribedSoil, _...)
    return soil_driver.α_PAR
end

"""
    ground_albedo_NIR(soil_driver::PrescribedSoil, _...)

Returns the soil albedo in the NIR for a PrescribedSoil driver.
"""
function ground_albedo_NIR(soil_driver::PrescribedSoil, _...)
    return soil_driver.α_NIR
end
