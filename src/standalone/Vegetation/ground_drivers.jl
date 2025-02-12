using StaticArrays
export AbstractGroundConditions, PrescribedGroundConditions, ground_albedo

"""
An abstract type of ground conditions for the canopy model;
we will support prescribed or prognostic ground conditions.
"""
abstract type AbstractGroundConditions end

"""
     PrescribedGroundConditions <: AbstractGroundConditions

A container for holding prescribed ground conditions needed by the canopy model
when running the canopy in standalone mode, including the soil pressure, surface
temperature, albedo, and emissivity.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedGroundConditions{
    FT,
    F1 <: Function,
    F2 <: Function,
    SD <: AbstractSpectralDiscretization,
    VEC <: AbstractArray{FT},
} <: AbstractGroundConditions
    "The depth of the root tips, in meters"
    root_depths::VEC
    "Prescribed soil potential (m) in the root zone as a function of time"
    ψ::F1
    "Prescribed ground surface temperature (K) as a function of time"
    T::F2
    "Spectral discretization"
    spectral_discretization::SD
    "Spectral ground albedo"
    α_ground::Tuple
    "Ground emissivity"
    ϵ::FT
end

"""
     function PrescribedGroundConditions(FT;
         root_depths::AbstractArray{FT},
         ψ::Function,
         T::Function,
         α_PAR::FT,
         α_NIR::FT,
         ϵ::FT
     ) where {FT}

An outer constructor for the PrescribedGroundConditions allowing the user to
specify the ground parameters by keyword arguments.
"""
function PrescribedGroundConditions(
    FT;
    root_depths::AbstractArray = (SVector{10, FT}(
        -(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0,
    )),
    ψ::Function = t -> 0.0,
    T::Function = t -> 298.0,
    spectral_discretization::AbstractSpectralDiscretization = ClimaLand.TwoBandSpectralDiscretization{FT}(),
    α_ground::Tuple = FT.((0.2, 0.4)),
    ϵ = FT(0.99),
)
    return PrescribedGroundConditions{
        FT,
        typeof(ψ),
        typeof(T),
        typeof(spectral_discretization),
        typeof(root_depths),
    }(
        root_depths,
        ψ,
        T,
        spectral_discretization,
        α_ground,
        ϵ,
    )
end

"""
    ground_albedo(prognostic_land_components::Val{(:canopy,)}, ground::PrescribedGroundConditions, λ::Float64, _...)

Returns the ground albedo for a PrescribedGroundConditions driver. In this case,
the prognostic_land_components only contain `:canopy`, because the canopy is being run in standalone
mode.
"""
function ground_albedo(
    prognostic_land_components::Val{(:canopy,)},
    ground::PrescribedGroundConditions,
    _...,
)
    return ground.α_ground
end
