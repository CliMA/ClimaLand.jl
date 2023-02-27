export BeerLambertParameters, BeerLambertModel

abstract type AbstractRadiationModel{FT} <: AbstractCanopyComponent{FT} end

"""
    BeerLambertParameters{FT <: AbstractFloat}

The required parameters for the Beer-Lambert radiative transfer model.
$(DocStringExtensions.FIELDS)
"""
struct BeerLambertParameters{FT <: AbstractFloat}
    "Leaf angle distribution function (unitless)"
    ld::FT
    "PAR canopy reflectance (unitless)"
    ρ_leaf::FT
    "Clumping index following Braghiere (2021) (unitless)"
    Ω::FT
    "Typical wavelength per PAR photon (m)"
    λ_γ::FT
end

"""
    function BeerLambertParameters{FT}(;
        ld = FT(0.5),    
        ρ_leaf = FT(0.1),
        Ω = FT(1),
        λ_γ = FT(5e-7)
    ) where {FT}

A constructor supplying default values for the BeerLambertParameters struct.
    """
function BeerLambertParameters{FT}(;
    ld = FT(0.5),
    ρ_leaf = FT(0.1),
    Ω = FT(1),
    λ_γ = FT(5e-7),
) where {FT}
    return BeerLambertParameters{FT}(ld, ρ_leaf, Ω, λ_γ)
end

struct BeerLambertModel{FT} <: AbstractRadiationModel{FT}
    parameters::BeerLambertParameters{FT}
end

ClimaLSM.name(model::AbstractRadiationModel) = :radiative_transfer
ClimaLSM.auxiliary_vars(model::BeerLambertModel) = (:apar,)
ClimaLSM.auxiliary_types(model::BeerLambertModel{FT}) where {FT} = (FT,)
