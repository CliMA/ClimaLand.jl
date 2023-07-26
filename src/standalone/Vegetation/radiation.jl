using ..ClimaLSM.Canopy: AbstractSoilDriver

export BeerLambertParameters,
    BeerLambertModel, TwoStreamParameters, TwoStreamModel

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

ClimaLSM.auxiliary_vars(model::BeerLambertModel) = (:apar, :par)
ClimaLSM.auxiliary_types(model::BeerLambertModel{FT}) where {FT} = (FT, FT)

"""
    TwoStreamParameters{FT <: AbstractFloat}

The required parameters for the two-stream radiative transfer model.
$(DocStringExtensions.FIELDS)
"""
struct TwoStreamParameters{FT <: AbstractFloat}
    "Leaf angle distribution function (unitless)"
    ld::FT
    "PAR canopy reflectance (unitless)"
    ρ_leaf::FT
    "leaf element transmittance"
    τ::FT
    "Clumping index following Braghiere (2021) (unitless)"
    Ω::FT
    "Typical wavelength per PAR photon (m)"
    λ_γ::FT
    "number of layers to simulate radiative transfer through"
    n_layers::UInt64
    "Proportion of diffuse radiation"
    diff_perc::FT
end

"""
    function TwoStreamParameters{FT}(;
        ld = FT(0.5),
        ρ_leaf = FT(0.3),
        τ = FT(0.2),
        Ω = FT(1),
        λ_γ = FT(5e-7),
        n_layers = UInt64(20),
        diff_perc = FT(0)
    ) where {FT}

A constructor supplying default values for the TwoStreamParameters struct.
"""
function TwoStreamParameters{FT}(;
    ld = FT(0.5),
    ρ_leaf = FT(0.3),
    τ = FT(0.2),
    Ω = FT(1),
    λ_γ = FT(5e-7),
    n_layers = UInt64(20),
    diff_perc = FT(0),
) where {FT}
    return TwoStreamParameters{FT}(ld, ρ_leaf, τ, Ω, λ_γ, n_layers, diff_perc)
end

struct TwoStreamModel{FT} <: AbstractRadiationModel{FT}
    parameters::TwoStreamParameters{FT}
end

"""
    compute_PAR(
        model::AbstractRadiationModel,
        solar_radiation::ClimaLSM.PrescribedRadiativeFluxes,
        t,
    )

Returns the estimated PAR (W/m^2) given the input solar radiation
for a radiative transfer model.

The estimated PAR is half of the incident shortwave radiation.
"""
function compute_PAR(
    model::AbstractRadiationModel,
    solar_radiation::ClimaLSM.PrescribedRadiativeFluxes,
    t,
)
    return solar_radiation.SW_d(t) / 2
end

# Make radiation models broadcastable
Base.broadcastable(RT::AbstractRadiationModel) = tuple(RT)

ClimaLSM.name(model::AbstractRadiationModel) = :radiative_transfer
ClimaLSM.auxiliary_vars(model::Union{BeerLambertModel, TwoStreamModel}) =
    (:apar, :par)
ClimaLSM.auxiliary_types(
    model::Union{BeerLambertModel{FT}, TwoStreamModel{FT}},
) where {FT} = (FT, FT)
