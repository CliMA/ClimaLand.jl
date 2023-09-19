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
    "PAR leaf reflectance (unitless)"
    α_PAR_leaf::FT
    "NIR leaf reflectance"
    α_NIR_leaf::FT
    "Clumping index following Braghiere (2021) (unitless)"
    Ω::FT
    "Typical wavelength per PAR photon (m)"
    λ_γ_PAR::FT
    "Typical wavelength per NIR photon (m)"
    λ_γ_NIR::FT
end

"""
    function BeerLambertParameters{FT}(;
        ld = FT(0.5),    
        α_PAR_leaf = FT(0.1),
        α_NIR_leaf = FT(0.4),
        Ω = FT(1),
        λ_γ_PAR = FT(5e-7),
        λ_γ_NIR = FT(1.65e-6),
    ) where {FT}

A constructor supplying default values for the BeerLambertParameters struct.
    """
function BeerLambertParameters{FT}(;
    ld = FT(0.5),
    α_PAR_leaf = FT(0.1),
    α_NIR_leaf = FT(0.4),
    Ω = FT(1),
    λ_γ_PAR = FT(5e-7),
    λ_γ_NIR = FT(1.65e-6),
) where {FT}
    return BeerLambertParameters{FT}(
        ld,
        α_PAR_leaf,
        α_NIR_leaf,
        Ω,
        λ_γ_PAR,
        λ_γ_NIR,
    )
end

struct BeerLambertModel{FT} <: AbstractRadiationModel{FT}
    parameters::BeerLambertParameters{FT}
end


"""
    TwoStreamParameters{FT <: AbstractFloat}

The required parameters for the two-stream radiative transfer model.
$(DocStringExtensions.FIELDS)
"""
struct TwoStreamParameters{FT <: AbstractFloat}
    "Leaf angle distribution function (unitless)"
    ld::FT
    "PAR leaf reflectance (unitless)"
    α_PAR_leaf::FT
    "PAR leaf element transmittance"
    τ_PAR_leaf::FT
    "NIR leaf reflectance"
    α_NIR_leaf::FT
    "NIR leaf element transmittance"
    τ_NIR_leaf::FT
    "Clumping index following Braghiere (2021) (unitless)"
    Ω::FT
    "Typical wavelength per PAR photon (m)"
    λ_γ_PAR::FT
    "Typical wavelength per NIR photon (m)"
    λ_γ_NIR::FT
    "number of layers to simulate radiative transfer through"
    n_layers::UInt64
end

"""
    function TwoStreamParameters{FT}(;
        ld = FT(0.5),
        α_PAR_leaf = FT(0.3),
        τ_PAR_leaf = FT(0.2),
        α_NIR_leaf = FT(0.4),
        τ_NIR_leaf = FT(0.25),
        Ω = FT(1),
        λ_γ_PAR = FT(5e-7),
        λ_γ_NIR = FT(1.65e-6),
        n_layers = UInt64(20),
        diff_perc = FT(0)
    ) where {FT}

A constructor supplying default values for the TwoStreamParameters struct.
"""
function TwoStreamParameters{FT}(;
    ld = FT(0.5),
    α_PAR_leaf = FT(0.3),
    τ_PAR_leaf = FT(0.2),
    α_NIR_leaf = FT(0.4),
    τ_NIR_leaf = FT(0.25),
    Ω = FT(1),
    λ_γ_PAR = FT(5e-7),
    λ_γ_NIR = FT(1.65e-6),
    n_layers = UInt64(20),
) where {FT}
    return TwoStreamParameters{FT}(
        ld,
        α_PAR_leaf,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        Ω,
        λ_γ_PAR,
        λ_γ_NIR,
        n_layers,
    )
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

"""
    compute_NIR(
        model::AbstractRadiationModel,
        solar_radiation::ClimaLSM.PrescribedRadiativeFluxes,
        t,
    )

Returns the estimated NIR (W/m^2) given the input solar radiation
for a radiative transfer model.

The estimated PNIR is half of the incident shortwave radiation.
"""
function compute_NIR(
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
    (:apar, :par, :anir, :nir)
ClimaLSM.auxiliary_types(
    model::Union{BeerLambertModel{FT}, TwoStreamModel{FT}},
) where {FT} = (FT, FT, FT, FT)
ClimaLSM.auxiliary_domain_names(::Union{BeerLambertModel, TwoStreamModel}) =
    (:surface, :surface, :surface, :surface)
