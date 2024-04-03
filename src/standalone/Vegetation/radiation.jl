using ..ClimaLand.Canopy: AbstractSoilDriver

export BeerLambertParameters,
    BeerLambertModel,
    TwoStreamParameters,
    TwoStreamModel,
    canopy_radiant_energy_fluxes!,
    ConstantGFunction,
    CLMGFunction

abstract type AbstractRadiationModel{FT} <: AbstractCanopyComponent{FT} end

abstract type AbstractGFunction{FT} end

"""
    ConstantGFunction

A type for a constant G function, which is used to represent the leaf angle
distribution function in the radiative transfer models.
"""
struct ConstantGFunction{FT} <: AbstractGFunction{FT}
    "Leaf angle distribution value (unitless)"
    ld::FT
end

# Make the ConstantGFunction broadcastable
Base.broadcastable(G::ConstantGFunction) = Ref(G)

"""
    CLMGFunction

A type for a G function that is parameterized by the solar zenith angle,
following the CLM approach to parameterizing the leaf angle distribution function.
"""
struct CLMGFunction{FT} <: AbstractGFunction{FT}
    "Leaf orientation index (unitless)"
    χl::FT
end

# Make the CLMGFunction broadcastable
Base.broadcastable(G::CLMGFunction) = Ref(G)

"""
    BeerLambertParameters{FT <: AbstractFloat}

The required parameters for the Beer-Lambert radiative transfer model.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct BeerLambertParameters{
    FT <: AbstractFloat,
    G <: AbstractGFunction{FT},
}
    "PAR leaf reflectance (unitless)"
    α_PAR_leaf::FT
    "NIR leaf reflectance"
    α_NIR_leaf::FT
    "Emissivity of the canopy"
    ϵ_canopy::FT
    "Clumping index following Braghiere (2021) (unitless)"
    Ω::FT
    "Typical wavelength per PAR photon (m)"
    λ_γ_PAR::FT
    "Typical wavelength per NIR photon (m)"
    λ_γ_NIR::FT
    "Leaf angle distribution function"
    G_Function::G
end

Base.eltype(::BeerLambertParameters{FT}) where {FT} = FT

struct BeerLambertModel{FT, BLP <: BeerLambertParameters{FT}} <:
       AbstractRadiationModel{FT}
    parameters::BLP
end

function BeerLambertModel{FT}(
    parameters::BeerLambertParameters{FT},
) where {FT <: AbstractFloat}
    return BeerLambertModel{eltype(parameters), typeof(parameters)}(parameters)
end

"""
    TwoStreamParameters{FT <: AbstractFloat}

The required parameters for the two-stream radiative transfer model.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct TwoStreamParameters{
    FT <: AbstractFloat,
    G <: AbstractGFunction{FT},
}
    "PAR leaf reflectance (unitless)"
    α_PAR_leaf::FT
    "PAR leaf element transmittance"
    τ_PAR_leaf::FT
    "NIR leaf reflectance"
    α_NIR_leaf::FT
    "NIR leaf element transmittance"
    τ_NIR_leaf::FT
    "Emissivity of the canopy"
    ϵ_canopy::FT
    "Clumping index following Braghiere 2021 (unitless)"
    Ω::FT
    "Typical wavelength per PAR photon (m)"
    λ_γ_PAR::FT
    "Typical wavelength per NIR photon (m)"
    λ_γ_NIR::FT
    "Number of layers to partition the canopy into when integrating the
    absorption over the canopy vertically. Unrelated to the number of layers in 
    the vertical discretization of the canopy for the plant hydraulics model.
    (Constant, and should eventually move to ClimaParams)"
    n_layers::UInt64
    "Leaf angle distribution function"
    G_Function::G
end

Base.eltype(::TwoStreamParameters{FT}) where {FT} = FT

struct TwoStreamModel{FT, TSP <: TwoStreamParameters{FT}} <:
       AbstractRadiationModel{FT}
    parameters::TSP
end

function TwoStreamModel{FT}(
    parameters::TwoStreamParameters{FT},
) where {FT <: AbstractFloat}
    return TwoStreamModel{eltype(parameters), typeof(parameters)}(parameters)
end

"""
    compute_PAR(
        model::AbstractRadiationModel,
        solar_radiation::ClimaLand.PrescribedRadiativeFluxes,
        p,
        t,
    )

Returns the estimated PAR (W/m^2) given the input solar radiation
for a radiative transfer model.

The estimated PAR is half of the incident shortwave radiation.
"""
function compute_PAR(
    model::AbstractRadiationModel,
    solar_radiation::ClimaLand.PrescribedRadiativeFluxes,
    p,
    t,
)
    return p.drivers.SW_d ./ 2
end

"""
    compute_NIR(
        model::AbstractRadiationModel,
        solar_radiation::ClimaLand.PrescribedRadiativeFluxes,
        p,
        t,
    )

Returns the estimated NIR (W/m^2) given the input solar radiation
for a radiative transfer model.

The estimated PNIR is half of the incident shortwave radiation.
"""
function compute_NIR(
    model::AbstractRadiationModel,
    solar_radiation::ClimaLand.PrescribedRadiativeFluxes,
    p,
    t,
)
    return p.drivers.SW_d ./ 2
end

# Make radiation models broadcastable
Base.broadcastable(RT::AbstractRadiationModel) = tuple(RT)

ClimaLand.name(model::AbstractRadiationModel) = :radiative_transfer
ClimaLand.auxiliary_vars(model::Union{BeerLambertModel, TwoStreamModel}) =
    (:apar, :par, :rpar, :tpar, :anir, :nir, :rnir, :tnir, :LW_n, :SW_n, :ϵ)
ClimaLand.auxiliary_types(
    model::Union{BeerLambertModel{FT}, TwoStreamModel{FT}},
) where {FT} = (FT, FT, FT, FT, FT, FT, FT, FT, FT, FT, FT)
ClimaLand.auxiliary_domain_names(::Union{BeerLambertModel, TwoStreamModel}) = (
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
)

"""
    canopy_radiant_energy_fluxes!(p::NamedTuple,
                                  s::PrescribedSoil,
                                  canopy,
                                  radiation::PrescribedRadiativeFluxes,
                                  earth_param_set::PSE,
                                  Y::ClimaCore.Fields.FieldVector,
                                  t,
                                 ) where {PSE}


Computes and stores the net long and short wave radition, in W/m^2,
absorbed by the canopy when the canopy is run in standalone mode,
with a PrescribedSoil conditions.

LW and SW net radiation are stored in `p.canopy.radiative_transfer.LW_n`
and `p.canopy.radiative_transfer.SW_n`.
"""
function canopy_radiant_energy_fluxes!(
    p::NamedTuple,
    s::PrescribedSoil,
    canopy,
    radiation::PrescribedRadiativeFluxes,
    earth_param_set::PSE,
    Y::ClimaCore.Fields.FieldVector,
    t,
) where {PSE}
    FT = eltype(earth_param_set)

    # Short wave makes use of precomputed APAR and ANIR
    # in moles of photons per m^2 per s
    c = FT(LP.light_speed(earth_param_set))
    h = FT(LP.planck_constant(earth_param_set))
    N_a = FT(LP.avogadro_constant(earth_param_set))
    (; α_PAR_leaf, λ_γ_PAR, λ_γ_NIR) = canopy.radiative_transfer.parameters
    APAR = p.canopy.radiative_transfer.apar
    ANIR = p.canopy.radiative_transfer.anir
    energy_per_photon_PAR = h * c / λ_γ_PAR
    energy_per_photon_NIR = h * c / λ_γ_NIR
    @. p.canopy.radiative_transfer.SW_n =
        (energy_per_photon_PAR * N_a * APAR) +
        (energy_per_photon_NIR * N_a * ANIR)
    ϵ_canopy = p.canopy.radiative_transfer.ϵ # this takes into account LAI/SAI
    # Long wave: use soil conditions from the PrescribedSoil driver
    T_soil::FT = s.T(t)
    ϵ_soil = s.ϵ
    _σ = FT(LP.Stefan(earth_param_set))
    LW_d = p.drivers.LW_d
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p, t)
    LW_d_canopy = @. (1 - ϵ_canopy) * LW_d + ϵ_canopy * _σ * T_canopy^4
    LW_u_soil = @. ϵ_soil * _σ * T_soil^4 + (1 - ϵ_soil) * LW_d_canopy
    @. p.canopy.radiative_transfer.LW_n =
        ϵ_canopy * LW_d - 2 * ϵ_canopy * _σ * T_canopy^4 + ϵ_canopy * LW_u_soil
end
