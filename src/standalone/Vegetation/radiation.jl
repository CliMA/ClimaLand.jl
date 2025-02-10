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
struct ConstantGFunction{F <: Union{AbstractFloat, ClimaCore.Fields.Field}} <:
       AbstractGFunction{F}
    "Leaf angle distribution value (unitless)"
    ld::F
end

# Make the ConstantGFunction broadcastable
Base.broadcastable(G::ConstantGFunction) = tuple(G)

"""
    CLMGFunction

A type for a G function that is parameterized by the solar zenith angle,
following the CLM approach to parameterizing the leaf angle distribution
function.
"""
struct CLMGFunction{F <: Union{AbstractFloat, ClimaCore.Fields.Field}} <:
       AbstractGFunction{F}
    "Leaf orientation index (unitless)"
    χl::F
end

# Make the CLMGFunction broadcastable
Base.broadcastable(G::CLMGFunction) = tuple(G)

"""
    BeerLambertParameters{FT <: AbstractFloat}

The required parameters for the Beer-Lambert radiative transfer model.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct BeerLambertParameters{
    FT <: AbstractFloat,
    SD <: AbstractSpectralDiscretization{FT},
    G <: AbstractGFunction,
    F <: Union{Tuple, ClimaCore.Fields.Field},
}
    "Discretization of shortwave spectrum"
    spectral_discretization::SD
    "Spectral leaf element reflectances"
    ρ_leaf::F
    "Emissivity of the canopy"
    ϵ_canopy::FT
    "Clumping index following Braghiere (2021) (unitless)"
    Ω::FT
    "Typical wavelength per PAR photon (m)"
    λ_γ_PAR::FT
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
    SD <: AbstractSpectralDiscretization{FT},
    G <: AbstractGFunction,
    F <: Union{FT, ClimaCore.Fields.Field},
    TF <: Union{Tuple, ClimaCore.Fields.Field},
}
    "Discretization of shortwave spectrum"
    spectral_discretization::SD
    "Spectral leaf element reflectances"
    ρ_leaf::TF
    "Spectral leaf element transmittances"
    τ_leaf::TF
    "Emissivity of the canopy"
    ϵ_canopy::FT
    "Clumping index following Braghiere 2021 (unitless)"
    Ω::F
    "Typical wavelength per PAR photon (m)"
    λ_γ_PAR::FT
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
    compute_spectral_sw!(sw_in,
        model::AbstractRadiationModel,
        solar_radiation::ClimaLand.PrescribedRadiativeFluxes,
        p,
        t,
    )

Update `sw_in` with the spectral shortwave radiation (W/m^2) given the input
solar radiation for a radiative transfer model.

The spectral shortwave radiation per each band is computed based on the relative
proportion of the total shortwave radiation from the solar irradiance curve.
"""
function compute_spectral_sw!(
    sw_in,
    model::AbstractRadiationModel,
    solar_radiation::ClimaLand.PrescribedRadiativeFluxes,
    p,
    t,
)
    @. sw_in = p.drivers.SW_d .* model.parameters.spectral_discretization.I
end

# Make radiation models broadcastable
Base.broadcastable(RT::AbstractRadiationModel) = tuple(RT)

ClimaLand.name(model::AbstractRadiationModel) = :radiative_transfer
ClimaLand.auxiliary_vars(model::Union{BeerLambertModel, TwoStreamModel}) =
    (:SW_d, :rt, :LW_n, :SW_n, :ϵ, :frac_diff, :G, :K)
ClimaLand.auxiliary_types(
    model::Union{BeerLambertModel{FT}, TwoStreamModel{FT}},
) where {FT} = (
    FT,
    NTuple{
        length(model.parameters.spectral_discretization.λ - 1),
        NamedTuple{(:abs, :refl, :trans), Tuple{FT, FT, FT}},
    },
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
)
ClimaLand.auxiliary_domain_names(::Union{BeerLambertModel, TwoStreamModel}) = (
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
                                  ground::PrescribedGroundConditions
                                  canopy,
                                  radiation::PrescribedRadiativeFluxes,
                                  earth_param_set::PSE,
                                  Y::ClimaCore.Fields.FieldVector,
                                  t,
                                 ) where {PSE}


Computes and stores the net long and short wave radiation, in W/m^2, over all
bands, absorbed by the canopy when the canopy is run in standalone mode, with
only a :canopy model as a prognostic component, with PrescribedGroundConditions.

LW and SW net radiation are stored in `p.canopy.radiative_transfer.LW_n`
and `p.canopy.radiative_transfer.SW_n`.
"""
function canopy_radiant_energy_fluxes!(
    p::NamedTuple,
    ground::PrescribedGroundConditions,
    canopy,
    radiation::PrescribedRadiativeFluxes,
    earth_param_set::PSE,
    Y::ClimaCore.Fields.FieldVector,
    t,
) where {PSE}
    FT = eltype(earth_param_set)
    SW_d = p.canopy.radiative_transfer.SW_d
    SW_abs = ntuple(
        i -> p.canopy.radiative_transfer.rt[i].abs,
        length(p.canopy.radiative_transfer.rt),
    )
    @. p.canopy.radiative_transfer.SW_n = sum(SW_abs .* SW_d)
    ϵ_canopy = p.canopy.radiative_transfer.ϵ # this takes into account LAI/SAI
    # Long wave: use ground conditions from the ground driver
    T_ground::FT = ground.T(t)
    ϵ_ground = ground.ϵ
    _σ = FT(LP.Stefan(earth_param_set))
    LW_d = p.drivers.LW_d
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p, t)
    LW_d_canopy = @. (1 - ϵ_canopy) * LW_d + ϵ_canopy * _σ * T_canopy^4
    LW_u_ground = @. ϵ_ground * _σ * T_ground^4 + (1 - ϵ_ground) * LW_d_canopy
    @. p.canopy.radiative_transfer.LW_n =
        ϵ_canopy * LW_d - 2 * ϵ_canopy * _σ * T_canopy^4 +
        ϵ_canopy * LW_u_ground
end
