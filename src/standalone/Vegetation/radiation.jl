export SpectralDiscretization,
    BeerLambertParameters,
    BeerLambertModel,
    TwoStreamParameters,
    TwoStreamModel,
    canopy_radiant_energy_fluxes!,
    ConstantGFunction,
    CLMGFunction

abstract type AbstractRadiationModel{FT} <: AbstractCanopyComponent{FT} end

abstract type AbstractGFunction{FT} end

"""
    SpectralDiscretization

A type used to represent the spectral discretization of radiation models.
Consists of the wavelength boundaries between each band in a discretization.
"""
struct SpectralDiscretization{FT <: AbstractFloat}
    "Wavelength boundaries between each band in the discretization (m)"
    λ::Tuple{Vararg{FT}}
    "Banded irradiance curve for the solar radiation (W/m^2)"
    I::Tuple{Vararg{FT}}
    "Proportions of PAR radiation in each band"
    PAR_proportions::Tuple{Vararg{FT}}
    "Proportions of NIR radiation in each band"
    NIR_proportions::Tuple{Vararg{FT}}
end

"""
    SolarIrradianceCurve(λ::FT) where {FT}

Function to approximate the solar irradiance curve at a given wavelength. The
solar irradiance curve gices the intensity of solar radiation over the spectrum.
"""
function SolarIrradianceCurve(λ::FT) where {FT}
    return (2 * 6.26e-34 * 3e8 / λ^5) * (1 / (exp(6.26e-34 * 3e8 / (λ * 1.381e-23 * 5778)) - 1))
end

"""
    ComputeBandedIrradianceCurve(
        λ::Tuple{Vararg{FT}},
    ) where {FT}

Assigns a proportion of total incoming solar radiation to each band in the
spectral discretization to allow the computation of the input radiation for
a discretization with an arbitrary number of bands.
"""
function ComputeBandedIrradianceCurve(λ::Tuple{Vararg{FT}}) where {FT}
    # Find the values of the solar irradiance curve at band separating points
    curve_points = solarIrradianceCurve.(λ)
    # Use the trapezoid method to compute the area under the curve between each
    trap_sums = MVector{length(λ) - 1, FT}(zeros(length(λ) - 1)...)
    for i in 1:length(λ) - 1
        trap_sums[i] = (curve_points[i] + curve_points[i + 1]) * (λ[i + 1] - λ[i]) / 2
    end
    # Now take each band's area as a proportion of the total area under the
    # curve
    return Tuple(trap_sums ./ sum(trap_sums))
end

"""
    ComputeBandedProportions(
        λ::Tuple{Vararg{FT}},
    ) where {FT}

Use the solar irradiance curve to compute the proportion of each band in a
spectral discretization that is within a specific range, used for PAR and NIR
radiation.
"""
function ComputeBandedProportions(λ::Tuple{Vararg{FT}}, bounds::Tuple{FT, FT}) where {FT}
    # PAR proportions for each band
    proportions = MVector{length(λ) - 1, FT}(zeros(length(λ) - 1)...)

    # In each band, determine whether the band is entirely within the PAR range,
    # partly in the PAR range, or entirely outside the PAR range
    for i in 1:length(λ) - 1
        # If the band is entirely within the PAR range, the proportion of PAR
        # radiation in the band is 1
        if λ[i] >= bounds[1] && λ[i + 1] <= bounds[2]
            proportions[i] = 1
        # If the band is entirely outside the PAR range, the proportion of PAR
        # radiation in the band is 0
        elseif λ[i] >= bounds[2] || λ[i + 1] <= bounds[1]
            proportions[i] = 0
        # If the band is partly in the PAR range, we need to know what part of 
        # the band is in the PAR range
        elseif λ[i] <= bounds[1]
            # trapezoidal area of the band that is in the PAR range
            bound_trap = (solarIrradianceCurve(λ[i]) + solarIrradianceCurve(bounds[2])) * (bounds[2] - λ[i]) / 2
            # trapezoidal area of the band
            trap = (solarIrradianceCurve(λ[i]) + solarIrradianceCurve(λ[i + 1])) * (λ[i + 1] - λ[i]) / 2
            # proportion of PAR radiation in the band
            proportions[i] = bound_trap / trap
        else
            # trapezoidal area of the band that is in the PAR range
            bound_trap = (solarIrradianceCurve(λ[i]) + solarIrradianceCurve(bounds[2])) * (bounds[2] - λ[i]) / 2
            # trapezoidal area of the band
            trap = (solarIrradianceCurve(λ[i]) + solarIrradianceCurve(λ[i + 1])) * (λ[i + 1] - λ[i]) / 2
            # proportion of PAR radiation in the band
            PAR_proportions[i] = bound_trap / trap
        end
    end
    return Tuple(proportions)
end

"""
    ComputeBandedPAR(
        λ::Tuple{Vararg{FT}},
    ) where {FT}

Use the solar irradiance curve to compute the proportion of each band in a
spectral discretization that is within the PAR range.
"""
function ComputeBandedPAR(λ::Tuple{Vararg{FT}}) where {FT}
    return ComputeBandedProportions(λ, (400e-9, 700e-9))
end

"""
    ComputeBandedNIR(
        λ::Tuple{Vararg{FT}},
    ) where {FT}
Use the solar irradiance curve to compute the proportion of each band in a
spectral discretization that is within the NIR range.
"""
function ComputeBandedNIR(λ::Tuple{Vararg{FT}}) where {FT}
    return ComputeBandedProportions(λ, (700e-9, 2500e-9))
end

"""
    SpectralDiscretization(
        λ::Tuple{Vararg{FT}},
    ) where {FT}
    
Creates a spectral discretization with the given wavelength boundaries.
"""
function SpectralDiscretization(λ::Tuple{Vararg{FT}}) where {FT}
    # Check that the given boundaries are in ascending order and cover the full
    # SW range from 100 nm to 3000 nm
    @assert λ[1] == 100e-9 && λ[end] == 3000e-9
    @assert all(λ[i] < λ[i + 1] for i in 1:length(λ) - 1)
    I = ComputeBandedIrradianceCurve(λ)
    PAR_proportions = ComputeBandedPAR(λ)
    NIR_proportions = ComputeBandedNIR(λ)
    return SpectralDiscretization{FT}(λ, I, PAR_proportions, NIR_proportions)
end

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
    SD <: SpectralDiscretization{FT},
    G <: AbstractGFunction,
    F <: Union{FT, ClimaCore.Fields.Field},
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
    SD <: SpectralDiscretization{FT},
    G <: AbstractGFunction,
    F <: Union{FT, ClimaCore.Fields.TF <: Union{Tuple, ClimaCore.Fields.Field}},
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
    NTuple{length(model.parameters.spectral_discretization.λ - 1), NamedTuple{(:abs, :refl, :trans), Tuple{FT, FT, FT}}},
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
    SW_abs = ntuple(i->p.canopy.radiative_transfer.rt[i].abs, length(p.canopy.radiative_transfer.rt))
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
