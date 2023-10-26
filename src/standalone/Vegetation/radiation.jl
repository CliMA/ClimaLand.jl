using ..ClimaLSM.Canopy: AbstractSoilDriver

export BeerLambertParameters,
    BeerLambertModel,
    TwoStreamParameters,
    TwoStreamModel,
    canopy_radiant_energy_fluxes!

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
    "Emissivity of the canopy"
    ϵ_canopy::FT
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
        ϵ_canopy = FT(0.98),
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
    ϵ_canopy = FT(0.98),
    Ω = FT(1),
    λ_γ_PAR = FT(5e-7),
    λ_γ_NIR = FT(1.65e-6),
) where {FT}
    return BeerLambertParameters{FT}(
        ld,
        α_PAR_leaf,
        α_NIR_leaf,
        ϵ_canopy,
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
    "Emissivity of the canopy"
    ϵ_canopy::FT
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
        ϵ_canopy = FT(0.98),
        Ω = FT(1),
        λ_γ_PAR = FT(5e-7),
        λ_γ_NIR = FT(1.65e-6),
        n_layers = UInt64(20),
    ) where {FT}

A constructor supplying default values for the TwoStreamParameters struct.
    """
function TwoStreamParameters{FT}(;
    ld = FT(0.5),
    α_PAR_leaf = FT(0.3),
    τ_PAR_leaf = FT(0.2),
    α_NIR_leaf = FT(0.4),
    τ_NIR_leaf = FT(0.25),
    ϵ_canopy = FT(0.98),
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
        ϵ_canopy,
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
    (:apar, :par, :rpar, :tpar, :anir, :nir, :rnir, :tnir, :LW_n, :SW_n)
ClimaLSM.auxiliary_types(
    model::Union{BeerLambertModel{FT}, TwoStreamModel{FT}},
) where {FT} = (FT, FT, FT, FT, FT, FT, FT, FT, FT, FT)
ClimaLSM.auxiliary_domain_names(::Union{BeerLambertModel, TwoStreamModel}) = (
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
                                  s::PrescribedSoil{FT},
                                  canopy,
                                  radiation::PrescribedRadiativeFluxes,
                                  earth_param_set::PSE,
                                  Y::ClimaCore.Fields.FieldVector,
                                  t,
                                 ) where {FT, PSE}


Computes and stores the net long and short wave radition, in W/m^2,
absorbed by the canopy when the canopy is run in standalone mode,
with a PrescribedSoil conditions.

LW and SW net radiation are stored in `p.canopy.radiative_transfer.LW_n`
and `p.canopy.radiative_transfer.SW_n`.
"""
function canopy_radiant_energy_fluxes!(
    p::NamedTuple,
    s::PrescribedSoil{FT},
    canopy,
    radiation::PrescribedRadiativeFluxes,
    earth_param_set::PSE,
    Y::ClimaCore.Fields.FieldVector,
    t,
) where {FT, PSE}

    # Short wave makes use of precomputed APAR and ANIR
    # in moles of photons per m^2 per s
    c = FT(LSMP.light_speed(earth_param_set))
    h = FT(LSMP.planck_constant(earth_param_set))
    N_a = FT(LSMP.avogadro_constant(earth_param_set))
    (; α_PAR_leaf, λ_γ_PAR, λ_γ_NIR, ϵ_canopy) =
        canopy.radiative_transfer.parameters
    APAR = p.canopy.radiative_transfer.apar
    ANIR = p.canopy.radiative_transfer.anir
    energy_per_photon_PAR = h * c / λ_γ_PAR
    energy_per_photon_NIR = h * c / λ_γ_NIR
    @. p.canopy.radiative_transfer.SW_n =
        (energy_per_photon_PAR * N_a * APAR) +
        (energy_per_photon_NIR * N_a * ANIR)

    # Long wave: use soil conditions from the PrescribedSoil driver
    T_soil = s.T(t)
    ϵ_soil = s.ϵ
    _σ = FT(LSMP.Stefan(earth_param_set))
    LW_d::FT = canopy.radiation.LW_d(t)

    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p, t)
    LW_d_canopy = @. (1 - ϵ_canopy) * LW_d + ϵ_canopy * _σ * T_canopy^4
    LW_u_soil = @. ϵ_soil * _σ * T_soil^4 + (1 - ϵ_soil) * LW_d_canopy
    @. p.canopy.radiative_transfer.LW_n =
        ϵ_canopy * LW_d - 2 * ϵ_canopy * _σ * T_canopy^4 + ϵ_canopy * LW_u_soil
end
