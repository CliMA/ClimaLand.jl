export SIFParameters, Lee2015SIFModel

abstract type AbstractSIFModel{FT} <: AbstractCanopyComponent{FT} end

"""
    SIFParameters{FT<:AbstractFloat}

The required parameters for the SIF parameterisation 
Lee et al, 2015. Global Change Biology 21, 3469-3477, doi:10.1111/gcb.12948.
$(DocStringExtensions.FIELDS)
"""
@kwdef struct SIFParameters{FT <: AbstractFloat}
    "The rate coefficient for florescence, unitless"
    kf::FT = FT(0.05)
    "Parameter used to compute the rate coefficient for heat loss in dark-adapted conditions, Tol et al. 2014, unitless"
    kd_p1::FT = FT(0.03)
    "Parameter used to compute the rate coefficient for heat loss in dark-adapted conditions, Tol et al. 2014, unitless"
    kd_p2::FT = FT(0.0273)
    "Parameter used to compute the rate coefficient for heat loss in dark-adapted conditions, Tol et al. 2014, unitless"
    min_kd::FT = FT(0.087)
    "Parameter used to compute the rate coefficient for heat loss in light-adapted conditions, Lee et al 2013 (unitless)"
    kn_p1::FT = FT(6.2473)
    "Parameter used to compute the rate coefficient for heat loss in light-adapted conditions, Lee et al 2013 (unitless)"
    kn_p2::FT = FT(0.5944)
    "Rate coefficient for photochemical quenching"
    kp::FT = FT(4.0)
    "Slope of line relating leaf-level fluorescence to spectrometer-observed fluorescence as a function of Vcmax 25. Lee et al 2015."
    kappa_p1::FT = FT(0.045)
    "Intercept of line relating leaf-level fluorescence to spectrometer-observed fluorescence as a function of Vcmax 25.  Lee et al 2015."
    kappa_p2::FT = FT(7.85)
end

Base.eltype(::SIFParameters{FT}) where {FT} = FT

struct Lee2015SIFModel{FT, SP <: SIFParameters{FT}} <: AbstractSIFModel{FT}
    parameters::SP
    function Lee2015SIFModel{FT}() where {FT}
        parameters = SIFParameters{FT}()
        new{FT, typeof(parameters)}(parameters)
    end
end

ClimaLand.name(model::AbstractSIFModel) = :sif
ClimaLand.auxiliary_vars(model::Lee2015SIFModel) = (:SIF,)
ClimaLand.auxiliary_types(model::Lee2015SIFModel{FT}) where {FT} = (FT,)
ClimaLand.auxiliary_domain_names(::Lee2015SIFModel) = (:surface,)


# 4 Solar Induced Fluorescence (SIF)

# call function below inside photosynthesis.jl p

"""
    update_SIF!(
        SIF::FT,
        sif_model::Lee2015SIFModel,
        APAR::FT,
        Tc::FT,
        Vcmax25::FT,
        R::FT,
        T_freeze::FT,
        photosynthesis_parameters,
        energy_per_mole_photon_PAR,
    ) where {FT}

Updates observed SIF at 755 nm in W/m^2. Note that Tc is in Kelvin, and photo
synthetic rates are in mol/m^2/s, and APAR is in PPFD.
Lee et al, 2015. Global Change Biology 21, 3469-3477, doi:10.1111/gcb.12948
"""
function update_SIF!(
    SIF,
    sif_model::Lee2015SIFModel,
    f_abs_par,
    Tc,
    Vcmax25,
    R,
    T_freeze,
    photosynthesis_parameters,
    energy_per_mole_photon_par,
    par_d,
)
    (; ΔHJmax, To, θj, ϕ) = photosynthesis_parameters
    sif_parameters = sif_model.parameters
    @. SIF = compute_SIF_at_a_point(
        par_d * f_abs_par / energy_per_mole_photon_par,
        Tc,
        Vcmax25,
        R,
        T_freeze,
        ΔHJmax,
        To,
        θj,
        ϕ,
        sif_parameters,
    )
end

Base.broadcastable(m::SIFParameters) = tuple(m)
function compute_SIF_at_a_point(
    APAR::FT,
    Tc::FT,
    Vcmax25::FT,
    R::FT,
    T_freeze::FT,
    ΔHJmax::FT,
    To::FT,
    θj::FT,
    ϕ::FT,
    sif_parameters::SIFParameters{FT},
) where {FT}
    Jmax = max_electron_transport(Vcmax25, ΔHJmax, Tc, To, R)
    J = electron_transport(APAR, Jmax, θj, ϕ)
    (; kf, kd_p1, kd_p2, min_kd, kn_p1, kn_p2, kp, kappa_p1, kappa_p2) =
        sif_parameters
    kd = max(kd_p1 * (Tc - T_freeze) + kd_p2, min_kd)
    x = 1 - J / Jmax
    kn = (kn_p1 * x - kn_p2) * x
    ϕp0 = kp / (kf + kp + kn)
    ϕp = J / Jmax * ϕp0
    ϕf = kf / (kf + kd + kn) * (1 - ϕp)
    κ = kappa_p1 * Vcmax25 * FT(1e6) + kappa_p2 # formula expects Vcmax25 in μmol/m^2/s
    F = APAR * ϕf
    SIF_755 = F / κ

    return SIF_755
end
