export SIFParameters, FarquharParameters, FarquharModel, C3, C4

abstract type AbstractPhotosynthesisModel{FT} <: AbstractCanopyComponent{FT} end

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

abstract type AbstractPhotosynthesisMechanism end
"""
    C3 <: AbstractPhotosynthesisMechanism

Helper struct for dispatching between C3 and C4 photosynthesis.
"""
struct C3 <: AbstractPhotosynthesisMechanism end

"""
    C4 <: AbstractPhotosynthesisMechanism

Helper struct for dispatching between C3 and C4 photosynthesis.
"""
struct C4 <: AbstractPhotosynthesisMechanism end

abstract type AbstractPhotosynthesisModel{FT} <: AbstractCanopyComponent{FT} end
"""
    FarquharParameters{FT<:AbstractFloat, MECH <: AbstractPhotosynthesisMechanism, SP <: SIFParameters}

The required parameters for the Farquhar photosynthesis model.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct FarquharParameters{
    FT <: AbstractFloat,
    MECH <: AbstractPhotosynthesisMechanism,
    SP <: SIFParameters,
}
    "Vcmax at 25 °C (mol CO2/m^2/s)"
    Vcmax25::FT
    "Γstar at 25 °C (mol/mol)"
    Γstar25::FT
    "Michaelis-Menten parameter for CO2 at 25 °C (mol/mol)"
    Kc25::FT
    "Michaelis-Menten parameter for O2 at 25 °C (mol/mol)"
    Ko25::FT
    "Energy of activation for CO2 (J/mol)"
    ΔHkc::FT
    "Energy of activation for oxygen (J/mol)"
    ΔHko::FT
    "Energy of activation for Vcmax (J/mol)"
    ΔHVcmax::FT
    "Energy of activation for Γstar (J/mol)"
    ΔHΓstar::FT
    "Energy of activation for Jmax (J/mol)"
    ΔHJmax::FT
    "Energy of activation for Rd (J/mol)"
    ΔHRd::FT
    "Reference temperature equal to 25 degrees Celsius (K)"
    To::FT
    "Intercelluar O2 concentration (mol/mol); taken to be constant"
    oi::FT
    "Quantum yield of photosystem II (Bernacchi, 2003; unitless)"
    ϕ::FT
    "Curvature parameter, a fitting constant to compute J, unitless"
    θj::FT
    "Constant factor appearing the dark respiration term, equal to 0.015."
    f::FT
    "Sensitivity to low water pressure, in the moisture stress factor, (Pa^{-1}) [Tuzet et al. (2003)]"
    sc::FT
    "Reference water pressure for the moisture stress factor (Pa) [Tuzet et al. (2003)]"
    pc::FT
    "Photosynthesis mechanism: C3 or C4"
    mechanism::MECH
    "Parameters of the sif model, Lee et al 2015"
    sif_parameters::SP
end

Base.eltype(::FarquharParameters{FT}) where {FT} = FT

struct FarquharModel{FT, FP <: FarquharParameters{FT}} <:
       AbstractPhotosynthesisModel{FT}
    parameters::FP
end

function FarquharModel{FT}(
    parameters::FarquharParameters{FT},
) where {FT <: AbstractFloat}
    return FarquharModel{eltype(parameters), typeof(parameters)}(parameters)
end

ClimaLand.name(model::AbstractPhotosynthesisModel) = :photosynthesis
ClimaLand.auxiliary_vars(model::FarquharModel) =
    (:An, :GPP, :Rd, :Vcmax25, :SIF)
ClimaLand.auxiliary_types(model::FarquharModel{FT}) where {FT} =
    (FT, FT, FT, FT, FT)
ClimaLand.auxiliary_domain_names(::FarquharModel) =
    (:surface, :surface, :surface, :surface, :surface)

function photosynthesis_at_a_point_Farquhar(
    T,
    β,
    Rd,
    APAR,
    c_co2,
    medlyn_factor,
    R,
    parameters,
)
    (;
        Vcmax25,
        Γstar25,
        ΔHJmax,
        ΔHVcmax,
        ΔHΓstar,
        f,
        ΔHRd,
        To,
        θj,
        ϕ,
        mechanism,
        oi,
        Kc25,
        Ko25,
        ΔHkc,
        ΔHko,
    ) = parameters
    Jmax = max_electron_transport(Vcmax25, ΔHJmax, T, To, R)
    J = electron_transport(APAR, Jmax, θj, ϕ)
    Vcmax = compute_Vcmax(Vcmax25, T, To, R, ΔHVcmax)
    Γstar = co2_compensation(Γstar25, ΔHΓstar, T, To, R)
    ci = intercellular_co2(c_co2, Γstar, medlyn_factor)
    Aj = light_assimilation(mechanism, J, ci, Γstar)
    Kc = MM_Kc(Kc25, ΔHkc, T, To, R)
    Ko = MM_Ko(Ko25, ΔHko, T, To, R)
    Ac = rubisco_assimilation(mechanism, Vcmax, ci, Γstar, Kc, Ko, oi)
    return net_photosynthesis(Ac, Aj, Rd, β)
end

"""
    update_photosynthesis!(Rd, An, Vcmax25, SIF,
        model::FarquharModel,
        T,
        APAR,
        β,
        medlyn_factor,
        c_co2,
        R,
    )

Computes the net photosynthesis rate `An` for the Farquhar model, along with the
dark respiration `Rd`, and updates them in place.

To do so, we require the canopy leaf temperature `T`, Medlyn factor, `APAR` in
photons per m^2 per second, CO2 concentration in the atmosphere,
moisture stress factor `β` (unitless), and the universal gas constant
`R`.
"""
function update_photosynthesis!(
    Rd,
    An,
    Vcmax25field,
    SIF,
    model::FarquharModel,
    T,
    APAR,
    β,
    medlyn_factor,
    c_co2,
    R,
)
    (; Vcmax25, f, ΔHRd, To) = model.parameters

    @. Rd = dark_respiration(Vcmax25, β, f, ΔHRd, T, To, R)
    @. An = photosynthesis_at_a_point_Farquhar(
        T,
        β,
        Rd,
        APAR,
        c_co2,
        medlyn_factor,
        R,
        model.parameters,
    )
    Vcmax25field .= Vcmax25
    @. SIF = compute_SIF_at_a_point(APAR, T, Vcmax25, R, model.parameters)

end
Base.broadcastable(m::AbstractPhotosynthesisMechanism) = tuple(m)
Base.broadcastable(m::FarquharParameters) = tuple(m)
Base.broadcastable(m::SIFParameters) = tuple(m)

include("./optimality_farquhar.jl")
