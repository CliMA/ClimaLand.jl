export FarquharParameters, FarquharModel, C3, C4

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
    FarquharParameters{FT<:AbstractFloat, MECH <: AbstractPhotosynthesisMechanism}

The required parameters for the Farquhar photosynthesis model.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct FarquharParameters{
    FT <: AbstractFloat,
    MECH <: AbstractPhotosynthesisMechanism,
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
ClimaLand.auxiliary_vars(model::FarquharModel) = (:An, :GPP, :Rd, :Vcmax25)
ClimaLand.auxiliary_types(model::FarquharModel{FT}) where {FT} =
    (FT, FT, FT, FT)
ClimaLand.auxiliary_domain_names(::FarquharModel) =
    (:surface, :surface, :surface, :surface)

"""
    update_photosynthesis!(Rd, An, Vcmax25,
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
    model::FarquharModel,
    T,
    APAR,
    β,
    medlyn_factor,
    c_co2,
    R,
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
    ) = model.parameters
    Jmax = max_electron_transport.(Vcmax25, ΔHJmax, T, To, R)
    J = electron_transport.(APAR, Jmax, θj, ϕ)
    Vcmax = compute_Vcmax.(Vcmax25, T, To, R, ΔHVcmax)
    Γstar = co2_compensation.(Γstar25, ΔHΓstar, T, To, R)
    ci = intercellular_co2.(c_co2, Γstar, medlyn_factor)
    Aj = light_assimilation.(mechanism, J, ci, Γstar)
    Kc = MM_Kc.(Kc25, ΔHkc, T, To, R)
    Ko = MM_Ko.(Ko25, ΔHko, T, To, R)
    Ac = rubisco_assimilation.(mechanism, Vcmax, ci, Γstar, Kc, Ko, oi)
    @. Rd = dark_respiration(Vcmax25, β, f, ΔHRd, T, To, R)
    @. An = net_photosynthesis(Ac, Aj, Rd, β)
    Vcmax25field .= Vcmax25
end
Base.broadcastable(m::AbstractPhotosynthesisMechanism) = tuple(m) # this is so that @. does not broadcast on Ref(canopy.autotrophic_respiration)

include("./optimality_farquhar.jl")
