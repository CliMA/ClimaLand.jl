export FarquharParameters, FarquharModel

abstract type AbstractPhotosynthesisModel{FT} <: AbstractCanopyComponent{FT} end

"""
    FarquharParameters{
        FT<:AbstractFloat,
        MECH <: Union{FT, ClimaCore.Fields.Field},
        VC <: Union{FT, ClimaCore.Fields.Field},
    }

The required parameters for the Farquhar photosynthesis model.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct FarquharParameters{
    FT <: AbstractFloat,
    MECH <: Union{FT, ClimaCore.Fields.Field},
    VC <: Union{FT, ClimaCore.Fields.Field},
}
    "Vcmax at 25 °C (mol CO2/m^2/s)"
    Vcmax25::VC
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
    "Constant factor appearing the dark respiration term for C3, equal to 0.015."
    fC3::FT
    "Constant factor appearing the dark respiration term for C4, equal to 0.025."
    fC4::FT
    "Sensitivity to low water pressure, in the moisture stress factor, (Pa^{-1}) [Tuzet et al. (2003)]"
    sc::FT
    "Reference water pressure for the moisture stress factor (Pa) [Tuzet et al. (2003)]"
    pc::FT
    "Q10 temperature parameter for Vcmax and Rd for C4 photosynthesis; unitless"
    Q10::FT
    "Parameter appearing in temperature dependence of C4 Vcmax; K^{-1}"
    s1::FT
    "Parameter appearing in temperature dependence of C4 Vcmax; K"
    s2::FT
    "Parameter appearing in temperature dependence of C4 Vcmax; K^{-1}"
    s3::FT
    "Parameter appearing in temperature dependence of C4 Vcmax; K"
    s4::FT
    "Parameter appearing in temperature dependence of C4 Rd; K^{-1}"
    s5::FT
    "Parameter appearing in temperature dependence of C4 Rd; K"
    s6::FT
    "Quantum yield for C4 photosynthesis; mol/mol"
    E::FT
    "Photosynthesis mechanism: 1.0 indicates C3, 0.0 indicates C4"
    is_c3::MECH
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


function photosynthesis_at_a_point_Farquhar(
    T,
    β,
    Rd,
    APAR,
    c_co2,
    medlyn_factor,
    R,
    Vcmax25,
    is_c3,
    Γstar25,
    ΔHJmax,
    ΔHVcmax,
    ΔHΓstar,
    fC3,
    fC4,
    ΔHRd,
    To,
    θj,
    ϕ,
    oi,
    Kc25,
    Ko25,
    ΔHkc,
    ΔHko,
    Q10,
    s1,
    s2,
    s3,
    s4,
    s5,
    s6,
    E,
)
    Jmax = max_electron_transport(Vcmax25, ΔHJmax, T, To, R)
    J = electron_transport(APAR, Jmax, θj, ϕ)
    Vcmax =
        compute_Vcmax(is_c3, Vcmax25, T, R, To, ΔHVcmax, Q10, s1, s2, s3, s4)
    Γstar = co2_compensation(Γstar25, ΔHΓstar, T, To, R)
    ci = intercellular_co2(c_co2, Γstar, medlyn_factor)
    Aj = light_assimilation(is_c3, J, ci, Γstar, APAR, E)
    Kc = MM_Kc(Kc25, ΔHkc, T, To, R)
    Ko = MM_Ko(Ko25, ΔHko, T, To, R)
    Ac = rubisco_assimilation(is_c3, Vcmax, ci, Γstar, Kc, Ko, oi)
    return net_photosynthesis(Ac, Aj, Rd, β)
end

"""
    update_photosynthesis!(
        Rd,
        An,
        Vcmax25field,
        model::FarquharModel,
        T,
        f_abs,
        β,
        medlyn_factor,
        c_co2,
        R,
        energy_per_mole_photon_par,
        par_d,
)

Computes the net photosynthesis rate `An` for the Farquhar model, along with the
dark respiration `Rd`, and updates them in place.

To do so, we require the canopy leaf temperature `T`, Medlyn factor, fraction
of `par_d` aborbed `f_abs`, CO2 concentration in the atmosphere,
moisture stress factor `β` (unitless), and the universal gas constant
`R`.

The typical `energy_per_mole_photon_par` is used to convert from an absorbed energy
flux to a flux of moles of photons, as needed by photosynthetic rate computations.
"""
function update_photosynthesis!(
    Rd,
    An,
    Vcmax25field,
    model::FarquharModel,
    T,
    f_abs,
    β,
    medlyn_factor,
    c_co2,
    R,
    energy_per_mole_photon_par,
    par_d,
)
    (;
        Vcmax25,
        is_c3,
        Γstar25,
        ΔHJmax,
        ΔHVcmax,
        ΔHΓstar,
        fC3,
        fC4,
        ΔHRd,
        To,
        θj,
        ϕ,
        oi,
        Kc25,
        Ko25,
        ΔHkc,
        ΔHko,
        Q10,
        s1,
        s2,
        s3,
        s4,
        s5,
        s6,
        E,
    ) = model.parameters
    Vcmax25field .= Vcmax25
    @. Rd = dark_respiration(
        is_c3,
        Vcmax25,
        β,
        T,
        R,
        To,
        fC3,
        ΔHRd,
        Q10,
        s5,
        s6,
        fC4,
    )
    @. An = photosynthesis_at_a_point_Farquhar(
        T,
        β,
        Rd,
        f_abs * par_d / energy_per_mole_photon_par, # This function requires flux in moles of photons, not J
        c_co2,
        medlyn_factor,
        R,
        Vcmax25,
        is_c3,
        Γstar25,
        ΔHJmax,
        ΔHVcmax,
        ΔHΓstar,
        fC3,
        fC4,
        ΔHRd,
        To,
        θj,
        ϕ,
        oi,
        Kc25,
        Ko25,
        ΔHkc,
        ΔHko,
        Q10,
        s1,
        s2,
        s3,
        s4,
        s5,
        s6,
        E,
    )

end
Base.broadcastable(m::FarquharParameters) = tuple(m)

include("./optimality_farquhar.jl")
