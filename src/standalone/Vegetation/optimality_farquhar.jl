export OptimalityFarquharParameters, OptimalityFarquharModel


"""
    OptimalityFarquharParameters{FT<:AbstractFloat}

The required parameters for the optimality Farquhar photosynthesis model.
Currently, only C3 photosynthesis is supported.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct OptimalityFarquharParameters{FT <: AbstractFloat}
    "Photosynthesis mechanism: C3 only"
    mechanism::C3
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
    "Intercellular O2 concentration (mol/mol); taken to be constant"
    oi::FT
    "Quantum yield of photosystem II (Bernacchi, 2003; unitless)"
    ϕ::FT
    "Curvature parameter, a fitting constant to compute J, unitless"
    θj::FT
    "Constant factor appearing the dark respiration term, equal to 0.015."
    f::FT
    "Fitting constant to compute the moisture stress factor (Pa^{-1})"
    sc::FT
    "Fitting constant to compute the moisture stress factor (Pa)"
    pc::FT
    "Constant describing cost of maintaining electron transport (unitless)"
    c::FT
end

Base.eltype(::OptimalityFarquharParameters{FT}) where {FT} = FT

"""
    OptimalityFarquharModel{FT,
                            OPFT <: OptimalityFarquharParameters{FT}
                            } <: AbstractPhotosynthesisModel{FT}

Optimality model of Smith et al. (2019) for estimating Vcmax, based on the assumption that Aj = Ac.
Smith et al. (2019). Global photosynthetic capacity is optimized to the environment. Ecology Letters, 22(3), 506–517. https://doi.org/10.1111/ele.13210
"""
struct OptimalityFarquharModel{FT, OPFT <: OptimalityFarquharParameters{FT}} <:
       AbstractPhotosynthesisModel{FT}
    "Required parameters for the Optimality based Farquhar model of Smith et al. (2019)"
    parameters::OPFT
end

function OptimalityFarquharModel{FT}(
    parameters::OptimalityFarquharParameters{FT},
) where {FT <: AbstractFloat}
    return OptimalityFarquharModel{eltype(parameters), typeof(parameters)}(
        parameters,
    )
end

ClimaLand.auxiliary_vars(model::OptimalityFarquharModel) =
    (:An, :GPP, :Rd, :Vcmax25)
ClimaLand.auxiliary_types(model::OptimalityFarquharModel{FT}) where {FT} =
    (FT, FT, FT, FT)
ClimaLand.auxiliary_domain_names(::OptimalityFarquharModel) =
    (:surface, :surface, :surface, :surface)

"""
    update_photosynthesis!(Rd, An, Vcmax25,
        model::OptimalityFarquharModel,
        T,
        APAR,
        β,
        medlyn_factor,
        c_co2,
        R,
    )

Computes the net photosynthesis rate `An` for the Optimality Farquhar model, along with the
dark respiration `Rd`, and the value of `Vcmax25`, and updates them in place.
        
 To do so, we require the canopy leaf temperature `T`, Medlyn factor, `APAR` in
photons per m^2 per second, CO2 concentration in the atmosphere,
moisture stress factor `β` (unitless), and the universal gas constant
`R`.
"""
function update_photosynthesis!(
    Rd,
    An,
    Vcmax25,
    model::OptimalityFarquharModel,
    T,
    APAR,
    β,
    medlyn_factor,
    c_co2,
    R,
)
    (;
        Γstar25,
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
        c,
    ) = model.parameters

    Γstar = co2_compensation.(Γstar25, ΔHΓstar, T, To, R)
    ci = intercellular_co2.(c_co2, Γstar, medlyn_factor)# may change?
    Kc = MM_Kc.(Kc25, ΔHkc, T, To, R)
    Ko = MM_Ko.(Ko25, ΔHko, T, To, R)
    rates =
        optimality_max_photosynthetic_rates.(
            APAR,
            θj,
            ϕ,
            oi,
            ci,
            Γstar,
            Kc,
            Ko,
            c,
        )
    Jmax = rates.:1
    Vcmax = rates.:2
    J = electron_transport.(APAR, Jmax, θj, ϕ)
    Aj = light_assimilation.(mechanism, J, ci, Γstar)
    Ac = rubisco_assimilation.(mechanism, Vcmax, ci, Γstar, Kc, Ko, oi)

    @. Vcmax25 = Vcmax / arrhenius_function(T, To, R, ΔHVcmax)
    @. Rd = dark_respiration(Vcmax25, β, f, ΔHRd, T, To, R)
    @. An = net_photosynthesis(Ac, Aj, Rd, β)
end
