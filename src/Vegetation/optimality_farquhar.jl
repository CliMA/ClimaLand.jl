export OptimalityFarquharParameters, OptimalityFarquharModel


"""
    OptimalityFarquharParameters{FT<:AbstractFloat}

The required parameters for the optimality Farquhar photosynthesis model.
$(DocStringExtensions.FIELDS)
"""
struct OptimalityFarquharParameters{FT <: AbstractFloat}
    "Photosynthesis mechanism: C3 or C4"
    mechanism::AbstractPhotosynthesisMechanism
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
    ψc::FT
    "Ratio of Jmax25 to Vcmax25 (unitless)"
    jv_ratio::FT
end

"""
    function OptimalityFarquharParameters{FT}(mechanism::AbstractPhotosynthesisMechanism;
        oi = FT(0.209),# mol/mol
        ϕ = FT(0.6), # unitless
        θj = FT(0.9), # unitless
        f = FT(0.015), # unitless
        sc = FT(5e-6),# Pa
        ψc = FT(-2e6), # Pa
        Γstar25 = FT(4.275e-5),  # converted from 42.75 μmol/mol to mol/mol
        Kc25 = FT(4.049e-4), # converted from 404.9 μmol/mol to mol/mol
        Ko25 = FT(0.2874), # converted from 278.4 mmol/mol to mol/mol
        To = FT(298.15), # 25 C
        ΔHkc = FT(79430), #J/mol, Table 11.2 Bonan
        ΔHko = FT(36380), #J/mol, Table 11.2 Bonan
        ΔHVcmax = FT(58520), #J/mol, Table 11.2 Bonan
        ΔHΓstar = FT(37830), #J/mol, 11.2 Bonan
        ΔHJmax = FT(43540), # J/mol, 11.2 Bonan
        ΔHRd = FT(43390), # J/mol, 11.2 Bonan
        jv_ratio = FT(1.67), # unitless, Table 11.5 Bonan
        ) where {FT}

A constructor supplying default values for the FarquharParameters struct.
"""
function OptimalityFarquharParameters{FT}(
    mechanism::AbstractPhotosynthesisMechanism;
    oi = FT(0.209),
    ϕ = FT(0.6),
    θj = FT(0.9),
    f = FT(0.015),
    sc = FT(5e-6),
    ψc = FT(-2e6),
    Γstar25 = FT(4.275e-5),
    Kc25 = FT(4.049e-4),
    Ko25 = FT(0.2874),
    To = FT(298.15),
    ΔHkc = FT(79430),
    ΔHko = FT(36380),
    ΔHVcmax = FT(58520),
    ΔHΓstar = FT(37830),
    ΔHJmax = FT(43540),
    ΔHRd = FT(46390),
    jv_ratio = FT(1.67),
) where {FT}
    return OptimalityFarquharParameters{FT}(
        mechanism,
        Γstar25,
        Kc25,
        Ko25,
        ΔHkc,
        ΔHko,
        ΔHVcmax,
        ΔHΓstar,
        ΔHJmax,
        ΔHRd,
        To,
        oi,
        ϕ,
        θj,
        f,
        sc,
        ψc,
        jv_ratio,
    )
end

"""
    OptimalityFarquharModel{FT} <: AbstractPhotosynthesisModel{FT}

Optimality model of Smith et al. (2019) for estimating Vcmax, based on the assumption that Aj = Ac.

Smith et al. (2019). Global photosynthetic capacity is optimized to the environment. Ecology Letters, 22(3), 506–517. https://doi.org/10.1111/ele.13210
"""
struct OptimalityFarquharModel{FT} <: AbstractPhotosynthesisModel{FT}
    "Required parameters for the Optimality based Farquhar model of Smith et al."
    parameters::OptimalityFarquharParameters{FT}
end

ClimaLSM.auxiliary_vars(model::OptimalityFarquharModel) = (:An, :GPP)
ClimaLSM.auxiliary_types(model::OptimalityFarquharModel{FT}) where {FT} = (FT, FT)

"""
    compute_photosynthesis(
        model::OptimalityFarquharModel,
        T,
        medlyn_factor,
        APAR,
        c_co2,
        β,
        R,
    )

Computes the net photosynthesis rate for the Farquhar optimality model,
given the canopy leaf temperature `T`, Medlyn factor, `APAR` in
photons per m^2 per second, CO2 concentration in the atmosphere,
moisture stress factor `beta` (unitless), and the universal gas constant
`R`.
"""
function compute_photosynthesis(
    model::OptimalityFarquharModel,
    T,
    medlyn_factor,
    APAR,
    c_co2,
    β,
    R,
)
    (;
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
        jv_ratio,
    ) = model.parameters
    Γstar = co2_compensation(Γstar25, ΔHΓstar, T, To, R)
    ci = intercellular_co2(c_co2, Γstar, medlyn_factor)
    Kc = MM_Kc(Kc25, ΔHkc, T, To, R)
    Ko = MM_Ko(Ko25, ΔHko, T, To, R)

    # Light utilization of APAR
    IPSII = ϕ * APAR / 2

    Jmax_arr = arrhenius_function(T, To, R, ΔHJmax)
    Vcmax_arr = arrhenius_function(T, To, R, ΔHVcmax)
    Vcmax_over_Jmax = Vcmax_arr/Jmax_arr/jv_ratio
    if typeof(mechanism) == C3
        C = Vcmax_over_Jmax*(4*(ci + 2*Γstar)/(ci + Kc*(1 + oi/Ko)))
    elseif typeof(mechanism) == C4
        C = Vcmax_over_Jmax
    end

    Jmax = (C - 1)*IPSII/(C*(C*θj - 1))
    J = electron_transport(APAR, Jmax, θj, ϕ)
    Jmax25 = Jmax/Jmax_arr
    Vcmax25 = compute_Vcmax25(Jmax25, jv_ratio)
    Vcmax = Vcmax_over_Jmax*Jmax
    Aj = light_assimilation(mechanism, J, ci, Γstar)
    Ac = rubisco_assimilation(mechanism, Vcmax, ci, Γstar, Kc, Ko, oi)
    @assert Ac ≈ Aj
    Rd = dark_respiration(Vcmax25, β, f, ΔHRd, T, To, R)
    return net_photosynthesis(Ac, Aj, Rd, β)
end