export PModelParameters,
    PModelDrivers, PModelConstants, compute_full_pmodel_outputs, PModel

"""
    PModelParameters{FT<:AbstractFloat}

The required parameters for P-model (Stocker et al. 2020). Parameters are typically
tunable with considerable uncertainty. 
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PModelParameters{FT <: AbstractFloat}
    "Constant describing cost of maintaining electron transport (unitless)
    Typical value = 0.41"
    cstar::FT
    "Ratio of unit costs of transpiration and carboxylation (unitless)
    Typical value = 146"
    β::FT
    "Scaling parameter for temp-dependent intrinsic quantum yield (unitless)
    Typical value = 0.087"
    ϕc::FT
    "Temp-independent intrinsic quantum yield. If provided, overrides ϕc. (unitless)
    Typical value = 0.05"
    ϕ0::FT
    """Constant term in temp-dependent intrinsic quantum yield (unitless)."""
    ϕa0::FT
    """First order term in temp-dependent intrinsic quantum yield (K^-1)."""
    ϕa1::FT
    """Second order term in temp-dependent intrinsic quantum yield (K^-2)."""
    ϕa2::FT
end

"""
    PModelDrivers{FT<:AbstractFloat}

The required drivers for P-model (Stocker et al. 2020). Drivers are defined as
external variables used to compute the optimal photosynthetic capacities. 
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PModelDrivers{FT <: AbstractFloat}
    "Canopy temperature (K)"
    T_canopy::FT
    "Absorbed PAR in moles of photons (mol m^-2 s^-1)"
    I_abs::FT
    "Ambient CO2 partial pressure (Pa)"
    ca::FT
    "Ambient air pressure (Pa)"
    P_air::FT
    "Vapor pressure deficit (Pa)"
    VPD::FT
    """Soil moisture stress factor (unitless)"""
    βm::FT
end


"""
    PModelConstants{FT<:AbstractFloat}

The required constants for P-model (Stocker et al. 2020). These are physical
or biochemical constants that are not expected to change with time or space.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PModelConstants{FT}
    """Gas constant (J mol^-1 K^-1)"""
    R::FT
    """Michaelis-Menten parameter for carboxylation at 25°C (μmol mol^-1)"""
    Kc25::FT
    """Michaelis-Menten parameter for oxygenation at 25°C (μmol mol^-1)"""
    Ko25::FT
    """Reference temperature equal to 25˚C (K)"""
    To::FT
    """Energy of activation for Kc (J mol^-1)"""
    ΔHkc::FT
    """Energy of activation for Ko (J mol^-1)"""
    ΔHko::FT
    """Relative diffusivity of CO2 in the stomatal pores, equal to 1.6."""
    Drel::FT
    """Effective energy of activation for Γstar (J mol^-1)"""
    ΔHΓstar::FT
    """Γstar at 25 °C (Pa)"""
    Γstar25::FT
    """Effective energy of activation for Vcmax (J mol^-1)"""
    Ha_Vcmax::FT
    """Effective energy of deactivation for Vcmax (J mol^-1)"""
    Hd_Vcmax::FT
    """Intercept term for dS in Vcmax deactivation factor (J K^-1 mol^-1)"""
    aS_Vcmax::FT
    """Slope term for dS in Vcmax deactivation factor (J K^-2 mol^-1)"""
    bS_Vcmax::FT
    """Effective energy of activation for Jmax (J mol^-1)"""
    Ha_Jmax::FT
    """Effective energy of deactivation for Jmax (J mol^-1)"""
    Hd_Jmax::FT
    """Intercept term for dS in Jmax deactivation factor (J K^-1 mol^-1)"""
    aS_Jmax::FT
    """Slope term for dS in Jmax deactivation factor (J K^-2 mol^-1)"""
    bS_Jmax::FT
    """Molar mass of carbon (kg mol^-1)"""
    Mc::FT
    """Intercellular O2 mixing ratio (unitless)"""
    oi::FT
    """First order coefficient for temp-dependent Rd (K^-1)"""
    aRd::FT
    """Second order coefficient for temp-dependent Rd (K^-2)"""
    bRd::FT
    """Constant factor appearing the dark respiration term for C3 plants (unitless)"""
    fC3::FT
end

Base.eltype(::PModelParameters{FT}) where {FT} = FT
Base.eltype(::PModelDrivers{FT}) where {FT} = FT
Base.eltype(::PModelConstants{FT}) where {FT} = FT

Base.broadcastable(m::PModelParameters) = tuple(m)
Base.broadcastable(m::PModelDrivers) = tuple(m)
Base.broadcastable(m::PModelConstants) = tuple(m)

"""
    PModelConstants(FT)

Creates a `PModelConstants` object with default values for the P-model constants.
"""
function PModelConstants(FT)
    return PModelConstants(
        R = LP.gas_constant(LP.LandParameters(FT)),
        Kc25 = FT(39.97),
        Ko25 = FT(27480),
        To = FT(298.15),
        ΔHkc = FT(79430),
        ΔHko = FT(36380),
        Drel = FT(1.6),
        ΔHΓstar = FT(37830),
        Γstar25 = FT(4.332),
        Ha_Vcmax = FT(71513),
        Hd_Vcmax = FT(200000),
        aS_Vcmax = FT(668.39),
        bS_Vcmax = FT(1.07),
        Ha_Jmax = FT(49884),
        Hd_Jmax = FT(200000),
        aS_Jmax = FT(659.70),
        bS_Jmax = FT(0.75),
        Mc = FT(0.0120107),
        oi = FT(0.2095),
        aRd = FT(0.1012),
        bRd = FT(-0.0005),
        fC3 = FT(0.015),
    )
end

"""
    PModel{FT,
                OPFT <: PModelParameters{FT},
                OPCT <: PModelConstants{FT}
                } <: AbstractPhotosynthesisModel{FT}
"""
struct PModel{FT, OPFT <: PModelParameters{FT}, OPCT <: PModelConstants{FT}} <:
       AbstractPhotosynthesisModel{FT}
    "Required parameters for the P-model of Stocker et al. (2020)"
    parameters::OPFT
    "Constants for the P-model"
    constants::OPCT
end

function PModel{FT}(
    parameters::PModelParameters{FT},
    constants::PModelConstants{FT} = PModelConstants(FT),
) where {FT <: AbstractFloat}
    return PModel{FT, typeof(parameters), typeof(constants)}(
        parameters,
        constants,
    )
end

ClimaLand.auxiliary_vars(model::PModel) = (:An, :GPP, :Rd, :Vcmax25)
ClimaLand.auxiliary_types(model::PModel{FT}) where {FT} = (FT, FT, FT, FT)
ClimaLand.auxiliary_domain_names(::PModel) =
    (:surface, :surface, :surface, :surface)


"""
    compute_full_pmodel_outputs(
        parameters::PModelParameters, 
        drivers::PModelDrivers, 
        constants::PModelConstants
    )

Performs the P-model computations as defined in Stocker et al. (2020) 
and returns a dictionary of outputs. See https://github.com/geco-bern/rpmodel
for a code reference.

Output name      Description (units)
    "gpp"           Gross primary productivity per leaf area (kg C m^-2 s^-1)
    "gammastar"     CO2 compensation point (Pa)
    "kmm"           Effective MM coefficient for Rubisco-limited photosynthesis (Pa)
    "ns_star"       Viscosity of water normalized to 25 deg C (unitless)
    "chi"           Optimal ratio of intercellular to ambient CO2 (unitless) 
    "xi"            Sensitivity of χ to VPD (Pa^1/2)
    "mj"            CO2 limitation factor for light-limited photosynthesis (unitless)
    "mc"            CO2 limitation factor for Rubisco-limited photosynthesis (unitless)
    "ci"            Intercellular CO2 concentration (Pa)
    "iwue"          Intrinsic water use efficiency (Pa)
    "gs"            Stomatal conductance (mol m^-2 s^-1 Pa^-1)
    "vcmax"         Maximum rate of carboxlation (mol m^-2 s^-1)
    "vcmax25"       Vcmax normalized to 25°C via modified-Arrhenius type function (mol m^-2 s^-1)
    "jmax"          Maximum rate of electron transport (mol m^-2 s^-1)
    "jmax25"        Jmax normalized to 25°C via modified-Arrhenius type function (mol m^-2 s^-1)

"""
function compute_full_pmodel_outputs(
    parameters::PModelParameters{FT},
    drivers::PModelDrivers{FT},
    constants::PModelConstants{FT},
) where {FT}
    # Unpack parameters
    (; cstar, β, ϕc, ϕ0, ϕa0, ϕa1, ϕa2) = parameters

    # Unpack drivers
    (; T_canopy, I_abs, ca, P_air, VPD, βm) = drivers

    # Unpack constants
    (;
        R,
        Kc25,
        Ko25,
        To,
        ΔHkc,
        ΔHko,
        Drel,
        ΔHΓstar,
        Γstar25,
        Ha_Vcmax,
        Hd_Vcmax,
        aS_Vcmax,
        bS_Vcmax,
        Ha_Jmax,
        Hd_Jmax,
        aS_Jmax,
        bS_Jmax,
        Mc,
        oi,
        aRd,
        bRd,
        fC3,
    ) = constants

    # Compute intermediate values
    ϕ0 = isnan(ϕ0) ? intrinsic_quantum_yield(T_canopy, ϕc, ϕa0, ϕa1, ϕa2) : ϕ0

    Γstar = co2_compensation_p(T_canopy, To, P_air, R, ΔHΓstar, Γstar25)
    ηstar = compute_viscosity_ratio(T_canopy, P_air, true)
    Kmm = compute_Kmm(T_canopy, P_air, Kc25, Ko25, ΔHkc, ΔHko, To, R, oi)
    χ, ξ, mj, mc = optimal_co2_ratio_c3(Kmm, Γstar, ηstar, ca, VPD, β, Drel)
    ci = χ * ca
    mprime = compute_mj_with_jmax_limitation(mj, cstar)

    Vcmax = βm * ϕ0 * I_abs * mprime / mc
    inst_temp_scaling_vcmax25 = inst_temp_scaling(
        T_canopy,
        T_canopy,
        To,
        Ha_Vcmax,
        Hd_Vcmax,
        aS_Vcmax,
        bS_Vcmax,
        R,
    )
    Vcmax25 = Vcmax / inst_temp_scaling_vcmax25

    Jmax = FT(4) * ϕ0 * I_abs / sqrt((mj / (βm * mprime))^2 - 1)
    Jmax25 =
        Jmax / inst_temp_scaling(
            T_canopy,
            T_canopy,
            To,
            Ha_Jmax,
            Hd_Jmax,
            aS_Jmax,
            bS_Jmax,
            R,
        )
    J = electron_transport_pmodel(ϕ0, I_abs, Jmax)

    Ac = Vcmax * mc
    Aj = J * mj / FT(4)

    LUE = compute_LUE(ϕ0, βm, mprime, Mc)
    GPP = I_abs * LUE

    # intrinsic water use efficiency (iWUE) and stomatal conductance (gs)
    iWUE = (ca - ci) / Drel
    gs = pmodel_gs(χ, ca, Ac)

    # dark respiration 
    rd = fC3 * inst_temp_scaling_rd(T_canopy, To, aRd, bRd) * Vcmax25

    return (;
        gpp = GPP,
        gammastar = Γstar,
        kmm = Kmm,
        ca = ca,
        ns_star = ηstar,
        chi = χ,
        xi = ξ,
        mj = mj,
        mc = mc,
        ci = ci,
        iwue = iWUE,
        gs = gs,
        vcmax = Vcmax,
        vcmax25 = Vcmax25,
        jmax = Jmax,
        jmax25 = Jmax25,
        rd = rd,
    )
end


get_Vcmax25(p, m::PModelParameters) = p.canopy.photosynthesis.Vcmax25
