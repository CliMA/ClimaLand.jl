export PModelParameters, 
    PModelDrivers, 
    PModelConstants, 
    create_pmodel_constants, 
    compute_pmodel_outputs, 
    PModelModel

"""
    PModelParameters{FT<:AbstractFloat}

The required parameters for P-model (Stocker et al. 2020).
Currently, only C3 photosynthesis is supported.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PModelParameters{
    FT <: AbstractFloat
}
    "Fitting constant to compute the moisture stress factor (Pa^{-1})"
    sc::FT
    "Fitting constant to compute the moisture stress factor (Pa)"
    pc::FT
    "Constant describing cost of maintaining electron transport (unitless)"
    cstar::FT
    "Ratio of unit costs of transpiration and carboxylation (unitless)"
    β::FT 
    "Scaling parameter for temp-dependent intrinsic quantum yield (unitless)"
    ϕc::FT
    "Temp-independent intrinsic quantum yield (if provided, overrides ϕc)" 
    ϕ0::FT 
    # "Use soil moisture stress in the photosynthesis model (default: true)"
    # use_soil_moisture_stress::Bool = true
end

Base.eltype(::PModelParameters{FT}) where {FT} = FT

Base.@kwdef struct PModelDrivers{
    FT <: AbstractFloat
}
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
    """Effective energy of activation for Vcmax (J mol^-1)"""
    ΔHVcmax::FT
    """Effective energy of activation for Jmax (J mol^-1)"""
    ΔHJmax::FT
    """Effective energy of activation for Rd (J mol^-1)"""
    ΔHRd::FT
    """Constant factor appearing in the dark respiration term for C3 plants, equal to 0.015."""
    fC3::FT
    """Relative diffusivity of CO2 in the stomatal pores, equal to 1.6."""
    Drel::FT
    """Effective energy of activation for Γstar (J mol^-1)"""
    ΔHΓstar::FT
    """Γstar at 25 °C (Pa)"""
    Γstar25::FT
    Ha_Vcmax::FT
    Hd_Vcmax::FT
    aS_Vcmax::FT
    bS_Vcmax::FT
    Ha_Jmax::FT
    Hd_Jmax::FT
    aS_Jmax::FT
    bS_Jmax::FT

end

function create_pmodel_constants(FT)
    return PModelConstants(
        R = LP.gas_constant(LP.LandParameters(FT)),
        Kc25 = FT(39.97),
        Ko25 = FT(27480),
        To = FT(298.15),
        ΔHkc = FT(79430),
        ΔHko = FT(36380),
        ΔHVcmax = FT(65330),
        ΔHJmax = FT(43990),
        ΔHRd = FT(46390),
        fC3 = FT(0.015),
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
        bS_Jmax = FT(0.75)
    )
end

"""
    PModelModel{FT,
                OPFT <: PModelParameters{FT},
                OPCT <: PModelConstants{FT}
                } <: AbstractPhotosynthesisModel{FT}

Optimality model of Smith et al. (2019) for estimating Vcmax, based on the assumption that Aj = Ac.
Smith et al. (2019). Global photosynthetic capacity is optimized to the environment. Ecology Letters, 22(3), 506–517. https://doi.org/10.1111/ele.13210
"""
struct PModelModel{FT, OPFT <: PModelParameters{FT}, OPCT <: PModelConstants{FT}} <:
       AbstractPhotosynthesisModel{FT}
    "Required parameters for the P-model of Stocker et al. (2020)"
    parameters::OPFT
    "Constants for the P-model"
    constants::OPCT
end

function PModelModel{FT}(
    parameters::PModelParameters{FT},
    constants::PModelConstants{FT} = create_pmodel_constants(FT),
) where {FT <: AbstractFloat}
    return PModelModel{FT, typeof(parameters), typeof(constants)}(
        parameters,
        constants,
    )
end

ClimaLand.auxiliary_vars(model::PModelModel) =
    (:An, :GPP, :Rd, :Vcmax25)
ClimaLand.auxiliary_types(model::PModelModel{FT}) where {FT} =
    (FT, FT, FT, FT)
ClimaLand.auxiliary_domain_names(::PModelModel) =
    (:surface, :surface, :surface, :surface)


"""
    compute_pmodel_outputs(
        parameters::PModelParameters, 
        drivers::PModelDrivers, 
        constants::PModelConstants
    )

Performs the P-model computations and returns a dictionary of outputs.
"""
function compute_pmodel_outputs(
    parameters::PModelParameters{FT}, 
    drivers::PModelDrivers{FT}, 
    constants::PModelConstants{FT}
) where {FT}
    # Unpack parameters
    (; sc, pc, cstar, β, ϕc, ϕ0) = parameters

    # Unpack drivers
    (; T_canopy, I_abs, ca, P_air, VPD, βm) = drivers

    # Unpack constants
    (; R, Kc25, Ko25, To, ΔHkc, ΔHko, 
        ΔHVcmax, ΔHJmax, ΔHRd, fC3, Drel, ΔHΓstar, Γstar25,
        Ha_Vcmax, Hd_Vcmax, aS_Vcmax, bS_Vcmax, 
        Ha_Jmax, Hd_Jmax, aS_Jmax, bS_Jmax) = constants

    # Compute intermediate values
    ϕ0 = isnan(ϕ0) ? intrinsic_quantum_yield(T_canopy, ϕc) : ϕ0

    Γstar = co2_compensation_p(T_canopy, To, P_air, R, ΔHΓstar, Γstar25)
    ηstar = compute_viscosity_ratio(T_canopy, P_air)
    Kmm = compute_Kmm(T_canopy, P_air, Kc25, Ko25, ΔHkc, ΔHko, To, R)
    χ, ξ, mj, mc = optimal_co2_ratio_c3(Kmm, Γstar, ηstar, ca, VPD, β, Drel)
    ci = χ * ca
    mprime = compute_mj_with_jmax_limitation(mj, cstar)

    Vcmax = βm * ϕ0 * I_abs * mprime / mc
    Vcmax25 = Vcmax / inst_temp_scaling(T_canopy, T_canopy, To, Ha_Vcmax, Hd_Vcmax, aS_Vcmax, bS_Vcmax, R)

    Jmaxlim = Vcmax * (ci + FT(2) * Γstar) / (ϕ0 * I_abs * (ci + Kmm))
    Jmax = FT(4) * ϕ0 * I_abs / sqrt((FT(1)/Jmaxlim)^2 - FT(1)) 
    Jmax25 = Jmax / inst_temp_scaling(T_canopy, T_canopy, To, Ha_Jmax, Hd_Jmax, aS_Jmax, bS_Jmax, R)
    J = electron_transport_pmodel(ϕ0, I_abs, Jmax)

    Ac = Vcmax * mc
    Aj = J * mj / FT(4)

    @assert isapprox(Ac, Aj; atol=1e-6) "Ac and Aj are not approximately equal: Ac = $Ac, Aj = $Aj"

    LUE = compute_LUE(ϕ0, βm, mprime)

    # I noticed that the GPP definition in optimality_farquhar.jl is different from 
    # the one in Stocker et al. (2020). It uses LAI, Ω, and extinction coefficients.
    # Here, we use the simpler definition from Stocker et al. (2020). but TODO is 
    # to figure out the discrepancy 
    GPP = I_abs * LUE 

    # intrinsic water use efficiency (iWUE) and stomatal conductance (gs)
    iWUE = (ca - ci) / Drel
    gs = Ac / (ca - ci)

    return Dict(
        "gpp" => GPP,
        "ca" => ca,
        "gammastar" => Γstar,
        "kmm" => Kmm,
        "ns_star" => ηstar,
        "chi" => χ,
        "xi" => ξ,
        "mj" => mj,
        "mc" => mc,
        "ci" => ci,
        "iwue" => iWUE,
        "gs" => gs,
        "vcmax" => Vcmax,
        "vcmax25" => Vcmax25,
        "jmax" => Jmax,
        "jmax25" => Jmax25
    )
end


"""
    update_photosynthesis!(p, Y, model::PModelModel, canopy)

Computes the net photosynthesis rate `An` (mol CO2/m^2/s) for the P-model, along with the
dark respiration `Rd` (mol CO2/m^2/s), the value of `Vcmax25` (mol CO2/m^2/s), and the gross primary 
productivity `GPP` (mol CO2/m^2/s), and updates them in place.
"""
function update_photosynthesis!(p, Y, model::PModelModel, canopy)
    # Unpack required fields from `p` and `canopy`
    earth_param_set = canopy.parameters.earth_param_set
    lightspeed = LP.light_speed(earth_param_set)
    planck_h = LP.planck_constant(earth_param_set)
    N_a = LP.avogadro_constant(earth_param_set)
    (; _, λ_γ_PAR, Ω) = canopy.radiative_transfer.parameters
    energy_per_mole_photon_par = planck_h * lightspeed / λ_γ_PAR * N_a

    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    f_abs = p.canopy.radiative_transfer.par.abs
    ca = p.drivers.c_co2
    P_air = p.drivers.P
    par_d = p.canopy.radiative_transfer.par_d
    VPD = ClimaLand.vapor_pressure_deficit(
        p.drivers.T, p.drivers.P, p.drivers.q, canopy.parameters.earth_param_set.thermo_params
    )

    # Calculate I_abs directly
    I_abs = f_abs * par_d / energy_per_mole_photon_par

    # Create PModelDrivers with I_abs
    drivers = PModelDrivers(
        T_canopy = T_canopy,
        I_abs = I_abs,
        ca = ca,
        P_air = P_air,
        VPD = VPD,
        βm = FT(1.0) # TODO: add a new file for soil moisture stress parameterizations
    )

    # Use model's parameters and constants
    parameters = model.parameters
    constants = model.constants

    outputs = compute_pmodel_outputs(parameters, drivers, constants)

    # Update the outputs in place
    p.canopy.photosynthesis.Vcmax25 .= outputs["vcmax25"]
    p.canopy.photosynthesis.Rd .= outputs["rd"]
    p.canopy.photosynthesis.An .= outputs["An"]
    p.canopy.photosynthesis.GPP .= outputs["gpp"]
end

get_Vcmax25(p, m::PModelParameters) =
    p.canopy.photosynthesis.Vcmax25
Base.broadcastable(m::PModelParameters) = tuple(m)
