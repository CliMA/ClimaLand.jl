export PModelParameters, PModelDrivers, PModelConstants, compute_pmodel_outputs, PModelModel

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
    "Absorbed photosynthetically active radiation (PAR) (mol m^-2 s^-1)"
    f_abs::FT
    "Ambient CO2 partial pressure (Pa)"
    ca::FT
    "Ambient air pressure (Pa)"
    P_air::FT
    "Downwelling PAR (mol m^-2 s^-1)"
    par_d::FT
    "Vapor pressure deficit (Pa)"
    VPD::FT
end

# TODO: add descriptions to these
Base.@kwdef struct PModelConstants{FT}
    R::FT
    Kc25::FT
    Ko25::FT
    To::FT
    ΔHkc::FT
    ΔHko::FT
    ΔHVcmax::FT
    ΔHRd::FT
    fC3::FT
    Drel::FT
    energy_per_mole_photon_par::FT
end

function create_pmodel_constants(canopy, FT)
    earth_param_set = canopy.parameters.earth_param_set
    lightspeed = LP.light_speed(earth_param_set)
    planck_h = LP.planck_constant(earth_param_set)
    N_a = LP.avogadro_constant(earth_param_set)
    (; _, λ_γ_PAR, Ω) = canopy.radiative_transfer.parameters
    energy_per_mole_photon_par = planck_h * lightspeed / λ_γ_PAR * N_a

    return PModelConstants(
        R = LP.gas_constant(earth_param_set),
        Kc25 = FT(39.97),
        Ko25 = FT(27480.0),
        To = FT(298.15),
        ΔHkc = FT(79430.0),
        ΔHko = FT(36380.0),
        ΔHVcmax = FT(58520.0),
        ΔHRd = FT(46390.0),
        fC3 = FT(0.015),
        Drel = FT(1.6),
        energy_per_mole_photon_par = FT(energy_per_mole_photon_par),
        Ω = FT(Ω)
    )
end

"""
    PModelModel{FT,
                OPFT <: PModelParameters{FT}
                } <: AbstractPhotosynthesisModel{FT}

Optimality model of Smith et al. (2019) for estimating Vcmax, based on the assumption that Aj = Ac.
Smith et al. (2019). Global photosynthetic capacity is optimized to the environment. Ecology Letters, 22(3), 506–517. https://doi.org/10.1111/ele.13210
"""
struct PModelModel{FT, OPFT <: PModelParameters{FT}} <:
       AbstractPhotosynthesisModel{FT}
    "Required parameters for the P-model of Stocker et al. (2020)"
    parameters::OPFT
end

function PModelModel{FT}(
    parameters::PModelParameters{FT},
) where {FT <: AbstractFloat}
    return PModelModel{eltype(parameters), typeof(parameters)}(
        parameters,
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
    (; T_canopy, f_abs, ca, P_air, par_d, VPD) = drivers

    # Unpack constants
    (; R, Kc25, Ko25, To, ΔHkc, ΔHko, 
        ΔHVcmax, ΔHRd, fC3, Drel, energy_per_mole_photon_par, Ω) = constants

    # Compute intermediate values
    ϕ0 = isnan(ϕ0) ? intrinsic_quantum_yield(T_canopy, ϕc) : ϕ0

    Γstar = co2_compensation_p(T_canopy, P_air, R)
    ηstar = compute_viscosity_ratio(T_canopy, P_air)
    Kmm = compute_Kmm(T_canopy, P_air, Kc25, Ko25, ΔHkc, ΔHko, To, R)
    I_abs = f_abs * par_d / energy_per_mole_photon_par
    χ, ξ, mj, mc = optimal_co2_ratio_c3(Kmm, Γstar, ηstar, ca, VPD, β, Drel)
    ci = χ * ca
    mprime = compute_mj_with_jmax_limitation(mj, cstar)
    Vcmax = pmodel_vcmax(ϕ0, I_abs, mprime, mc)
    Ac = Vcmax * mc

    Aj = Ac # Coordination hypothesis
    LUE = compute_LUE(ϕ0, FT(1.0), mprime)
    Vcmax25 = Vcmax / arrhenius_function(T_canopy, To, R, ΔHVcmax)
    Rd = dark_respiration(FT(1.0), Vcmax25, FT(1.0), T_canopy, R, To, fC3, ΔHRd)
    An = net_photosynthesis(Ac, Aj, Rd, FT(1.0))

    # I noticed that the GPP definition in optimality_farquhar.jl is different from 
    # the one in Stocker et al. (2020). It uses LAI, Ω, and extinction coefficients.
    # Here, we use the simpler definition from Stocker et al. (2020). but TODO is 
    # to figure out the discrepancy 
    GPP = I_abs * LUE 

    iWUE = (ca - ci) / Drel
    gs = An / (ca - ci)
    Jmax = pmodel_jmax(ϕ0, I_abs, mprime)

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
        "rd" => Rd,
        "An" => An
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
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    f_abs = p.canopy.radiative_transfer.par.abs
    ca = p.drivers.c_co2
    P_air = p.drivers.P
    par_d = p.canopy.radiative_transfer.par_d
    VPD = ClimaLand.vapor_pressure_deficit(
        p.drivers.T, p.drivers.P, p.drivers.q, canopy.parameters.earth_param_set.thermo_params
    )

    # Create PModelDrivers
    drivers = PModelDrivers(
        T_canopy = T_canopy,
        f_abs = f_abs,
        ca = ca,
        P_air = P_air,
        par_d = par_d,
        VPD = VPD
    )

    # Create PModelConstants
    constants = create_pmodel_constants(canopy, eltype(p.canopy.photosynthesis.Rd))

    parameters = model.parameters

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
