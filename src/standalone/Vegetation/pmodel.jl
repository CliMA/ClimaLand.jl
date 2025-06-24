export PModelParameters, PModelModel

"""
    PModelParameters{FT<:AbstractFloat}

The required parameters for P-model (Stocker et al. 2020).
Currently, only C3 photosynthesis is supported.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PModelParameters{
    FT <: AbstractFloat,
    MECH <: Union{FT, ClimaCore.Fields.Field},
}
    "Photosynthesis mechanism: C3 only"
    is_c3::MECH
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
    "Energy of activation for Rd (J/mol)"
    ΔHRd::FT
    "Reference temperature equal to 25 degrees Celsius (K)"
    To::FT
    "Relative conductivity of water vapor compared to CO2 (unitless), equal to 1.6"
    Drel::FT
    "Constant factor appearing the dark respiration term for C3 plants, equal to 0.015."
    fC3::FT
    "Fitting constant to compute the moisture stress factor (Pa^{-1})"
    sc::FT
    "Fitting constant to compute the moisture stress factor (Pa)"
    pc::FT
    "Constant describing cost of maintaining electron transport (unitless)"
    c::FT
    "Ratio of unit costs of transpiration and carboxylation (unitless)"
    β::FT 
    "Scaling parameter for intrinsic quantum yield (unitless)"
    ϕc::FT
    "Use soil moisture stress in the photosynthesis model (default: true)"
    use_soil_moisture_stress::Bool = true
end

Base.eltype(::PModelParameters{FT}) where {FT} = FT

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
    update_photosynthesis!(p, Y, model::PModelModel,canopy)

Computes the net photosynthesis rate `An` (mol CO2/m^2/s)for the Optimality Farquhar model, along with the
dark respiration `Rd` (mol CO2/m^2/s), the value of `Vcmax25`(mol CO2/m^2/s) , and the gross primary 
productivity `GPP` (mol CO2/m^2/s), and updates them in place.
"""
function update_photosynthesis!(p, Y, model::PModelModel, canopy)
    # unpack a bunch of stuff from p and params
    Rd = p.canopy.photosynthesis.Rd
    An = p.canopy.photosynthesis.An
    GPP = p.canopy.photosynthesis.GPP
    Vcmax25 = p.canopy.photosynthesis.Vcmax25
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    f_abs = p.canopy.radiative_transfer.par.abs
    ψ = p.canopy.hydraulics.ψ
    ca = p.drivers.c_co2
    cosθs = p.drivers.cosθs
    P_air = p.drivers.P
    T_air = p.drivers.T
    q_air = p.drivers.q
    earth_param_set = canopy.parameters.earth_param_set
    lightspeed = LP.light_speed(earth_param_set)
    planck_h = LP.planck_constant(earth_param_set)
    N_a = LP.avogadro_constant(earth_param_set)
    grav = LP.grav(earth_param_set)
    ρ_l = LP.ρ_cloud_liq(earth_param_set)
    R = LP.gas_constant(earth_param_set)
    thermo_params = earth_param_set.thermo_params
    (; G_Function, λ_γ_PAR, Ω) = canopy.radiative_transfer.parameters
    energy_per_mole_photon_par = planck_h * lightspeed / λ_γ_PAR * N_a
    (; sc, pc) = canopy.photosynthesis.parameters
    # (; g1,) = canopy.conductance.parameters
    n_stem = canopy.hydraulics.n_stem
    n_leaf = canopy.hydraulics.n_leaf
    i_end = n_stem + n_leaf
    par_d = p.canopy.radiative_transfer.par_d
    area_index = p.canopy.hydraulics.area_index
    LAI = area_index.leaf

    (;
        ΔHVcmax,
        ΔHRd,
        fC3,
        To,
        Drel,
        is_c3,
        Kc25,
        Ko25,
        ΔHkc,
        ΔHko,
        c,
        β,
        ϕc,
        use_soil_moisture_stress # TODO: ideally this should take different functional form options 
    ) = model.parameters

    if use_soil_moisture_stress
        βm = @. lazy(quadratic_soil_moisture_stress())
    else
        βm = 1.0
    end

    # βm = @. lazy(moisture_stress(ψ.:($$i_end) * ρ_l * grav, sc, pc))

    # TODO: merge with co2_compensation used in optimality_farquhar.jl
    Γstar = @. lazy(co2_compensation_p(T_canopy, P_air, R))

    # amount of light absorbed [mol/m^2/s]
    I_abs = f_abs * par_d / energy_per_mole_photon_par

    # note that ϕc is called \hat{C}_L in Stocker et al. (2020) -- this is a calibratable parameter
    ϕ0 = @. lazy(intrinsic_quantum_yield(T_canopy, ϕc)) 
    ηstar = @. lazy(compute_viscosity_ratio(T_canopy, P_air))

    # compute effective Michaelis-Menten coefficient of Rubisco limited photosynthesis
    Kmm = @. lazy(compute_Kmm(T_canopy, P_air, Kc25, Ko25, ΔHkc, ΔHko, To, R)) 

    VPD = ClimaLand.vapor_pressure_deficit(T_air, P_air, q_air, thermo_params)
    χ, ξ, mj, mc = @. lazy(optimal_co2_ratio_c3(Kmm, Γstar, ηstar, ca, VPD, β, Drel))
    ci = χ * ca

    # compute Vcmax assuming optimality and coordination Aj = Ac 
    mprime = @. lazy(compute_mj_with_jmax_limitation(mj, 4.0 * c))
    Vcmax = @. lazy(pmodel_vcmax(ϕ0, I_abs, mprime, mc))

    # Equation 6 in Stocker et al. (2020)
    Ac = Vcmax * mc 

    Aj = Ac # coordination hypothesis 

    LUE =  @. lazy(compute_LUE(ϕ0, βm, mprime))

    @. Vcmax25 = Vcmax / arrhenius_function(T, To, R, ΔHVcmax)
    @. Rd = dark_respiration(is_c3, Vcmax25, βm, T, R, To, fC3, ΔHRd)
    @. An = net_photosynthesis(Ac, Aj, Rd, βm)
    # # Compute GPP: TODO - move to diagnostics only
    @. GPP = compute_GPP(An, extinction_coeff(G_Function, cosθs), LAI, Ω)
end

get_Vcmax25(p, m::PModelParameters) =
    p.canopy.photosynthesis.Vcmax25
Base.broadcastable(m::PModelParameters) = tuple(m)
