export PModelParameters, 
    PModelDrivers, 
    PModelConstants, 
    compute_pmodel_outputs, 
    set_historical_cache!,
    PModel,
    make_PModel_callback

"""
    PModelParameters{FT<:AbstractFloat}

The required parameters for P-model (Stocker et al. 2020). Parameters are typically
tunable with considerable uncertainty. 
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PModelParameters{
    FT <: AbstractFloat
}
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
    """Timescale parameter used in EMA for acclimation of optimal photosynthetic capacities (unitless)."""
    α::FT 

    sc::FT
    pc::FT 
end

"""
    PModelDrivers{FT<:AbstractFloat}

The required drivers for P-model (Stocker et al. 2020). Drivers are defined as
external variables used to compute the optimal photosynthetic capacities. 
$(DocStringExtensions.FIELDS)
"""
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

Base.broadcastable(x::PModelParameters) = tuple(x)
Base.broadcastable(x::PModelDrivers) = tuple(x)
Base.broadcastable(x::PModelConstants) = tuple(x)

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

"""
New cache variables:
- `OptVars`: a NamedTuple with keys `:ξ_opt`, `:Vcmax25_opt`, and `:Jmax25_opt` 
    containing the acclimated optimal values of ξ, Vcmax25, and Jmax25, respectively. These are updated
    using an exponential moving average (EMA) at local noon.
- `IntVars`: a NamedTuple with keys `:Γstar`, `:Kmm`, and `:ci` containing the common intermediate variables
    computed each timestep for instantaneous assimilation 
- `Jmax`: the instantaneous (temperature-adjusted) maximum electron transport rate (mol m^-2 s^-1)
- `J`: the instantaneous electron transport rate (mol m^-2 s^-1)
"""
ClimaLand.auxiliary_vars(model::PModel) =
    (:An, 
    :GPP, 
    :Rd, 
    :OptVars,
    :IntVars,
    :Jmax,
    :J)
ClimaLand.auxiliary_types(model::PModel{FT}) where {FT} =
    (FT,
    FT,
    FT,
    NamedTuple{(:ξ_opt, :Vcmax25_opt, :Jmax25_opt), Tuple{FT, FT, FT}},
    NamedTuple{(:Γstar, :Kmm, :ci), Tuple{FT, FT, FT}},
    FT,
    FT)
ClimaLand.auxiliary_domain_names(::PModel) =
    (:surface, :surface, :surface, :surface, :surface, :surface, :surface)


"""
    compute_pmodel_outputs(
        parameters::PModelParameters, 
        drivers::PModelDrivers, 
        constants::PModelConstants
    )

Performs the P-model computations as defined in Stocker et al. (2020) 
and returns a dictionary of outputs. See https://github.com/geco-bern/rpmodel
for a code reference.

Output name      Description (units)
    "gpp"           Gross primary productivity (kg m^-2 s^-1)
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
    "Rd"            Dark respiration rate (mol m^-2 s^-1)
"""
function compute_pmodel_outputs(
    parameters::PModelParameters{FT}, 
    drivers::PModelDrivers{FT}, 
    constants::PModelConstants{FT}
) where {FT}
    # Unpack parameters
    (; cstar, β, ϕc, ϕ0, ϕa0, ϕa1, ϕa2) = parameters

    # Unpack drivers
    (; T_canopy, I_abs, ca, P_air, VPD, βm) = drivers

    # Unpack constants
    (; R, Kc25, Ko25, To, ΔHkc, ΔHko, 
        Drel, ΔHΓstar, Γstar25,
        Ha_Vcmax, Hd_Vcmax, aS_Vcmax, bS_Vcmax, 
        Ha_Jmax, Hd_Jmax, aS_Jmax, bS_Jmax, Mc, oi, aRd, bRd, fC3) = constants

    # Compute intermediate values
    ϕ0 = isnan(ϕ0) ? intrinsic_quantum_yield(T_canopy, ϕc, ϕa0, ϕa1, ϕa2) : ϕ0

    Γstar = co2_compensation_p(T_canopy, To, P_air, R, ΔHΓstar, Γstar25)
    ηstar = compute_viscosity_ratio(T_canopy, P_air, true)
    Kmm = compute_Kmm(T_canopy, P_air, Kc25, Ko25, ΔHkc, ΔHko, To, R, oi)
    χ, ξ, mj, mc = optimal_co2_ratio_c3(Kmm, Γstar, ηstar, ca, VPD, β, Drel)
    ci = χ * ca
    mprime = compute_mj_with_jmax_limitation(mj, cstar)

    Vcmax = βm * ϕ0 * I_abs * mprime / mc
    inst_temp_scaling_vcmax25 = inst_temp_scaling(T_canopy, T_canopy, To, Ha_Vcmax, Hd_Vcmax, aS_Vcmax, bS_Vcmax, R)
    Vcmax25 = Vcmax / inst_temp_scaling_vcmax25

    Jmaxlim = Vcmax * (ci + FT(2) * Γstar) / (ϕ0 * I_abs * (ci + Kmm))
    Jmax = FT(4) * ϕ0 * I_abs / sqrt((FT(1)/Jmaxlim)^2 - FT(1)) 
    Jmax25 = Jmax / inst_temp_scaling(T_canopy, T_canopy, To, Ha_Jmax, Hd_Jmax, aS_Jmax, bS_Jmax, R)
    J = electron_transport_pmodel(ϕ0, I_abs, Jmax)

    Ac = Vcmax * mc
    Aj = J * mj / FT(4)

    LUE = compute_LUE(ϕ0, βm, mprime, Mc)
    GPP = I_abs * LUE 

    # intrinsic water use efficiency (iWUE) and stomatal conductance (gs)
    iWUE = (ca - ci) / Drel
    gs = pmodel_gs(χ, ca, Ac) 

    # dark respiration 
    rd = fC3 * (inst_temp_scaling_rd(T_canopy, To, aRd, bRd) / inst_temp_scaling_vcmax25) * Vcmax

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
        rd = rd
    )
end


function update_optimal_EMA(
    parameters::PModelParameters{FT}, 
    constants::PModelConstants{FT},
    OptVars::NamedTuple{(:ξ_opt, :Vcmax25_opt, :Jmax25_opt), Tuple{FT, FT, FT}}, 
    T_canopy::FT,
    P_air::FT,
    VPD::FT,
    ca::FT, 
    βm::FT,
    I_abs::FT,
    local_noon_mask::FT,
    verbose::Bool = false, 
) where {FT} 
    if local_noon_mask == FT(1.0)
        # Unpack parameters
        (; cstar, β, ϕc, ϕ0, ϕa0, ϕa1, ϕa2, α) = parameters

        # Unpack constants
        (; R, Kc25, Ko25, To, ΔHkc, ΔHko, 
            Drel, ΔHΓstar, Γstar25,
            Ha_Vcmax, Hd_Vcmax, aS_Vcmax, bS_Vcmax, 
            Ha_Jmax, Hd_Jmax, aS_Jmax, bS_Jmax, Mc, oi, aRd, bRd, fC3) = constants

        # Compute intermediate values
        ϕ0 = isnan(ϕ0) ? intrinsic_quantum_yield(T_canopy, ϕc, ϕa0, ϕa1, ϕa2) : ϕ0

        Γstar = co2_compensation_p(T_canopy, To, P_air, R, ΔHΓstar, Γstar25)
        ηstar = compute_viscosity_ratio(T_canopy, P_air, true)
        Kmm = compute_Kmm(T_canopy, P_air, Kc25, Ko25, ΔHkc, ΔHko, To, R, oi)
        
        verbose && print("Drivers: \n",
            "T_canopy: $T_canopy, P_air: $P_air, VPD: $VPD, ca: $ca, βm: $βm, I_abs: $I_abs\n\n"
        )
        # convert ca from mol/mol to Pa 
        ca_pp = ca * P_air 

        ξ = sqrt(β * (Kmm + Γstar) / (Drel * ηstar))
        χ = Γstar / ca_pp + (1 - Γstar / ca_pp) * ξ / (ξ + sqrt(VPD)) 
        γ = Γstar / ca_pp 
        κ = Kmm / ca_pp 
        mj = (χ - γ) / (χ + 2 * γ) # eqn 11 in Stocker et al. (2020)
        mc = (χ - γ) / (χ + κ) # eqn 7 in Stocker et al. (2020)
        mprime = compute_mj_with_jmax_limitation(mj, cstar)

        # print intermediate values 
        verbose && print("Intermediate values:\n",
            "Γstar: $Γstar, Kmm: $Kmm, ηstar: $ηstar, ca (Pa): $ca_pp\n",
            "ξ: $ξ, χ: $χ, γ: $γ, κ: $κ\n",
            "mprime: $mprime, mc: $mc, mj: $mj\n",
            "βm: $βm, ϕ0: $ϕ0, I_abs: $I_abs\n\n"
        )

        Vcmax = βm * ϕ0 * I_abs * mprime / mc
        Vcmax25 = Vcmax / inst_temp_scaling(T_canopy, T_canopy, To, Ha_Vcmax, Hd_Vcmax, aS_Vcmax, bS_Vcmax, R)

        Jmax = 4 * ϕ0 * I_abs / sqrt((mj / (βm * mprime)^2) - 1)
        Jmax25 = Jmax / inst_temp_scaling(T_canopy, T_canopy, To, Ha_Jmax, Hd_Jmax, aS_Jmax, bS_Jmax, R)
        verbose && print(
            "Vcmax: $Vcmax, Vcmax25: $Vcmax25, Jmax: $Jmax, Jmax25: $Jmax25\n
            ----------------------------------\n"
        )
        return (;
            ξ_opt = α * OptVars.ξ_opt + (1 - α) * ξ,
            Vcmax25_opt = α * OptVars.Vcmax25_opt + (1 - α) * Vcmax25,
            Jmax25_opt = α * OptVars.Jmax25_opt + (1 - α) * Jmax25
        )
    else 
        return OptVars
    end
end 


function update_intermediate_vars(
    constants::PModelConstants{FT},
    ξ_opt::FT,
    T_canopy::FT, 
    P_air::FT,
    VPD::FT, 
    ca::FT,
) where {FT}
    # Unpack constants
    (; R, Kc25, Ko25, To, ΔHkc, ΔHko, 
        Drel, ΔHΓstar, Γstar25,
        Ha_Vcmax, Hd_Vcmax, aS_Vcmax, bS_Vcmax, 
        Ha_Jmax, Hd_Jmax, aS_Jmax, bS_Jmax, Mc, oi, aRd, bRd, fC3) = constants

    # convert ca from ppm to Pa 
    ca_pp = ca * P_air 
    Γstar = co2_compensation_p(T_canopy, To, P_air, R, ΔHΓstar, Γstar25)
    ci = (ξ_opt * ca_pp + Γstar * sqrt(VPD)) / (ξ_opt + sqrt(VPD))

    return (;
        Γstar = Γstar,
        Kmm = compute_Kmm(T_canopy, P_air, Kc25, Ko25, ΔHkc, ΔHko, To, R, oi),
        ci = ci,
    )
end 

"""
function get_local_noon_mask(
    t::FT,
    dt::FT, 
    lon::FT 
) where {FT}    

This function determines whether the current time `t` (seconds UTC) is within a local noon window of width 
`dt` (seconds) centered around the local noon time for a given longitude `lon` (degrees, -180 to 180). 
Currently we neglect the correction due to obliquity and eccentricity. 

See https://clima.github.io/Insolation.jl/dev/ZenithAngleEquations/#Hour-Angle 
"""
function get_local_noon_mask(
    t,
    dt, 
    lon::FT
) where {FT}    
    day = IP.day(IP.InsolationParameters(FT)) 

    # local noon in seconds UTC 
    local_noon = day * (1/2 - lon/360)
    return FT(FT(t) >= local_noon - FT(dt) / 2 && FT(t) <= local_noon + FT(dt) / 2) 
end

"""
function set_historical_cache_pmodel!(p, Y0, model::PModel, canopy) 

The P-model requires a cache of optimal Vcmax25, Jmax25, and ξ that represent past acclimated values.
To set the initial condition, we assume that the acclimation is to the initial conditions of the simulation.
"""
function set_historical_cache!(p, Y0, model::PModel, canopy) 
    print("setting historical cache for p model\n")

    drivers = extract_pmodel_drivers(p, Y0, canopy)
    parameters = model.parameters
    constants = model.constants
    FT = eltype(parameters)

    # TODO: replace this with modular soil moisture stress parameterization
    βm = FT(1.0) #@. lazy(moisture_stress(drivers.ψ.:($$drivers.i_end) * drivers.ρ_l * drivers.grav, parameters.sc, parameters.pc)) 

    # Initialize OptVars with dummy values (will be overwritten since α=0)
    fill!(p.canopy.photosynthesis.OptVars, (;
        ξ_opt = FT(0),
        Vcmax25_opt = FT(0), 
        Jmax25_opt = FT(0)
    ))

    parameters_init = PModelParameters(
        cstar = parameters.cstar,
        β = parameters.β,
        ϕc = parameters.ϕc,
        ϕ0 = parameters.ϕ0,
        ϕa0 = parameters.ϕa0,
        ϕa1 = parameters.ϕa1,
        ϕa2 = parameters.ϕa2,
        α = FT(0),  # Override α to 0
        sc = parameters.sc,
        pc = parameters.pc
    )
    
    local_noon_mask = FT(1)  # Force update for initialization

    @. p.canopy.photosynthesis.OptVars = update_optimal_EMA(
        parameters_init, 
        constants, 
        p.canopy.photosynthesis.OptVars, 
        drivers.T_canopy, 
        drivers.P_air, 
        drivers.VPD,
        drivers.ca, 
        βm,
        drivers.I_abs,
        local_noon_mask,
    )
end


function call_update_optimal_EMA(p, Y, t; canopy, dt, longitude=nothing) 
    # update local noon mask 
    if isnothing(longitude) 
        try 
            longitude = ClimaCore.Fields.coordinate_field(axes(p.drivers.T)).long
        catch e
            error("Longitude must be provided explicitly if the domain you are working on does not \
                    have axes that specify longitude $e")
        end
    end 
    local_noon_mask = @. lazy(get_local_noon_mask(t, dt, longitude))
    
    # update the acclimated Vcmax25, Jmax25, ξ using EMA 
    drivers = extract_pmodel_drivers(p, Y, canopy)
    parameters = canopy.photosynthesis.parameters
    constants = canopy.photosynthesis.constants
    
    FT = eltype(parameters)
    # TODO: replace this with modular soil moisture stress parameterization
    βm = FT(1.0) # @. lazy(moisture_stress(ψ.:($$i_end) * ρ_l * grav, parameters.sc, parameters.pc)) 

    @. p.canopy.photosynthesis.OptVars = update_optimal_EMA(
        parameters, 
        constants, 
        p.canopy.photosynthesis.OptVars, 
        drivers.T_canopy, 
        drivers.P_air, 
        drivers.VPD,
        drivers.ca, 
        βm,
        drivers.I_abs,
        local_noon_mask,
    )
end 

"""
function make_PModel_callback(FT, start_date, dt, canopy, longitude=nothing) 

This constructs a FrequencyBasedCallback for the P-model that updates the optimal photosynthetic capacities
using the EMA equation at local noon every 10 minutes. 

Args
- `FT`: The floating-point type used in the model (e.g., `Float32`, `Float64`).
- `start_date`: datetime object for the start of the simulation (UTC).
- `dt`: timestep in seconds (this will get cast to type FT)
- `canopy`: the canopy object containing the P-model parameters and constants.
- `longitude`: optional longitude in degrees for local noon calculation (default is `nothing`). If we are on 
    a ClimaLand.Domains.Point, this will need to be supplied explicitly. Otherwise, if we are on a field, then 
    the longitude at each point can be automatically extracted from the field axes (THIS STILL NEEDS TO BE TESTED)
"""
function make_PModel_callback(FT, start_date, dt, canopy, longitude=nothing) 

    function seconds_after_midnight(t)
        return Hour(t).value * 3600 +
            Minute(t).value * 60
    end 

    day = IP.day(IP.InsolationParameters(FT)) 
    start_t = FT(seconds_after_midnight(start_date))

    return FrequencyBasedCallback(
        FT(600.0), # 10 minutes
        start_date, # start datetime, UTC 
        dt; # timestep, in seconds
        func = (integrator; ) -> call_update_optimal_EMA(
            integrator.p, 
            integrator.u,
            (integrator.t + start_t) % (day), # current time in seconds UTC; 
            canopy=canopy, dt=dt, longitude=longitude
        ),
    )
end


"""
    update_photosynthesis!(p, Y, model::PModel, canopy)

Computes the net photosynthesis rate `An` (mol CO2/m^2/s) for the P-model, along with the
dark respiration `Rd` (mol CO2/m^2/s), the value of `Vcmax25` (mol CO2/m^2/s), and the gross primary 
productivity `GPP` (mol CO2/m^2/s), and updates them in place.
"""
function update_photosynthesis!(p, Y, model::PModel, canopy)
    drivers = extract_pmodel_drivers(p, Y, canopy)
    parameters = model.parameters
    constants = model.constants

    # compute instantaneous max photosynthetic rates and assimilation rates 
    Vcmax_inst = @. lazy(
        p.canopy.photosynthesis.OptVars.Vcmax25_opt * inst_temp_scaling(
            drivers.T_canopy, 
            drivers.T_canopy, 
            constants.To, 
            constants.Ha_Vcmax, 
            constants.Hd_Vcmax, 
            constants.aS_Vcmax, 
            constants.bS_Vcmax, 
            constants.R
        )
    )

    Jmax_inst = @. lazy(
        p.canopy.photosynthesis.OptVars.Jmax25_opt * inst_temp_scaling(
            drivers.T_canopy, 
            drivers.T_canopy, 
            constants.To, 
            constants.Ha_Jmax, 
            constants.Hd_Jmax, 
            constants.aS_Jmax, 
            constants.bS_Jmax, 
            constants.R
        )
    )

    @. p.canopy.photosynthesis.IntVars = update_intermediate_vars(
        constants, 
        p.canopy.photosynthesis.OptVars.ξ_opt, 
        drivers.T_canopy, 
        drivers.P_air, 
        drivers.VPD, 
        drivers.ca
    )

    # Note: this is different than the Smith 2019 formulation used in optimality_farquhar.jl
    J_inst = @. lazy(
        electron_transport_pmodel(
            isnan(parameters.ϕ0) ? intrinsic_quantum_yield(
                drivers.T_canopy, 
                parameters.ϕc, 
                parameters.ϕa0, 
                parameters.ϕa1, 
                parameters.ϕa2) : parameters.ϕ0, 
            drivers.I_abs, 
            Jmax_inst
        )
    )

    # rubisco limited assimilation rate
    Ac = @. lazy(c3_rubisco_assimilation(
        Vcmax_inst,
        p.canopy.photosynthesis.IntVars.ci,
        p.canopy.photosynthesis.IntVars.Γstar,
        p.canopy.photosynthesis.IntVars.Kmm
    ))

    # light limited assimilation rate 
    Aj = @. lazy(c3_light_assimilation(
        J_inst, 
        p.canopy.photosynthesis.IntVars.ci,
        p.canopy.photosynthesis.IntVars.Γstar
    ))

    # dark respiration 
    Rd = @. lazy(
        constants.fC3 * (
            inst_temp_scaling_rd(
                drivers.T_canopy, 
                constants.To, 
                constants.aRd,
                constants.bRd) / 
            inst_temp_scaling(
                drivers.T_canopy, 
                drivers.T_canopy, 
                constants.To, 
                constants.Ha_Vcmax, 
                constants.Hd_Vcmax, 
                constants.aS_Vcmax, 
                constants.bS_Vcmax, 
                constants.R
            )) * Vcmax_inst
    )

    @. p.canopy.photosynthesis.Rd = Rd
    @. p.canopy.photosynthesis.An = min(Aj, Ac) - p.canopy.photosynthesis.Rd 
    @. p.canopy.photosynthesis.GPP = compute_GPP(p.canopy.photosynthesis.An, extinction_coeff(drivers.G_Function, p.drivers.cosθs), drivers.LAI, drivers.Ω)
    @. p.canopy.photosynthesis.Jmax = Jmax_inst
    @. p.canopy.photosynthesis.J = J_inst
end

get_Vcmax25(p, m::PModel) =
    p.canopy.photosynthesis.OptVars.Vcmax25_opt

"""
    extract_pmodel_drivers(p, Y, canopy)

Unified function to extract common drivers and parameters needed for P-model computations.
Returns a NamedTuple with all the commonly used variables across P-model functions.
"""
function extract_pmodel_drivers(p, Y, canopy)
    # Extract earth parameters
    earth_param_set = canopy.parameters.earth_param_set
    lightspeed = LP.light_speed(earth_param_set)
    planck_h = LP.planck_constant(earth_param_set)
    N_a = LP.avogadro_constant(earth_param_set)
    ρ_l = LP.ρ_cloud_liq(earth_param_set)
    grav = LP.grav(earth_param_set)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    
    # Extract radiative transfer parameters
    (; G_Function, λ_γ_PAR, Ω) = canopy.radiative_transfer.parameters
    energy_per_mole_photon_par = planck_h * lightspeed / λ_γ_PAR * N_a
    
    # Extract hydraulic parameters
    i_end = canopy.hydraulics.n_stem + canopy.hydraulics.n_leaf
    ψ = p.canopy.hydraulics.ψ
    
    # Compute derived variables
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    f_abs = p.canopy.radiative_transfer.par.abs
    ca = p.drivers.c_co2
    P_air = p.drivers.P
    par_d = p.canopy.radiative_transfer.par_d
    VPD = @. lazy(ClimaLand.vapor_pressure_deficit(
        p.drivers.T, p.drivers.P, p.drivers.q, thermo_params
    ))
    LAI = p.canopy.hydraulics.area_index.leaf
    I_abs = @. lazy(f_abs * par_d / energy_per_mole_photon_par)
    
    return (;
        earth_param_set,
        energy_per_mole_photon_par,
        i_end,
        ψ,
        ρ_l,
        grav,
        T_canopy,
        f_abs,
        ca,
        P_air,
        par_d,
        thermo_params,
        VPD,
        LAI,
        I_abs,
        G_Function,
        Ω
    )
end