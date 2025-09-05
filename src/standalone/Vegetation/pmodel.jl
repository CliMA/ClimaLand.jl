export PModelParameters,
    PModelConstants,
    PModel,
    compute_full_pmodel_outputs,
    set_historical_cache!,
    update_optimal_EMA,
    make_PModel_callback

"""
    PModelParameters{FT<:AbstractFloat}

The required parameters for P-model (Stocker et al. 2020). Parameters are typically
tunable with considerable uncertainty.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PModelParameters{FT <: AbstractFloat}
    "Constant describing cost of maintaining electron transport (unitless)"
    cstar::FT
    "Ratio of unit costs of transpiration and carboxylation (unitless)"
    β::FT
    "A boolean flag indicating if the quantum yield is a function of temperature or not"
    temperature_dep_yield::Bool
    "Temp-independent intrinsic quantum yield. (unitless); C3"
    ϕ0_c3::FT
    "Temp-independent intrinsic quantum yield. (unitless); C4"
    ϕ0_c4::FT
    """Constant term in temp-dependent intrinsic quantum yield (unitless); C3."""
    ϕa0_c3::FT
    """First order term in temp-dependent intrinsic quantum yield (K^-1); C3."""
    ϕa1_c3::FT
    """Second order term in temp-dependent intrinsic quantum yield (K^-2); C3."""
    ϕa2_c3::FT
    """Constant term in temp-dependent intrinsic quantum yield (unitless); C4."""
    ϕa0_c4::FT
    """First order term in temp-dependent intrinsic quantum yield (K^-1); C4."""
    ϕa1_c4::FT
    """Second order term in temp-dependent intrinsic quantum yield (K^-2); C4."""
    ϕa2_c4::FT
    """Timescale parameter used in EMA for acclimation of optimal photosynthetic capacities (unitless).
        Setting this to 0 represents no incorporation of past values. Since we update the EMA equation
        once per day, α = 1 - 1 day/τ where τ is the acclimation timescale in days."""
    α::FT
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
    """Michaelis-Menten parameter for carboxylation at 25°C (Pa)"""
    Kc25::FT
    """Michaelis-Menten parameter for oxygenation at 25°C (Pa)"""
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
    """Planck constant (J s)"""
    planck_h::FT
    """Speed of light (m s^-1)"""
    lightspeed::FT
    """Avogadro constant (mol^-1)"""
    N_a::FT
    """Density of water (kg m^-3)"""
    ρ_water::FT
end

Base.eltype(::PModelParameters{FT}) where {FT} = FT
Base.eltype(::PModelConstants{FT}) where {FT} = FT

# make these custom structs broadcastable as tuples
Base.broadcastable(x::PModelParameters) = tuple(x)
Base.broadcastable(x::PModelConstants) = tuple(x)

"""
    PModelConstants(FT)()

Creates a `PModelConstants` object with default values for the P-model constants.
See Stocker et al. (2020) Table A2 and references within for more information.
"""
function PModelConstants{FT}(;
    Kc25 = FT(39.97),
    Ko25 = FT(27480),
    ΔHkc = FT(79430),
    ΔHko = FT(36380),
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
    oi = FT(0.2095),
    aRd = FT(0.1012),
    bRd = FT(-0.0005),
    fC3 = FT(0.015),
) where {FT <: AbstractFloat}
    # Note: physical constants are not exposed to the user
    return PModelConstants{FT}(
        LP.get_default_parameter(FT, :universal_gas_constant),
        Kc25,
        Ko25,
        LP.get_default_parameter(FT, :kelvin_25C),
        ΔHkc,
        ΔHko,
        FT(1.6),
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
        FT(0.0120107),
        oi,
        aRd,
        bRd,
        fC3,
        LP.get_default_parameter(FT, :planck_constant),
        LP.get_default_parameter(FT, :light_speed),
        LP.get_default_parameter(FT, :avogadro_constant),
        LP.get_default_parameter(FT, :density_liquid_water),
    )
end

"""
    PModel{FT,
           OPFT <: PModelParameters{FT},
           OPCT <: PModelConstants{FT},
           F <: Union{FT, ClimaCore.Fields.Field},
           } <: AbstractPhotosynthesisModel{FT}

An implementation of the optimality photosynthesis model "P-model v1.0" of Stocker et al. (2020).

Stocker, B. D., Wang, H., Smith, N. G., Harrison, S. P., Keenan, T. F., Sandoval, D., Davis, T.,
    and Prentice, I. C.: P-model v1.0: an optimality-based light use efficiency model for simulating
    ecosystem gross primary production, Geosci. Model Dev., 13, 1545–1581,
    https://doi.org/10.5194/gmd-13-1545-2020, 2020.

The P-model computes photosynthesis rates at the canopy level, and ci, Γstar, Ko, Kc are in
units of Pa.
"""
struct PModel{
    FT,
    OPFT <: PModelParameters{FT},
    OPCT <: PModelConstants{FT},
    F <: Union{FT, ClimaCore.Fields.Field},
} <: AbstractPhotosynthesisModel{FT}
    "Required parameters for the P-model of Stocker et al. (2020)"
    parameters::OPFT
    "Constants for the P-model"
    constants::OPCT
    "Photosynthesis mechanism - 1 indicates C3, 0 indicates C4"
    is_c3::F
end

Base.eltype(::PModel{FT, OPFT, OPCT, F}) where {FT, OPFT, OPCT, F} = FT

"""
PModel{FT}(is_c3,
    parameters::PModelParameters{FT},
    constants::PModelConstants{FT} = PModelConstants(FT)
)

Outer constructor for the PModel struct. This takes a PModelParameters struct which includes
parameters with considerable uncertainty. PModelConstants is constructed by default to the
default values, but if you know what you are doing, you can override with your own constants.
"""
function PModel{FT}(
    is_c3,
    parameters::PModelParameters{FT};
    constants::PModelConstants{FT} = PModelConstants{FT}(),
) where {FT <: AbstractFloat}
    # if is_c3 is a field, is_c3 may contain values between 0.0 and 1.0 after regridding
    # this deals with that possibility by rounding to the closest int
    is_c3 = max.(min.(is_c3, FT(1)), FT(0)) # placeholder
    is_c3 = round.(is_c3)
    F = typeof(is_c3)
    return PModel{FT, typeof(parameters), typeof(constants), F}(
        parameters,
        constants,
        is_c3,
    )
end

"""
    ClimaLand.auxiliary_vars(model::PModel)
    ClimaLand.auxiliary_types(model::PModel)
    ClimaLand.auxiliary_domain_names(model::PModel)

Defines the auxiliary vars of the Pmode: canopy level net photosynthesis,
 canopy-level gross photosynthesis (`GPP`),
and dark respiration at the canopy level (`Rd`), and

- `OptVars`: a NamedTuple with keys `:ξ_opt`, `:Vcmax25_opt`, and `:Jmax25_opt`
    containing the acclimated optimal values of ξ, Vcmax25, and Jmax25, respectively. These are updated
    using an exponential moving average (EMA) at local noon.
"""
ClimaLand.auxiliary_vars(model::PModel) = (:An, :GPP, :Rd, :ci, :OptVars)
ClimaLand.auxiliary_types(model::PModel{FT}) where {FT} = (
    FT,
    FT,
    FT,
    FT,
    NamedTuple{(:ξ_opt, :Vcmax25_opt, :Jmax25_opt), Tuple{FT, FT, FT}},
)
ClimaLand.auxiliary_domain_names(::PModel) =
    (:surface, :surface, :surface, :surface, :surface)


"""
    compute_full_pmodel_outputs(
        parameters::PModelParameters{FT},
        constants::PModelConstants{FT},
        T_canopy::FT,
        P_air::FT,
        VPD::FT,
        ca::FT,
        βm::FT,
        APAR::FT;
        is_c3 = FT(1)
    ) where {FT}

Performs the P-model computations as defined in Stocker et al. (2020)
and returns a dictionary of full outputs. See https://github.com/geco-bern/rpmodel
for a code reference. This should replicate the behavior of the `rpmodel` package.

This is for C3 only, since the comparison data for C3.

Args:
- `parameters`:     PModelParameters object containing the model parameters.
- `constants`:      PModelConstants object containing the model constants.
- `T_canopy`:       Canopy temperature (K).
- `P_air`:          Ambient air pressure (Pa).
- `VPD`:            Vapor pressure deficit (Pa).
- `ca`:             Ambient CO2 concentration (mol/mol).
- `βm`:             Soil moisture stress factor (unitless).
- `APAR`:          Absorbed photosynthetically active radiation (mol photons m^-2 s^-1).

Returns: named tuple with the following keys and descriptions:
Output name         Description (units)
- "gpp"           Gross primary productivity (kg m^-2 s^-1)
- "gammastar"     CO2 compensation point (Pa)
- "kmm"           Effective MM coefficient for Rubisco-limited photosynthesis (Pa)
- "ns_star"       Viscosity of water normalized to 25 deg C (unitless)
- "chi"           Optimal ratio of intercellular to ambient CO2 (unitless)
- "xi"            Sensitivity of χ to VPD (Pa^1/2)
- "mj"            CO2 limitation factor for light-limited photosynthesis (unitless)
- "mc"            CO2 limitation factor for Rubisco-limited photosynthesis (unitless)
- "ci"            Intercellular CO2 concentration (Pa)
- "iwue"          Intrinsic water use efficiency (Pa)
- "gs"            Stomatal conductance (mol m^-2 s^-1 Pa^-1)
- "vcmax"         Maximum rate of carboxlation (mol m^-2 s^-1)
- "vcmax25"       Vcmax normalized to 25°C via modified-Arrhenius type function (mol m^-2 s^-1)
- "jmax"          Maximum rate of electron transport (mol m^-2 s^-1)
- "jmax25"        Jmax normalized to 25°C via modified-Arrhenius type function (mol m^-2 s^-1)
- "Rd"            Dark respiration rate (mol m^-2 s^-1)
"""
function compute_full_pmodel_outputs(
    parameters::PModelParameters{FT},
    constants::PModelConstants{FT},
    T_canopy::FT,
    P_air::FT,
    VPD::FT,
    ca::FT,
    βm::FT,
    APAR::FT;
    is_c3 = FT(1),
) where {FT}
    # Unpack parameters
    (; cstar, β) = parameters

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
        planck_h,
        lightspeed,
        N_a,
        ρ_water,
    ) = constants

    # Convert ca from mol/mol to a partial pressure (Pa)
    ca_pp = ca * P_air

    # Compute intermediate values
    ϕ0 = intrinsic_quantum_yield(is_c3, T_canopy, parameters)
    Γstar = co2_compensation_pmodel(T_canopy, To, P_air, R, ΔHΓstar, Γstar25)
    ηstar = compute_viscosity_ratio(T_canopy, To, ρ_water)
    Kmm = compute_Kmm(T_canopy, P_air, Kc25, Ko25, ΔHkc, ΔHko, To, R, oi)
    ξ = sqrt(β * (Kmm + Γstar) / (Drel * ηstar))
    ci = intercellular_co2_pmodel(ξ, ca_pp, Γstar, VPD)
    mj = compute_mj(is_c3, Γstar, ca_pp, ci, VPD)
    mc = compute_mc(is_c3, Γstar, ca_pp, ci, VPD, Kmm)
    mprime = compute_mj_with_jmax_limitation(mj, cstar)

    Vcmax = βm * ϕ0 * APAR * mprime / mc
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

    # check for negative arg before taking sqrt
    arg = (mj / (βm * mprime))^2 - 1
    Jmax = 4 * ϕ0 * APAR / (sqrt(max(arg, 0)) + eps(FT))
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
    J = electron_transport_pmodel(ϕ0, APAR, Jmax)

    Ac = Vcmax * mc
    Aj = J * mj / FT(4)

    LUE = compute_LUE(ϕ0, βm, mprime, Mc)
    GPP = APAR * LUE

    # intrinsic water use efficiency (iWUE) and stomatal conductance (gs)
    iWUE = (ca_pp - ci) / Drel
    χ = ci / ca_pp
    gs = gs_co2_pmodel(χ, ca, Ac)

    # dark respiration
    rd =
        fC3 *
        (
            inst_temp_scaling_rd(T_canopy, To, aRd, bRd) /
            inst_temp_scaling_vcmax25
        ) *
        Vcmax

    return (;
        gpp = GPP,
        gammastar = Γstar,
        kmm = Kmm,
        ca = ca_pp,
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


"""
    update_optimal_EMA(is_c3::FT,
        parameters::PModelParameters{FT},
        constants::PModelConstants{FT},
        OptVars::NamedTuple{(:ξ_opt, :Vcmax25_opt, :Jmax25_opt), Tuple{FT, FT, FT}},
        T_canopy::FT,
        P_air::FT,
        VPD::FT,
        ca::FT,
        βm::FT,
        APAR::FT,
        local_noon_mask::FT,
    ) where {FT}

This function updates the optimal photosynthetic capacities Vcmax25, Jmax25 and sensitivity of
stomatal conductance to dryness (ξ) using an exponential moving average (EMA) that computes new
optimal values at local noon, following Mengoli et al. (2022).

Args:
- `is_c3`:      Photosynthesis mechanism
- `parameters`: PModelParameters object containing the model parameters.
- `constants`: PModelConstants object containing the model constants.
- `OptVars`: NamedTuple containing the current optimal values of ξ, Vcmax25, and Jmax25.
- `T_canopy`: Canopy temperature (K).
- `P_air`: Ambient air pressure (Pa).
- `VPD`: Vapor pressure deficit (Pa).
- `ca`: Ambient CO2 concentration (mol/mol).
- `βm`: Soil moisture stress factor (unitless).
- `APAR`: Absorbed photosynthetically active radiation (mol photons m^-2 s^-1).
- `local_noon_mask`: A mask (0 or 1) indicating whether the current time is within the local noon window.

Returns:
- NamedTuple with updated optimal values of ξ, Vcmax25, and Jmax25

Reference:
Mengoli, G., Agustí-Panareda, A., Boussetta, S., Harrison, S. P., Trotta, C., & Prentice, I. C. (2022).
Ecosystem photosynthesis in land-surface models: A first-principles approach incorporating acclimation.
Journal of Advances in Modeling Earth Systems, 14, e2021MS002767. https://doi.org/10.1029/2021MS002767
"""
function update_optimal_EMA(
    is_c3::FT,
    parameters::PModelParameters{FT},
    constants::PModelConstants{FT},
    OptVars::NamedTuple{(:ξ_opt, :Vcmax25_opt, :Jmax25_opt), Tuple{FT, FT, FT}},
    T_canopy::FT,
    P_air::FT,
    VPD::FT,
    ca::FT,
    βm::FT,
    APAR_canopy_moles::FT,
    local_noon_mask::FT,
) where {FT}
    if local_noon_mask == FT(1.0)
        # Unpack parameters
        (; cstar, β, α) = parameters

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
            planck_h,
            lightspeed,
            N_a,
            ρ_water,
        ) = constants

        # Compute intermediate values
        ϕ0 = intrinsic_quantum_yield(is_c3, T_canopy, parameters)

        Γstar =
            co2_compensation_pmodel(T_canopy, To, P_air, R, ΔHΓstar, Γstar25)
        ηstar = compute_viscosity_ratio(T_canopy, To, ρ_water)
        Kmm = compute_Kmm(T_canopy, P_air, Kc25, Ko25, ΔHkc, ΔHko, To, R, oi)

        # convert ca from mol/mol to Pa
        ca_pp = ca * P_air

        ξ = sqrt(β * (Kmm + Γstar) / (Drel * ηstar))
        ci = intercellular_co2_pmodel(ξ, ca_pp, Γstar, VPD)
        mj = compute_mj(is_c3, Γstar, ca_pp, ci, VPD)
        mc = compute_mc(is_c3, Γstar, ca_pp, ci, VPD, Kmm)
        mprime = compute_mj_with_jmax_limitation(mj, cstar)

        Vcmax = βm * ϕ0 * APAR_canopy_moles * mprime / mc
        Vcmax25 =
            Vcmax / inst_temp_scaling(
                T_canopy,
                T_canopy,
                To,
                Ha_Vcmax,
                Hd_Vcmax,
                aS_Vcmax,
                bS_Vcmax,
                R,
            )

        Jmax = 4 * ϕ0 * APAR_canopy_moles / sqrt((mj / (βm * mprime))^2 - 1)
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
        return (;
            ξ_opt = α * OptVars.ξ_opt + (1 - α) * ξ,
            Vcmax25_opt = α * OptVars.Vcmax25_opt + (1 - α) * Vcmax25,
            Jmax25_opt = α * OptVars.Jmax25_opt + (1 - α) * Jmax25,
        )
    else
        return OptVars
    end
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
function get_local_noon_mask(t, dt, local_noon::FT) where {FT}
    strict_noon_mask =
        FT(t) >= local_noon - FT(dt) / 2 && FT(t) <= local_noon + FT(dt) / 2
    return FT(strict_noon_mask)
end

"""
    function set_historical_cache!(p, Y0, model::PModel, canopy)

The P-model requires a cache of optimal Vcmax25, Jmax25, and ξ that represent past acclimated values.
Before the simulation, we need to have some physically meaningful initial values for these variables,
which live in p.canopy.photosynthesis.OptVars.

This method assumes that the acclimation is to the initial conditions of the simulation. Note that
if the initial condition is e.g., nighttime, then initially the optimal Vcmax and Jmax are
zero, so it will take ~1 month (two e-folding timescales for α corresponding to 2 week acclimation
timescale) for the model to reach a physically meaningful state.

An alternative to this approach is to initialize the initial optimal values to some reasonable values
based on a spun-up simulation.
"""
function set_historical_cache!(p, Y0, model::PModel, canopy)
    parameters = model.parameters
    constants = model.constants

    # drivers
    FT = eltype(parameters)
    earth_param_set = canopy.parameters.earth_param_set

    ψ = p.canopy.hydraulics.ψ
    n = canopy.hydraulics.n_leaf + canopy.hydraulics.n_stem
    grav = LP.grav(earth_param_set)
    ρ_water = LP.ρ_cloud_liq(earth_param_set)
    βm = p.canopy.soil_moisture_stress.βm
    T_canopy = canopy_temperature(canopy.energy, canopy, Y0, p)
    # The Pmodel divides by sqrt(VPD); clip here to prevent numerical issues
    VPD = @. lazy(
        max(
            ClimaLand.vapor_pressure_deficit(
                p.drivers.T,
                p.drivers.P,
                p.drivers.q,
                LP.thermodynamic_parameters(canopy.parameters.earth_param_set),
            ),
            sqrt(eps(FT)),
        ),
    )
    APAR_canopy_moles = @. lazy(
        compute_APAR_canopy_moles(
            p.canopy.radiative_transfer.par.abs,
            p.canopy.radiative_transfer.par_d,
            canopy.radiative_transfer.parameters.λ_γ_PAR,
            constants.lightspeed,
            constants.planck_h,
            constants.N_a,
        ),
    )

    # Initialize OptVars with dummy values which will be overwritten
    fill!(
        p.canopy.photosynthesis.OptVars,
        (; ξ_opt = FT(0), Vcmax25_opt = FT(0), Jmax25_opt = FT(0)),
    )

    parameters_init = PModelParameters(
        cstar = parameters.cstar,
        β = parameters.β,
        ϕ0_c4 = parameters.ϕ0_c4,
        ϕ0_c3 = parameters.ϕ0_c3,
        ϕa0_c3 = parameters.ϕa0_c3,
        ϕa1_c3 = parameters.ϕa1_c3,
        ϕa2_c3 = parameters.ϕa2_c3,
        ϕa0_c4 = parameters.ϕa0_c4,
        ϕa1_c4 = parameters.ϕa1_c4,
        ϕa2_c4 = parameters.ϕa2_c4,
        α = FT(0),  # this allows us to use the initial values directly
    )

    local_noon_mask = FT(1)  # Force update for initialization

    @. p.canopy.photosynthesis.OptVars = update_optimal_EMA(
        model.is_c3,
        parameters_init,
        constants,
        p.canopy.photosynthesis.OptVars,
        T_canopy,
        p.drivers.P,
        VPD,
        p.drivers.c_co2,
        βm,
        APAR_canopy_moles,
        local_noon_mask,
    )
end

"""
    call_update_optimal_EMA(
        p,
        Y,
        t;
        canopy,
        dt,
        local_noon,
    )

Updates the optimal parameters according to conditions at local noon.
"""
function call_update_optimal_EMA(p, Y, t; canopy, dt, local_noon)
    local_noon_mask = @. lazy(get_local_noon_mask(t, dt, local_noon))

    # update the acclimated Vcmax25, Jmax25, ξ using EMA
    parameters = canopy.photosynthesis.parameters
    constants = canopy.photosynthesis.constants
    earth_param_set = canopy.parameters.earth_param_set

    # drivers
    FT = eltype(parameters)
    βm = p.canopy.soil_moisture_stress.βm
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    # The Pmodel divides by sqrt(VPD); clip here to prevent numerical issues
    VPD = @. lazy(
        max(
            ClimaLand.vapor_pressure_deficit(
                p.drivers.T,
                p.drivers.P,
                p.drivers.q,
                LP.thermodynamic_parameters(earth_param_set),
            ),
            sqrt(eps(FT)),
        ),
    )
    APAR_canopy_moles = @. lazy(
        compute_APAR_canopy_moles(
            p.canopy.radiative_transfer.par.abs,
            p.canopy.radiative_transfer.par_d,
            canopy.radiative_transfer.parameters.λ_γ_PAR,
            constants.lightspeed,
            constants.planck_h,
            constants.N_a,
        ),
    )

    @. p.canopy.photosynthesis.OptVars = update_optimal_EMA(
        canopy.photosynthesis.is_c3,
        parameters,
        constants,
        p.canopy.photosynthesis.OptVars,
        T_canopy,
        p.drivers.P,
        VPD,
        p.drivers.c_co2,
        βm,
        APAR_canopy_moles,
        local_noon_mask,
    )
end

"""
    function make_PModel_callback(
        ::Type{FT},
        start_date::Dates.DateTime,
        dt::Union{AbstractFloat, Dates.Period},
        canopy,
        longitude = nothing
    ) where {FT <: AbstractFloat}

This constructs a FrequencyBasedCallback for the P-model that updates the optimal photosynthetic capacities
using an exponential moving average at local noon.

We check for local noon using the provided `longitude` (once passing
in lat/lon for point/column domains, this can be automatically extracted from the domain axes) every dt.
The time of local noon is expressed in seconds UTC and neglects the effects of obliquity and eccentricity, so
it is constant throughout the year. This presents an error of up to ~20 minutes, but it is sufficient for our
application here (since meteorological drivers are often updated at coarser time intervals anyway).

Args
- `FT`: The floating-point type used in the model (e.g., `Float32`, `Float64`).
- `start_date`: datetime object for the start of the simulation (UTC).
- `dt`: timestep
- `canopy`: the canopy object containing the P-model parameters and constants.
- `longitude`: optional longitude in degrees for local noon calculation (default is `nothing`, which means
    that it will be inferred from the canopy domain).
"""
function make_PModel_callback(
    ::Type{FT},
    start_date::Dates.DateTime,
    dt::Union{AbstractFloat, Dates.Period},
    canopy,
    longitude = nothing,
) where {FT <: AbstractFloat}
    function seconds_after_midnight(date)
        return Float64(
            Hour(date).value * 3600 +
            Minute(date).value * 60 +
            Second(date).value,
        )
    end

    if isnothing(longitude)
        try
            longitude = get_long(canopy.domain.space.surface)
        catch e
            error(
                "Longitude must be provided explicitly if the domain you are working on does not \
                  have axes that specify longitude $e",
            )
        end
    end

    # this computes the time of local noon in seconds UTC without considering the
    # effects of obliquity and orbital eccentricity, so it is constant throughout the year
    # the max error is on the order of 20 minutes
    seconds_in_a_day = IP.day(IP.InsolationParameters(FT))
    start_t = seconds_after_midnight(start_date)
    local_noon = @. seconds_in_a_day * (FT(1 / 2) - longitude / 360) # allocates, but only on init

    return FrequencyBasedCallback(
        dt,         # period of this callback
        start_date, # initial datetime, UTC
        dt;         # integration timestep, in seconds
        func = (integrator) -> call_update_optimal_EMA(
            integrator.p,
            integrator.u,
            (float(integrator.t) + start_t) % (seconds_in_a_day), # current time in seconds UTC;
            canopy = canopy,
            dt = dt,
            local_noon = local_noon,
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
    parameters = model.parameters
    constants = model.constants
    FT = eltype(parameters)

    # Unpack preallocated variables to short names
    An = p.canopy.photosynthesis.An
    GPP = p.canopy.photosynthesis.GPP
    Rd = p.canopy.photosynthesis.Rd
    ci = p.canopy.photosynthesis.ci
    OptVars = p.canopy.photosynthesis.OptVars

    # drivers
    P_air = p.drivers.P
    ca_pp = @. lazy(p.drivers.c_co2 * P_air) # partial pressure of co2
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    # The Pmodel divides by sqrt(VPD); clip here to prevent numerical issues
    VPD = @. lazy(
        max(
            ClimaLand.vapor_pressure_deficit(
                p.drivers.T,
                p.drivers.P,
                p.drivers.q,
                LP.thermodynamic_parameters(canopy.parameters.earth_param_set),
            ),
            sqrt(eps(FT)),
        ),
    )
    APAR_canopy_moles = @. lazy(
        compute_APAR_canopy_moles(
            p.canopy.radiative_transfer.par.abs,
            p.canopy.radiative_transfer.par_d,
            canopy.radiative_transfer.parameters.λ_γ_PAR,
            constants.lightspeed,
            constants.planck_h,
            constants.N_a,
        ),
    )

    # compute intermediate vars
    Γstar = @. lazy(
        co2_compensation_pmodel(
            T_canopy,
            constants.To,
            P_air,
            constants.R,
            constants.ΔHΓstar,
            constants.Γstar25,
        ),
    )
    @. ci = intercellular_co2_pmodel(OptVars.ξ_opt, ca_pp, Γstar, VPD)
    Kmm = @. lazy(
        compute_Kmm(
            T_canopy,
            P_air,
            constants.Kc25,
            constants.Ko25,
            constants.ΔHkc,
            constants.ΔHko,
            constants.To,
            constants.R,
            constants.oi,
        ),
    )
    # compute instantaneous max photosynthetic rates and assimilation rates
    Jmax = @. lazy(
        p.canopy.photosynthesis.OptVars.Jmax25_opt * inst_temp_scaling(
            T_canopy,
            T_canopy,
            constants.To,
            constants.Ha_Jmax,
            constants.Hd_Jmax,
            constants.aS_Jmax,
            constants.bS_Jmax,
            constants.R,
        ),
    )

    J = @. lazy(
        electron_transport_pmodel(
            intrinsic_quantum_yield(model.is_c3, T_canopy, parameters),
            APAR_canopy_moles,
            Jmax,
        ),
    )

    Vcmax = @. lazy(
        p.canopy.photosynthesis.OptVars.Vcmax25_opt * inst_temp_scaling(
            T_canopy,
            T_canopy,
            constants.To,
            constants.Ha_Vcmax,
            constants.Hd_Vcmax,
            constants.aS_Vcmax,
            constants.bS_Vcmax,
            constants.R,
        ),
    )

    # rubisco limited assimilation rate
    Ac = @. lazy(Vcmax * compute_mc(model.is_c3, Γstar, ca_pp, ci, VPD, Kmm)) # c3 or c4 is reflected in the value of mc

    # light limited assimilation rate
    Aj = @. lazy(J / 4 * compute_mj(model.is_c3, Γstar, ca_pp, ci, VPD)) # c3 or c4 is reflected in the value of mj

    # dark respiration
    # Here we make an assumption about how to relate Rd25 to Vcmax25_opt
    # To extend to C4, defined `compute_dark_respiration_pmodel() which dispatches off of the is_c3 field
    # This function below would become c3_dark_respiration_pmodel
    @. Rd =
        constants.fC3 *
        p.canopy.photosynthesis.OptVars.Vcmax25_opt *
        inst_temp_scaling_rd(
            T_canopy,
            constants.To,
            constants.aRd,
            constants.bRd,
        )

    # Note: net_photosynthesis applies the moisture stress to GPP, but since the P-model already applies
    # this factor to Vcmax and Jmax, we do not apply it again here
    @. GPP = gross_photosynthesis(Ac, Aj)
    @. An = net_photosynthesis(GPP, Rd)

end

get_Vcmax25_leaf(p, m::PModel) = @. lazy(
    p.canopy.photosynthesis.OptVars.Vcmax25_opt / max(
        p.canopy.hydraulics.area_index.leaf,
        sqrt(eps(eltype(m.constants))),
    ),
)
get_Rd_leaf(p, m::PModel) = @. lazy(
    p.canopy.photosynthesis.Rd / max(
        p.canopy.hydraulics.area_index.leaf,
        sqrt(eps(eltype(m.constants))),
    ),
)
get_An_leaf(p, m::PModel) = @.lazy(
    p.canopy.photosynthesis.An / max(
        p.canopy.hydraulics.area_index.leaf,
        sqrt(eps(eltype(m.constants))),
    ),
)

get_GPP_canopy(p, m::PModel) = p.canopy.photosynthesis.GPP

function get_J_over_Jmax(Y, p, canopy, m::PModel)
    Jmax = compute_Jmax_canopy(Y, p, canopy, m) # lazy
    J = compute_J_canopy(Y, p, canopy, m) # lazy
    FT = eltype(m.constants)
    return @. lazy(J / max(Jmax, sqrt(eps(FT))))
end

function compute_Jmax_canopy(Y, p, canopy, m::PModel) # used internally to pmodel photosynthesis as a helper function
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    constants = m.constants
    return @. lazy(
        p.canopy.photosynthesis.OptVars.Jmax25_opt * inst_temp_scaling(
            T_canopy,
            T_canopy,
            constants.To,
            constants.Ha_Jmax,
            constants.Hd_Jmax,
            constants.aS_Jmax,
            constants.bS_Jmax,
            constants.R,
        ),
    )
end

function compute_J_canopy(Y, p, canopy, m::PModel) # used internally to pmodel photosynthesis as a helper function
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    earth_param_set = canopy.parameters.earth_param_set
    f_abs_par = p.canopy.radiative_transfer.par.abs
    par_d = p.canopy.radiative_transfer.par_d
    (; λ_γ_PAR,) = canopy.radiative_transfer.parameters
    c = LP.light_speed(earth_param_set)
    planck_h = LP.planck_constant(earth_param_set)
    N_a = LP.avogadro_constant(earth_param_set)
    APAR_canopy_moles = @. lazy(
        compute_APAR_canopy_moles(f_abs_par, par_d, λ_γ_PAR, c, planck_h, N_a),
    )

    Jmax_canopy = compute_Jmax_canopy(Y, p, canopy, m)
    parameters = m.parameters
    constants = m.constants
    return @. lazy(
        electron_transport_pmodel(
            intrinsic_quantum_yield(m.is_c3, T_canopy, parameters),
            APAR_canopy_moles,
            Jmax_canopy,
        ),
    )
end

"""
    compute_Kmm(
        T::FT,
        p::FT,
        Kc25::FT,
        Ko25::FT,
        ΔHkc::FT,
        ΔHko::FT,
        To::FT,
        R::FT,
        oi::FT
    ) where {FT}

Computes the effective Michaelis-Menten coefficient for Rubisco-limited photosynthesis (`Kmm`),
in units Pa, as a function of temperature T (K), atmospheric pressure p (Pa), and constants:
Kc25 (Michaelis-Menten coefficient for CO2 at 25 °C), Ko25 (Michaelis-Menten coefficient for O2 at 25 °C),
ΔHkc (effective enthalpy of activation for Kc), ΔHko (effective enthalpy of activation for Ko),
To (reference temperature, typically 298.15 K), R (universal gas constant), and oi (O2 mixing ratio,
typically 0.209).
"""
function compute_Kmm(
    T::FT,
    p::FT,
    Kc25::FT,
    Ko25::FT,
    ΔHkc::FT,
    ΔHko::FT,
    To::FT,
    R::FT,
    oi::FT,
) where {FT}
    Kc = MM_Kc(Kc25, ΔHkc, T, To, R)
    Ko = MM_Ko(Ko25, ΔHko, T, To, R)

    return Kc * (1 + po2(p, oi) / Ko)
end

"""
    intrinsic_quantum_yield(
        is_c3::FT, T::FT, parameters) where {FT}
Computes the intrinsic quantum yield of photosystem II.
"""
function intrinsic_quantum_yield(is_c3::FT, T::FT, parameters) where {FT}
    if is_c3 > 0.5
        parameters.temperature_dep_yield ?
        quadratic_intrinsic_quantum_yield(
            T,
            parameters.ϕa0_c3,
            parameters.ϕa1_c3,
            parameters.ϕa2_c3,
        ) : parameters.ϕ0_c3
    else
        parameters.temperature_dep_yield ?
        quadratic_intrinsic_quantum_yield(
            T,
            parameters.ϕa0_c4,
            parameters.ϕa1_c4,
            parameters.ϕa2_c4,
        ) : parameters.ϕ0_c4
    end
end

"""
    quadratic_intrinsic_quantum_yield(
        T::FT,
        ϕa0::FT,
        ϕa1::FT,
        ϕa2::FT
    ) where {FT}

Computes the intrinsic quantum yield of photosynthesis ϕ (mol/mol)
as a function of temperature T (K); appropriate for C3 or C4 plants depending on the values of ϕa0 (unitless), ϕa1 (1/degrees C), ϕa2 (1/degrees C^2). 

The functional form given in Bernacchi et al (2003) and used in Stocker
et al. (2020) is a second order polynomial in T (deg C) with coefficients ϕa0,
ϕa1, and ϕa2.
"""
function quadratic_intrinsic_quantum_yield(
    T::FT,
    ϕa0::FT,
    ϕa1::FT,
    ϕa2::FT,
) where {FT}
    # convert to C
    T = T - FT(273.15)
    ϕ = ϕa0 + ϕa1 * T + ϕa2 * T^2
    return min(max(ϕ, FT(0)), FT(1)) # Clip to [0,1]
end


"""
    viscosity_h2o(
        T::FT,
        ρ_water::FT
    ) where {FT}

Computes the viscosity of water in Pa s given temperature T (K) and density ρ_water (kg/m^3)
according to Huber et al. (2009) [https://doi.org/10.1063/1.3088050].

Can consider simplifying if this level of precision is not needed
"""
function viscosity_h2o(T::FT, ρ::FT) where {FT <: AbstractFloat}
    # Reference constants
    tk_ast = FT(647.096)     # K
    ρ_ast = FT(322.0)       # kg m⁻³
    μ_ast = FT(1e-6)        # Pa s

    # Dimensionless variables
    tbar = T / tk_ast
    tbarx = sqrt(tbar)
    tbar2 = tbar * tbar
    tbar3 = tbar2 * tbar
    ρbar = ρ / ρ_ast

    # μ0 (Eq. 11 & Table 2)
    μ0 = (
        FT(1.67752) + FT(2.20462) / tbar + FT(0.6366564) / tbar2 -
        FT(0.241605) / tbar3
    )
    μ0 = FT(100) * tbarx / μ0

    # Using tuples keeps everything isbits and avoids heap allocation on device
    # blame the formatter
    h = (
        (
            FT(0.520094),
            FT(0.0850895),
            FT(-1.08374),
            FT(-0.289555),
            FT(0.0),
            FT(0.0),
        ),
        (
            FT(0.222531),
            FT(0.999115),
            FT(1.88797),
            FT(1.26613),
            FT(0.0),
            FT(0.120573),
        ),
        (
            FT(-0.281378),
            FT(-0.906851),
            FT(-0.772479),
            FT(-0.489837),
            FT(-0.257040),
            FT(0.0),
        ),
        (FT(0.161913), FT(0.257399), FT(0.0), FT(0.0), FT(0.0), FT(0.0)),
        (FT(-0.0325372), FT(0.0), FT(0.0), FT(0.0698452), FT(0.0), FT(0.0)),
        (FT(0.0), FT(0.0), FT(0.0), FT(0.0), FT(0.00872102), FT(0.0)),
        (FT(0.0), FT(0.0), FT(0.0), FT(-0.00435673), FT(0.0), FT(-0.000593264)),
    )

    # μ1 (Eq. 12 & Table 3), with Horner evaluation and iterative powers
    ctbar = inv(tbar) - one(FT)
    δρ = ρbar - one(FT)
    μ1 = zero(FT)
    coef1 = one(FT)
    @inbounds for i in 1:6
        coef2 = h[7][i]
        coef2 = muladd(δρ, coef2, h[6][i])
        coef2 = muladd(δρ, coef2, h[5][i])
        coef2 = muladd(δρ, coef2, h[4][i])
        coef2 = muladd(δρ, coef2, h[3][i])
        coef2 = muladd(δρ, coef2, h[2][i])
        coef2 = muladd(δρ, coef2, h[1][i])

        μ1 = muladd(coef1, coef2, μ1)  # accumulate ctbar^(i-1) * coef2
        coef1 *= ctbar # update ctbar power for next i
    end
    μ1 = exp(ρbar * μ1)

    μ_bar = μ0 * μ1
    μ = μ_bar * μ_ast
    return μ
end


"""
    compute_viscosity_ratio(
        T::FT,
        To::FT,
        ρ_water::FT
    ) where {FT}

Computes η*, the ratio of the viscosity of water at temperature T to that at To = 25˚C.
"""
function compute_viscosity_ratio(T::FT, To::FT, ρ_water::FT) where {FT}
    η25 = viscosity_h2o(To, ρ_water)
    ηstar = viscosity_h2o(T, ρ_water) / η25
    return FT(ηstar)
end


"""
    po2(
        P_air::FT,
        oi::FT
    ) where {FT}

Computes the partial pressure of O2 in the air (Pa) given atmospheric pressure (`P_air`)
and a constant mixing ratio of O2 (`oi`), typically 0.209.
"""
function po2(P_air::FT, oi::FT) where {FT}
    return oi * P_air
end

"""
    co2_compensation_pmode(
        T::FT,
        To::FT,
        p::FT,
        R::FT,
        ΔHΓstar::FT
        Γstar25::FT
    ) where {FT}

Computes the CO2 compensation point (`Γstar`), in units Pa, as a function of temperature T (K)
and pressure p (Pa). See Equation B5 of Stocker et al. (2020).
"""
function co2_compensation_pmodel(
    T::FT,
    To::FT,
    p::FT,
    R::FT,
    ΔHΓstar::FT,
    Γstar25::FT,
) where {FT}
    Γstar = Γstar25 * p / FT(101325.0) * arrhenius_function(T, To, R, ΔHΓstar)
    return Γstar
end

"""
    intercellular_co2_pmodel(ξ::FT, ca_pp::FT, Γstar::FT, VPD::FT) where {FT}

Computes the intercellular co2 concentration (`ci`) as a function of the
optimal `ξ` (sensitivity to dryness), `ca_pp` (ambient CO2 partial pressure),
`Γstar` (CO2 compensation point), and `VPD` (vapor pressure deficit).
"""
function intercellular_co2_pmodel(
    ξ::FT,
    ca_pp::FT,
    Γstar::FT,
    VPD::FT,
) where {FT}
    # VPD has been regularized already (VPD >= eps)
    return (ξ * ca_pp + Γstar * sqrt(VPD)) / (ξ + sqrt(VPD))
end



"""
    gs_co2_pmodel(
        χ::FT,
        ca::FT,
        A::FT
    ) where {FT}

Computes the stomatal conductance of CO2 (`gs_co2`), in units of mol CO2/m^2/s
via Fick's law. Parameters are the ratio of intercellular to ambient CO2
concentration (`χ`), the ambient CO2 concentration (`ca`, in mol/mol), and the
assimilation rate (`A`, mol m^-2 s^-1). This is related to the conductance of water by a
factor Drel (default value = 1.6).
"""
function gs_co2_pmodel(χ::FT, ca::FT, A::FT) where {FT}
    return A / (ca * (1 - χ) + eps(FT))
end

"""
    gs_h2o_pmodel(
        χ::FT,
        ca::FT,
        A::FT,
        Drel::FT
    ) where {FT}

Computes the stomatal conductance of H2O (`gs_h2o`), in units of mol H2O/m^2/s
via Fick's law. Parameters are the ratio of intercellular to ambient CO2
concentration (`χ`), the ambient CO2 concentration (`ca`, in mol/mol), the
assimilation rate (`A`, mol m^-2 s^-1), and the relative conductivity ratio `Drel` (unitless).
"""
function gs_h2o_pmodel(χ::FT, ca::FT, A::FT, Drel::FT) where {FT}
    return Drel * gs_co2_pmodel(χ, ca, A)
end

"""
    compute_mj_with_jmax_limitation(
        mj::FT,
        cstar::FT
    ) where {FT}

Computes m' such that Aj = ϕ0 APAR * m' (a LUE model) by assuming that dA/dJmax = c
is constant. cstar is defined as 4c, a free parameter. Wang etal (2017) derive cstar = 0.412
at STP and using Vcmax/Jmax = 1.88.
"""
function compute_mj_with_jmax_limitation(mj::FT, cstar::FT) where {FT}
    arg = cstar / mj
    arg = arg < 0 ? FT(0) : arg
    arg = 1 - arg^(FT(2 / 3))
    sqrt_arg = arg < 0 ? FT(0) : sqrt(arg) # avoid complex numbers
    return FT(mj * sqrt_arg)
end


"""
    compute_LUE(
        ϕ0::FT,
        β::FT,
        mprime::FT,
        Mc::FT
    ) where {FT}

Computes light use efficiency (LUE) in kg C/mol from intrinsic quantum yield (`ϕ0`),
moisture stress factor (`β`), and a Jmax modified capacity (`mprime`); see Eqn 17 and 19
in Stocker et al. (2020). Mc is the molar mass of carbon (kg/mol) = 0.0120107 kg/mol.
"""
function compute_LUE(ϕ0::FT, β::FT, mprime::FT, Mc::FT) where {FT}
    return ϕ0 * β * mprime * Mc
end


"""
    vcmax_pmodel(
        ϕ0::FT,
        APAR::FT,
        mprime::FT,
        mc::FT
        βm::FT
    ) where {FT}

Computes the maximum rate of carboxylation assuming optimality and Aj = Ac using
the intrinsic quantum yield (`ϕ0`), absorbed photosynthetically active radiation (`APAR`),
Jmax-adjusted capacity (`mprime`), a Rubisco-limited capacity (`mc`), and empirical
soil moisture stress factor (`βm`). See Eqns 16 and 6 in Stocker et al. (2020).
"""
function vcmax_pmodel(ϕ0::FT, APAR::FT, mprime::FT, mc::FT, βm::FT) where {FT}
    Vcmax = βm * ϕ0 * APAR * mprime / mc
    return Vcmax
end


"""
    electron_transport_pmodel(
        ϕ0::FT,
        APAR::FT,
        Jmax::FT
    ) where {FT}

Computes the rate of electron transport (`J`) in mol electrons/m^2/s for the pmodel.
"""
function electron_transport_pmodel(ϕ0::FT, APAR::FT, Jmax::FT) where {FT}
    J = 4 * ϕ0 * APAR / sqrt(1 + (4 * ϕ0 * APAR / max(Jmax, eps(FT)))^2)
    return J
end



"""
    inst_temp_scaling(
        T_canopy::FT,
        T_acclim::FT = T_canopy,
        To::FT,
        Ha::FT,
        Hd::FT,
        aS::FT,
        bS::FT,
        R::FT
    ) where {FT}

Given Vcmax or Jmax that have acclimated according to T_acclim, this function computes
the instantaneous temperature scaling factor f ∈ [0, ∞) for these maximum rates at the
instantaneous current temperature T_canopy. To is a reference temperature for the constants
and should be set to 298.15 K (25 °C). By default we assume that T_acclim = T_canopy.

The parameters (`Ha`, `Hd`, `aS`, `bS`) come from Kattge & Knorr (2007)
"""
function inst_temp_scaling(
    T_canopy::FT,
    T_acclim::FT,
    To::FT,
    Ha::FT,
    Hd::FT,
    aS::FT,
    bS::FT,
    R::FT,
) where {FT}
    T_acclim = T_acclim - FT(273.15)    # convert to C
    ΔS = aS - bS * T_acclim             # entropy term (J mol^-1 K^-1)

    # Arrhenius-type activation scaling factor
    f_act = arrhenius_function(T_canopy, To, R, Ha)

    # high temperature deactivation scaling factor
    num = 1 + exp((To * ΔS - Hd) / (R * To))
    den = 1 + exp((T_canopy * ΔS - Hd) / (R * T_canopy))
    f_deact = num / den

    return f_act * f_deact
end


"""
    inst_temp_scaling_rd(
        T_canopy::FT,
        To::FT,
        aRd::FT,
        bRd::FT
    ) where {FT}

Computes the instantaneous temperature scaling factor for dark respiration (Rd)
at canopy temperature `T_canopy` given reference temperature `To`, the first order
coefficient `aRd`, and the second order coefficient `bRd`.

Uses the log-quadratic functional form of Heskel et al. (2016)
https://www.pnas.org/doi/full/10.1073/pnas.1520282113
"""
function inst_temp_scaling_rd(T_canopy::FT, To::FT, aRd::FT, bRd::FT) where {FT}
    return exp(
        aRd * (T_canopy - To) +
        bRd * ((T_canopy - FT(273.15))^2 - (To - FT(273.15))^2),
    )
end

"""
    get_model_callbacks(component::AbstractCanopyComponent, canopy; kwargs...)

Creates the pmodel callback and returns it as a single element tuple of model callbacks;
we add the callback for the photosynthesis component,
and not for the conductance component(PModelConductance).

Note that the Δt passed here is an ITime because it is the Δt used in the simulation.
"""
function get_model_callbacks(
    component::PModel{FT},
    canopy;
    start_date,
    Δt,
) where {FT}
    pmodel_cb = make_PModel_callback(FT, start_date, float(Δt), canopy)
    return (pmodel_cb,)
end

"""
    compute_mj(
        is_c3::FT, T::FT, parameters) where {FT}

Computes the unitless factor `mj = (ci - Γstar)/(ci+2Γstar)` (for C3 plants)
and `mj = 1` for C4 plants, where the rubisco assimilation rate is Ac = Vcmax*mj.
"""
function compute_mj(is_c3::AbstractFloat, args...)
    return is_c3 > 0.5 ? c3_compute_mj(args...) : c4_compute_mj(args...)
end


function c3_compute_mj(Γstar::FT, ca_pp::FT, ci::FT, VPD::FT) where {FT}
    mj = (ci - Γstar) / (ci + 2 * Γstar) # eqn 11 in Stocker et al. (2020)
    return mj
end

function c4_compute_mj(::FT, ::FT, ::FT, ::FT) where {FT}
    return FT(1.0)
end

"""
    compute_mc(
        is_c3::FT, T::FT, parameters) where {FT}

Computes the unitless factor `mc = (ci - Γstar)/(ci+Kmm)` (for C3 plants)
and `mj = 1` for C4 plants, where the light assimilation rate is Aj = J/4 mj.
"""
function compute_mc(is_c3::AbstractFloat, args...)
    return is_c3 > 0.5 ? c3_compute_mc(args...) : c4_compute_mc(args...)
end

function c3_compute_mc(
    Γstar::FT,
    ca_pp::FT,
    ci::FT,
    VPD::FT,
    Kmm::FT,
) where {FT}
    mc = (ci - Γstar) / (ci + Kmm) # eqn 7 in Stocker et al. (2020)
    return mc
end

function c4_compute_mc(::FT, ::FT, ::FT, ::FT, ::FT) where {FT}
    return FT(1)
end
