export PModelParameters,
    PModelConstants,
    PModel,
    compute_full_pmodel_outputs,
    compute_optimal_capacities,
    compute_A0_and_χ

"""
    PModelParameters{FT<:AbstractFloat}

The required parameters for P-model (Stocker et al. 2020). Parameters are typically
tunable with considerable uncertainty.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PModelParameters{FT <: AbstractFloat}
    "Constant describing cost of maintaining electron transport (unitless)"
    cstar::FT
    "Ratio of unit costs of transpiration and carboxylation for C3 plants (unitless)"
    β_c3::FT
    "Ratio of unit costs of transpiration and carboxylation for C4 plants (unitless)"
    β_c4::FT
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
        once per day, α = 1 day/τ where τ is the acclimation timescale in days."""
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
    "Minimum ratio of √VPD/ξ"
    vpd_ratio_min::FT
    "Maximum ratio of Γ*/ca"
    Γ_ratio_max::FT
end

Base.eltype(::PModelParameters{FT}) where {FT} = FT
Base.eltype(::PModelConstants{FT}) where {FT} = FT

# make these custom structs broadcastable as tuples
Base.broadcastable(x::PModelParameters) = tuple(x)
Base.broadcastable(x::PModelConstants) = tuple(x)

"""
    PModelConstants(toml_dict::CP.ParamDict;
                    Kc25 = toml_dict["pmodel_Kc25"],
                    Ko25 = toml_dict["pmodel_Ko25"],
                    ΔHkc = toml_dict["pmodel_ΔHkc"],
                    ΔHko = toml_dict["pmodel_ΔHko"],
                    ΔHΓstar = toml_dict["pmodel_ΔHΓstar"],
                    Γstar25 = toml_dict["pmodel_Γstar25"],
                    Ha_Vcmax = toml_dict["pmodel_Ha_Vcmax"],
                    Hd_Vcmax = toml_dict["pmodel_Hd_Vcmax"],
                    aS_Vcmax = toml_dict["pmodel_aS_Vcmax"],
                    bS_Vcmax = toml_dict["pmodel_bS_Vcmax"],
                    Ha_Jmax = toml_dict["pmodel_Ha_Jmax"],
                    Hd_Jmax = toml_dict["pmodel_Hd_Jmax"],
                    aS_Jmax = toml_dict["pmodel_aS_Jmax"],
                    bS_Jmax = toml_dict["pmodel_bS_Jmax"],
                    oi = toml_dict["pmodel_oi"],
                    aRd = toml_dict["pmodel_aRd"],
                    bRd = toml_dict["pmodel_bRd"],
                    fC3 = toml_dict["pmodel_fC3"],
                    vpd_ratio_min = toml_dict["pmodel_vpd_ratio"],
                    Γ_ratio_max = toml_dict["pmodel_gamma_ratio"]
                )

Creates a `PModelConstants` object with default values for the P-model constants.
See Stocker et al. (2020) Table A2 and references within for more information.
"""
function PModelConstants(
    toml_dict::CP.ParamDict;
    Kc25 = toml_dict["pmodel_Kc25"],
    Ko25 = toml_dict["pmodel_Ko25"],
    ΔHkc = toml_dict["pmodel_ΔHkc"],
    ΔHko = toml_dict["pmodel_ΔHko"],
    ΔHΓstar = toml_dict["pmodel_ΔHΓstar"],
    Γstar25 = toml_dict["pmodel_Γstar25"],
    Ha_Vcmax = toml_dict["pmodel_Ha_Vcmax"],
    Hd_Vcmax = toml_dict["pmodel_Hd_Vcmax"],
    aS_Vcmax = toml_dict["pmodel_aS_Vcmax"],
    bS_Vcmax = toml_dict["pmodel_bS_Vcmax"],
    Ha_Jmax = toml_dict["pmodel_Ha_Jmax"],
    Hd_Jmax = toml_dict["pmodel_Hd_Jmax"],
    aS_Jmax = toml_dict["pmodel_aS_Jmax"],
    bS_Jmax = toml_dict["pmodel_bS_Jmax"],
    oi = toml_dict["pmodel_oi"],
    aRd = toml_dict["pmodel_aRd"],
    bRd = toml_dict["pmodel_bRd"],
    fC3 = toml_dict["pmodel_fC3"],
    vpd_ratio_min = toml_dict["pmodel_vpd_ratio"],
    Γ_ratio_max = toml_dict["pmodel_gamma_ratio"],
)
    # Note: physical constants are not exposed to the user
    FT = CP.float_type(toml_dict)
    return PModelConstants{FT}(
        toml_dict["universal_gas_constant"],
        Kc25,
        Ko25,
        toml_dict["kelvin_25C"],
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
        toml_dict["planck_constant"],
        toml_dict["light_speed"],
        toml_dict["avogadro_constant"],
        toml_dict["density_liquid_water"],
        vpd_ratio_min,
        Γ_ratio_max,
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
    T,
} <: AbstractPhotosynthesisModel{FT}
    "Required parameters for the P-model of Stocker et al. (2020)"
    parameters::OPFT
    "Constants for the P-model"
    constants::OPCT
    "Photosynthesis mechanism - 1 indicates C3, 0 indicates C4"
    fractional_c3::F
    "Time integrated prognostic vars"
    time_integrated_vars::T
end

Base.eltype(::PModel{FT, OPFT, OPCT, F, T}) where {FT, OPFT, OPCT, F, T} = FT

"""
    PModel{FT}(
        fractional_c3,
        parameters::PModelParameters{FT};
        constants::PModelConstants{FT} = PModelConstants(FT)
        binarize = false,
    )

Outer constructor for the PModel struct. This takes a PModelParameters struct which includes
parameters with considerable uncertainty. PModelConstants is constructed by default to the
default values, but if you know what you are doing, you can override with your own
constants.

If `binarize` is `true`, then `fractional_c3` is constrained to be either zero
or one.
"""
function PModel{FT}(
    fractional_c3,
    toml_dict::CP.ParamDict,
    parameters::PModelParameters{FT};
    constants::PModelConstants{FT} = PModelConstants(toml_dict),
    binarize = false,
) where {FT <: AbstractFloat}
    if binarize
        fractional_c3 = max.(min.(fractional_c3, FT(1)), FT(0))
        fractional_c3 = round.(fractional_c3)
    end
    F = typeof(fractional_c3)
    tiv = ClimaLand.time_integrated_variables(
        ClimaLand.TimeIntegratedVariable{FT}(;
            name = :acclimated,
            reduction = ClimaLand.RunningMean(
                IP.day(IP.InsolationParameters(FT)) / parameters.α,
            ),
            element_type = _pmodel_capacities_type(FT),
        ),
    )
    return PModel{FT, typeof(parameters), typeof(constants), F, typeof(tiv)}(
        parameters,
        constants,
        fractional_c3,
        tiv,
    )
end

"""
    ClimaLand.auxiliary_vars(model::PModel)
    ClimaLand.auxiliary_types(model::PModel)
    ClimaLand.auxiliary_domain_names(model::PModel)

Defines the auxiliary vars of the P-model:

- `instantaneous`: a NamedTuple with the canopy-level net photosynthesis (`An`),
    gross photosynthesis (`GPP`), dark respiration (`Rd`) and stomatal conductance
    to CO2 (`gs_co2`), computed each step from the acclimated capacities.
- `optimal`: a NamedTuple with keys `:ξ_c3`, `:ξ_c4`, `:Vcmax25_c3`,
    `:Vcmax25_c4`, `Jmax25_c3`, and `:Jmax25_c4`, holding the instantaneous optimal
    capacities — the target that the prognostic acclimated capacities
    `Y.canopy.photosynthesis.acclimated` relax toward (see `prognostic_vars`).
"""
# Element type of the optimal / acclimated capacity variables.
_pmodel_capacities_type(::Type{FT}) where {FT} = NamedTuple{
    (:ξ_c3, :ξ_c4, :Vcmax25_c3, :Vcmax25_c4, :Jmax25_c3, :Jmax25_c4),
    NTuple{6, FT},
}

ClimaLand.auxiliary_vars(model::PModel) = (:instantaneous, :optimal)
ClimaLand.auxiliary_types(model::PModel{FT}) where {FT} = (
    NamedTuple{(:Rd, :GPP, :An, :gs_co2), Tuple{FT, FT, FT, FT}},
    _pmodel_capacities_type(FT),
)
ClimaLand.auxiliary_domain_names(::PModel) = (:surface, :surface)

# The P-model's prognostic variable is the acclimated optimal capacities, a
# `RunningMean` time-integrated variable held in `Y` and advanced smoothly by the
# time-stepper, relaxing toward the instantaneous optimum `p.canopy.photosynthesis.optimal`
# weighted onto local solar noon (see `make_compute_exp_tendency`).

ClimaLand.prognostic_vars(m::PModel) =
    ClimaLand.time_integrated_prognostic_vars(m.time_integrated_vars)
ClimaLand.prognostic_types(m::PModel) =
    ClimaLand.time_integrated_prognostic_types(m.time_integrated_vars)
ClimaLand.prognostic_domain_names(m::PModel) =
    ClimaLand.time_integrated_prognostic_domain_names(m.time_integrated_vars)


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
        fractional_c3 = FT(1)
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
- "gs"            Stomatal conductance (mol m^-2 s^-1 Pa^-1)
- "vcmax25"       Vcmax normalized to 25°C via modified-Arrhenius type function (mol m^-2 s^-1)
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
    fractional_c3 = FT(1),
) where {FT}
    # Unpack parameters
    (; cstar, β_c3, β_c4) = parameters

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
        vpd_ratio_min,
        Γ_ratio_max,
    ) = constants

    # Convert ca from mol/mol to a partial pressure (Pa)
    ca_pp = ca * P_air

    # Compute intermediate values
    ϕ0_c3, ϕ0_c4 = intrinsic_quantum_yield(T_canopy, parameters)
    Γstar = co2_compensation_pmodel(T_canopy, To, P_air, R, ΔHΓstar, Γstar25)
    ηstar = compute_viscosity_ratio(T_canopy, To)
    Kmm = compute_Kmm(T_canopy, P_air, Kc25, Ko25, ΔHkc, ΔHko, To, R, oi)
    ξ_opt_c3 = sqrt(β_c3 * (Kmm + Γstar) / (Drel * ηstar))
    ξ_opt_c4 = sqrt(β_c4 * (Kmm + Γstar) / (Drel * ηstar))
    ci_c3 = intercellular_co2_pmodel(
        ξ_opt_c3,
        ca_pp,
        Γstar,
        VPD,
        vpd_ratio_min,
        Γ_ratio_max,
    )
    ci_c4 = intercellular_co2_pmodel(
        ξ_opt_c4,
        ca_pp,
        Γstar,
        VPD,
        vpd_ratio_min,
        Γ_ratio_max,
    )
    mj_c3 = c3_compute_mj(Γstar, ci_c3)
    mj_c4 = c4_compute_mj(Γstar, ci_c4) # This is misleading, because we approximate mj = 1
    mc_c3 = c3_compute_mc(Γstar, ci_c3, Kmm)
    mc_c4 = c4_compute_mc(Γstar, ci_c4, Kmm) # This is misleading, because we approximate mc = 1

    optimal_capacities = compute_optimal_capacities(
        parameters,
        constants,
        T_canopy,
        P_air,
        VPD,
        ca,
        βm,
        APAR,
    )
    (; ξ_c3, ξ_c4, Jmax25_c3, Jmax25_c4, Vcmax25_c3, Vcmax25_c4) =
        optimal_capacities

    blended_output = compute_blended_pmodel_photosynthesis(
        optimal_capacities,
        fractional_c3,
        P_air,
        VPD,
        ca,
        T_canopy,
        APAR,
        parameters,
        constants,
    )
    (; Rd, GPP, An, gs_co2) = blended_output

    return (;
        gpp = GPP*Mc,
        gammastar = Γstar,
        kmm = Kmm,
        ca = ca_pp,
        ns_star = ηstar,
        xi = blend(ξ_opt_c3, ξ_opt_c4, fractional_c3),
        mj = blend(mj_c3, mj_c4, fractional_c3),
        mc = blend(mc_c3, mc_c4, fractional_c3),
        ci = blend(ci_c3, ci_c4, fractional_c3),
        gs = gs_co2,
        vcmax25 = blend(Vcmax25_c3, Vcmax25_c4, fractional_c3),
        jmax25 = blend(Jmax25_c3, Jmax25_c4, fractional_c3),
        rd = Rd,
    )
end


"""
    compute_optimal_capacities(
        parameters::PModelParameters{FT},
        constants::PModelConstants{FT},
        T_canopy::FT,
        P_air::FT,
        VPD::FT,
        ca::FT,
        βm::FT,
        APAR_canopy_moles::FT,
    ) where {FT}

    compute_optimal_capacities(
        parameters::PModelParameters{FT},
        constants::PModelConstants{FT},
        thermo_params,
        T_canopy::FT,
        T_air::FT,
        P_air::FT,
        q_air::FT,
        ca::FT,
        βm::FT,
        APAR_canopy_moles::FT,
    ) where {FT}

Compute the instantaneous optimal photosynthetic capacities — the sensitivity of
stomatal conductance to dryness `ξ`, `Vcmax25`, and `Jmax25` (C3 and C4 variants) —
from the current environment, following Mengoli et al. (2022). These are the target
the acclimated capacities relax toward: the P-model `compute_exp_tendency!` advances
`Y.canopy.photosynthesis.acclimated` as a `RunningMean` of this instantaneous optimum,
applying the acclimation lag continuously through the time-stepper.

Args:
- `parameters`: PModelParameters object containing the model parameters.
- `constants`: PModelConstants object containing the model constants.
- `thermo_params`: Thermodynamic parameters, used to compute the VPD.
- `T_canopy`: Canopy temperature (K).
- `T_air`: Air temperature (K), used for the VPD.
- `P_air`: Ambient air pressure (Pa).
- `VPD`: Vapor pressure deficit (Pa).
- `q_air`: Specific humidity of the air (kg/kg).
- `ca`: Ambient CO2 concentration (mol/mol).
- `βm`: Soil moisture stress factor (unitless).
- `APAR_canopy_moles`: Absorbed photosynthetically active radiation (mol photons m^-2 s^-1).

Returns:
- NamedTuple with the instantaneous optimal `ξ`, `Vcmax25`, and `Jmax25` (C3/C4).

Reference:
Mengoli, G., Agustí-Panareda, A., Boussetta, S., Harrison, S. P., Trotta, C., & Prentice, I. C. (2022).
Ecosystem photosynthesis in land-surface models: A first-principles approach incorporating acclimation.
Journal of Advances in Modeling Earth Systems, 14, e2021MS002767. https://doi.org/10.1029/2021MS002767
"""
function compute_optimal_capacities(
    parameters::PModelParameters{FT},
    constants::PModelConstants{FT},
    thermo_params,
    T_canopy::FT,
    T_air::FT,
    P_air::FT,
    q_air::FT,
    ca::FT,
    βm::FT,
    APAR_canopy_moles::FT,
) where {FT}
    VPD = Thermodynamics.vapor_pressure_deficit(
        thermo_params,
        T_air,
        P_air,
        q_air,
    )
    return compute_optimal_capacities(
        parameters,
        constants,
        T_canopy,
        P_air,
        VPD,
        ca,
        βm,
        APAR_canopy_moles,
    )
end

function compute_optimal_capacities(
    parameters::PModelParameters{FT},
    constants::PModelConstants{FT},
    T_canopy::FT,
    P_air::FT,
    VPD::FT,
    ca::FT,
    βm::FT,
    APAR_canopy_moles::FT,
) where {FT}
    # Unpack parameters
    (; cstar, β_c3, β_c4) = parameters

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
        oi,
        ρ_water,
        vpd_ratio_min,
        Γ_ratio_max,
    ) = constants
    # Compute intermediate values
    Γstar = co2_compensation_pmodel(T_canopy, To, P_air, R, ΔHΓstar, Γstar25)
    ηstar = compute_viscosity_ratio(T_canopy, To)
    Kmm = compute_Kmm(T_canopy, P_air, Kc25, Ko25, ΔHkc, ΔHko, To, R, oi)

    # convert ca from mol/mol to Pa
    ca_pp = ca * P_air

    ξ_opt_c3 = sqrt(β_c3 * (Kmm + Γstar) / (Drel * ηstar))
    ξ_opt_c4 = sqrt(β_c4 * (Kmm + Γstar) / (Drel * ηstar))
    ci_c3 = intercellular_co2_pmodel(
        ξ_opt_c3,
        ca_pp,
        Γstar,
        VPD,
        vpd_ratio_min,
        Γ_ratio_max,
    )
    ci_c4 = intercellular_co2_pmodel(
        ξ_opt_c4,
        ca_pp,
        Γstar,
        VPD,
        vpd_ratio_min,
        Γ_ratio_max,
    )

    ϕ0_c3, ϕ0_c4 = intrinsic_quantum_yield(T_canopy, parameters)

    mj_c3 = c3_compute_mj(Γstar, ci_c3)
    mj_c4 = c4_compute_mj(Γstar, ci_c4) # This is misleading, because we approximate mj = 1
    L_c3 = sqrt(1-min(cstar/max(mj_c3, eps(FT)), 1)^FT(2/3)) # If mj < cstar, L = 0
    L_c4 = sqrt(1-min(cstar/max(mj_c4, eps(FT)), 1)^FT(2/3)) # If mj < cstar, L = 0
    mj_over_mc_c3 = (ci_c3 + Kmm)/(ci_c3+2*Γstar)
    mj_over_mc_c4 = 1
    # Vcmax opt is from Equation 10a of Ren et al.
    Vcmax_opt_c3 = βm * ϕ0_c3 * APAR_canopy_moles * mj_over_mc_c3 * L_c3
    Vcmax_opt_c4 = βm * ϕ0_c4 * APAR_canopy_moles * mj_over_mc_c4 * L_c4

    inst_temp_scaling_factor_Vcmax = inst_temp_scaling(
        T_canopy,
        T_canopy,
        To,
        Ha_Vcmax,
        Hd_Vcmax,
        aS_Vcmax,
        bS_Vcmax,
        R,
    )
    Vcmax25_opt_c3 = Vcmax_opt_c3 / inst_temp_scaling_factor_Vcmax
    Vcmax25_opt_c4 = Vcmax_opt_c4 / inst_temp_scaling_factor_Vcmax

    Jmax_opt_c3 =
        βm * 4 * ϕ0_c3 * APAR_canopy_moles * L_c3/sqrt(1 - (βm*L_c3)^FT(2))
    Jmax_opt_c4 =
        βm * 4 * ϕ0_c4 * APAR_canopy_moles * L_c4/sqrt(1 - (βm*L_c4)^FT(2))
    inst_temp_scaling_factor_Jmax = inst_temp_scaling(
        T_canopy,
        T_canopy,
        To,
        Ha_Jmax,
        Hd_Jmax,
        aS_Jmax,
        bS_Jmax,
        R,
    )
    Jmax25_opt_c3 = Jmax_opt_c3 / inst_temp_scaling_factor_Jmax
    Jmax25_opt_c4 = Jmax_opt_c4 / inst_temp_scaling_factor_Jmax
    return (;
        ξ_c3 = ξ_opt_c3,
        ξ_c4 = ξ_opt_c4,
        Vcmax25_c3 = Vcmax25_opt_c3,
        Vcmax25_c4 = Vcmax25_opt_c4,
        Jmax25_c3 = Jmax25_opt_c3,
        Jmax25_c4 = Jmax25_opt_c4,
    )
end

# Normalization for the solar-noon acclimation window `exp(κ (cosθ - 1))`: its daily
# mean, ⟨exp(κ (cosθ - 1))⟩_θ = exp(-κ) I₀(κ). Dividing the window by this makes its
# daily mean exactly 1, so the noon weighting reshapes the sub-daily forcing without
# changing the acclimation timescale τ. Evaluated once at setup by quadrature over the
# day, which avoids a SpecialFunctions dependency for I₀.
function _noon_window_norm(κ::FT) where {FT}
    n = 2048
    acc = zero(FT)
    for i in 0:(n - 1)
        acc += exp(κ * (cos(FT(2π) * i / n) - 1))
    end
    return acc / n
end

# Current time as seconds since UTC midnight, used to place the solar-noon window.
# Mirrors the ITime/epoch handling in `drivers.jl`: `date(t)` is used when `t` carries
# an epoch (the usual case), otherwise `start_date` provides the calendar reference.
function _seconds_of_day_utc(t, start_date, ::Type{FT}) where {FT}
    current = if t isa ITime
        isnothing(t.epoch) ? start_date + t.counter * t.period : date(t)
    else
        start_date + Second(round(Int, float(t)))
    end
    return FT(
        Hour(current).value * 3600 +
        Minute(current).value * 60 +
        Second(current).value,
    )
end

"""
    ClimaLand.make_compute_exp_tendency(component::PModel, canopy)

Advances the acclimated optimal capacities `Y.canopy.photosynthesis.acclimated` as a
`RunningMean` time-integrated variable, relaxing toward the instantaneous optimal
capacities `p.canopy.photosynthesis.optimal` (recomputed every step from the current
environment) but with the relaxation rate weighted by a smooth window centered on
local solar noon:

    dacclimated/dt = w(t) (optimal - acclimated) / τ,    τ = 1 day / α,

    w(t) = exp(κ (cos θ - 1)) / ⟨exp(κ (cos θ - 1))⟩,    θ = 2π (t_UTC - noon) / day.

`w` is the von Mises "midday window": it peaks at solar noon (computed per column
from longitude, neglecting the equation of time) and is ≈ 0 at night, so the
acclimation samples midday conditions as the P-model of Mengoli et al. (2022) and
`main` intended. Since `optimal ∝ APAR` vanishes at night, an unweighted
relaxation would instead average the optimum over the whole diurnal cycle and bias
the acclimated capacities low in a daylength-dependent way. `w` is normalized to unit
daily mean, so weighting reshapes the sub-daily forcing without changing the
acclimation timescale `τ`: `α = 1` (τ = 1 day) gives fast acclimation and `α → 0`
(τ → ∞) freezes the acclimated capacities. Because they live in `Y`, the acclimation
is advanced smoothly by the time-stepper and is checkpoint/restart-safe.
"""
function ClimaLand.make_compute_exp_tendency(
    component::PModel{FT},
    canopy,
) where {FT}
    reduction = component.time_integrated_vars.acclimated.reduction
    seconds_in_a_day = IP.day(IP.InsolationParameters(FT))
    # Per-column solar-noon time (seconds UTC) from longitude; neglects the equation of
    # time (≤ ~20 min error), matching the removed local-noon sampling of `main`.
    longitude = get_long(canopy.domain.space.surface)
    local_noon = @. seconds_in_a_day * (FT(1 / 2) - longitude / 360)
    # Gaussian-equivalent half-width of the midday window (s). At 1 h the acclimation
    # target tracks the exact solar-noon optimum to ~0.95 (vs the daylength-diluted
    # ~0.2-0.4 of an unweighted diurnal mean), while the ~16 min error from neglecting
    # the equation of time shifts it by <0.3%. κ matches exp(-κ θ²/2) near noon to a
    # Gaussian of this width.
    noon_window_seconds = FT(3600)
    κ = (seconds_in_a_day / FT(2π))^2 / noon_window_seconds^2
    inv_norm = 1 / _noon_window_norm(κ)
    start_date = try
        canopy.boundary_conditions.atmos.start_date
    catch
        nothing
    end
    ω = FT(2π) / seconds_in_a_day
    function compute_exp_tendency!(dY, Y, p, t)
        tod = _seconds_of_day_utc(t, start_date, FT)
        @. dY.canopy.photosynthesis.acclimated =
            apply_time_reduction(
                p.canopy.photosynthesis.optimal,
                Y.canopy.photosynthesis.acclimated,
                reduction,
            ) ⊠ (exp(κ * (cos(ω * (tod - local_noon)) - 1)) * inv_norm)
    end
    return compute_exp_tendency!
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

    # drivers
    P_air = p.drivers.P
    T_air = p.drivers.T
    q_air = p.drivers.q
    c_co2_air = p.drivers.c_co2
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    thermo_params = LP.thermodynamic_parameters(canopy.earth_param_set)
    λ_γ_PAR = canopy.radiative_transfer.parameters.λ_γ_PAR
    fAPAR = p.canopy.radiative_transfer.par.abs
    PAR = p.canopy.radiative_transfer.par_d
    βm = p.canopy.soil_moisture_stress.βm
    APAR_canopy_moles = @. lazy(
        compute_APAR_canopy_moles(
            fAPAR,
            PAR,
            λ_γ_PAR,
            constants.lightspeed,
            constants.planck_h,
            constants.N_a,
        ),
    )
    @. p.canopy.photosynthesis.optimal = compute_optimal_capacities(
        parameters,
        constants,
        thermo_params,
        T_canopy,
        T_air,
        P_air,
        q_air,
        c_co2_air,
        βm,
        APAR_canopy_moles,
    )

    @. p.canopy.photosynthesis.instantaneous =
        compute_blended_pmodel_photosynthesis(
            Y.canopy.photosynthesis.acclimated,
            model.fractional_c3,
            P_air,
            T_air,
            q_air,
            c_co2_air,
            T_canopy,
            APAR_canopy_moles,
            parameters,
            constants,
            thermo_params,
        )
end

function compute_blended_pmodel_photosynthesis(
    acclimated,
    fractional_c3::FT,
    P_air::FT,
    T_air::FT,
    q_air::FT,
    c_co2_air::FT,
    T_canopy::FT,
    APAR_canopy_moles::FT,
    parameters,
    constants,
    thermo_params,
) where {FT}

    VPD = Thermodynamics.vapor_pressure_deficit(
        thermo_params,
        T_air,
        P_air,
        q_air,
    )
    return compute_blended_pmodel_photosynthesis(
        acclimated,
        fractional_c3,
        P_air,
        VPD,
        c_co2_air,
        T_canopy,
        APAR_canopy_moles,
        parameters,
        constants,
    )
end

function compute_blended_pmodel_photosynthesis(
    acclimated,
    fractional_c3::FT,
    P_air::FT,
    VPD::FT,
    c_co2_air::FT,
    T_canopy::FT,
    APAR_canopy_moles::FT,
    parameters,
    constants,
) where {FT}
    (; Vcmax25_c3, Vcmax25_c4, Jmax25_c3, Jmax25_c4, ξ_c3, ξ_c4) = acclimated
    ca_pp = c_co2_air * P_air # partial pressure of co2
    Γstar = co2_compensation_pmodel(
        T_canopy,
        constants.To,
        P_air,
        constants.R,
        constants.ΔHΓstar,
        constants.Γstar25,
    )
    ci_c3 = intercellular_co2_pmodel(
        ξ_c3,
        ca_pp,
        Γstar,
        VPD,
        constants.vpd_ratio_min,
        constants.Γ_ratio_max,
    )
    ci_c4 = intercellular_co2_pmodel(
        ξ_c4,
        ca_pp,
        Γstar,
        VPD,
        constants.vpd_ratio_min,
        constants.Γ_ratio_max,
    )
    Kmm = compute_Kmm(
        T_canopy,
        P_air,
        constants.Kc25,
        constants.Ko25,
        constants.ΔHkc,
        constants.ΔHko,
        constants.To,
        constants.R,
        constants.oi,
    )
    # Fast responses: scale the acclimated capacities to the current canopy
    # temperature to get the instantaneous Vcmax and Jmax, and the assimilation rates
    inst_temp_scaling_Jmax_factor = inst_temp_scaling(
        T_canopy,
        T_canopy,
        constants.To,
        constants.Ha_Jmax,
        constants.Hd_Jmax,
        constants.aS_Jmax,
        constants.bS_Jmax,
        constants.R,
    )
    Jmax_c3 = Jmax25_c3 * inst_temp_scaling_Jmax_factor
    Jmax_c4 = Jmax25_c4 * inst_temp_scaling_Jmax_factor
    J_c3 = electron_transport_pmodel(
        c3_intrinsic_quantum_yield(T_canopy, parameters),
        APAR_canopy_moles,
        Jmax_c3,
    )
    J_c4 = electron_transport_pmodel(
        c4_intrinsic_quantum_yield(T_canopy, parameters),
        APAR_canopy_moles,
        Jmax_c4,
    )

    inst_temp_scaling_Vcmax_factor = inst_temp_scaling(
        T_canopy,
        T_canopy,
        constants.To,
        constants.Ha_Vcmax,
        constants.Hd_Vcmax,
        constants.aS_Vcmax,
        constants.bS_Vcmax,
        constants.R,
    )
    Vcmax_c3 = Vcmax25_c3 * inst_temp_scaling_Vcmax_factor
    Vcmax_c4 = Vcmax25_c4 * inst_temp_scaling_Vcmax_factor
    Ac_c3 = Vcmax_c3 * c3_compute_mc(Γstar, ci_c3, Kmm)
    Ac_c4 = Vcmax_c4 * c4_compute_mc(Γstar, ci_c4, Kmm)
    # light limited assimilation rate
    # c3 or c4 is reflected in the value of mj and J
    Aj_c3 = J_c3 * FT(0.25) * c3_compute_mj(Γstar, ci_c3)
    Aj_c4 = J_c4 * FT(0.25) * c4_compute_mj(Γstar, ci_c4)
    # dark respiration
    # Here we make an assumption about how to relate Rd25 to the acclimated Vcmax25
    # To extend to C4, defined `compute_dark_respiration_pmodel() which dispatches off of the is_c3 field
    # This function below would become c3_dark_respiration_pmodel
    Rd = blend(
        constants.fC3 *
        Vcmax25_c3 *
        inst_temp_scaling_rd(
            T_canopy,
            constants.To,
            constants.aRd,
            constants.bRd,
        ),
        constants.fC3 *
        Vcmax25_c4 *
        inst_temp_scaling_rd(
            T_canopy,
            constants.To,
            constants.aRd,
            constants.bRd,
        ),
        fractional_c3,
    )

    # Note: net_photosynthesis applies the moisture stress to GPP, but since the P-model already applies
    # this factor to Vcmax and Jmax, we do not apply it again here
    GPP = blend(
        gross_photosynthesis(Ac_c3, Aj_c3),
        gross_photosynthesis(Ac_c4, Aj_c4),
        fractional_c3,
    )
    An = net_photosynthesis(GPP, Rd)
    # gs blend
    χ_c3 = clamp(ci_c3 / ca_pp, FT(0), FT(1))
    χ_c4 = clamp(ci_c4 / ca_pp, FT(0), FT(1))
    gs_co2_c3 =
        gs_co2_pmodel(χ_c3, c_co2_air, gross_photosynthesis(Ac_c3, Aj_c3))
    gs_co2_c4 =
        gs_co2_pmodel(χ_c4, c_co2_air, gross_photosynthesis(Ac_c4, Aj_c4))
    gs_co2 = blend(gs_co2_c3, gs_co2_c4, fractional_c3) # Assumes C3 and C4 plants act in parallel
    return (; Rd, GPP, An, gs_co2)
end

get_Vcmax25_canopy(Y, p, m::PModel) = @. lazy(
    blend(
        Y.canopy.photosynthesis.acclimated.Vcmax25_c3,
        Y.canopy.photosynthesis.acclimated.Vcmax25_c4,
        m.fractional_c3,
    ),
)

get_Vcmax25_leaf(Y, p, m::PModel) = @. lazy(
    blend(
        Y.canopy.photosynthesis.acclimated.Vcmax25_c3,
        Y.canopy.photosynthesis.acclimated.Vcmax25_c4,
        m.fractional_c3,
    ) /
    max(p.canopy.biomass.area_index.leaf, sqrt(eps(eltype(m.constants)))),
)
get_Rd_canopy(p, m::PModel) = p.canopy.photosynthesis.instantaneous.Rd
get_Rd_leaf(p, m::PModel) = @. lazy(
    p.canopy.photosynthesis.instantaneous.Rd /
    max(p.canopy.biomass.area_index.leaf, sqrt(eps(eltype(m.constants)))),
)
get_An_canopy(p, m::PModel) = p.canopy.photosynthesis.instantaneous.An
get_An_leaf(p, m::PModel) = @. lazy(
    p.canopy.photosynthesis.instantaneous.An /
    max(p.canopy.biomass.area_index.leaf, sqrt(eps(eltype(m.constants)))),
)

get_GPP(p, m::PModel) = p.canopy.photosynthesis.instantaneous.GPP

function get_J_over_Jmax(Y, p, canopy, m::PModel)
    Jmax_c3, Jmax_c4 = compute_Jmax_canopy(Y, p, canopy, m) # lazy
    J_c3, J_c4 = compute_J_canopy(Y, p, canopy, m) # lazy
    FT = eltype(m.constants)
    return @. lazy(
        blend(
            J_c3 / max(Jmax_c3, sqrt(eps(FT))),
            J_c4 / max(Jmax_c4, sqrt(eps(FT))),
            m.fractional_c3,
        ),
    )
end

function compute_Jmax_canopy(Y, p, canopy, m::PModel) # used internally to pmodel photosynthesis as a helper function
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    constants = m.constants
    inst_temp_scaling_factor = @. lazy(
        inst_temp_scaling(
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
    return @. (
        lazy(
            Y.canopy.photosynthesis.acclimated.Jmax25_c3 *
            inst_temp_scaling_factor,
        ),
        lazy(
            Y.canopy.photosynthesis.acclimated.Jmax25_c4 *
            inst_temp_scaling_factor,
        ),
    )
end

function compute_J_canopy(Y, p, canopy, m::PModel) # used internally to pmodel photosynthesis as a helper function
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    earth_param_set = canopy.earth_param_set
    f_abs_par = p.canopy.radiative_transfer.par.abs
    par_d = p.canopy.radiative_transfer.par_d
    (; λ_γ_PAR,) = canopy.radiative_transfer.parameters
    c = LP.light_speed(earth_param_set)
    planck_h = LP.planck_constant(earth_param_set)
    N_a = LP.avogadro_constant(earth_param_set)
    APAR_canopy_moles = @. lazy(
        compute_APAR_canopy_moles(f_abs_par, par_d, λ_γ_PAR, c, planck_h, N_a),
    )

    Jmax_canopy_c3, Jmax_canopy_c4 = compute_Jmax_canopy(Y, p, canopy, m)
    parameters = m.parameters
    constants = m.constants
    return @. (
        lazy(
            electron_transport_pmodel(
                c3_intrinsic_quantum_yield(T_canopy, parameters),
                APAR_canopy_moles,
                Jmax_canopy_c3,
            ),
        ),
        lazy(
            electron_transport_pmodel(
                c4_intrinsic_quantum_yield(T_canopy, parameters),
                APAR_canopy_moles,
                Jmax_canopy_c4,
            ),
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
    intrinsic_quantum_yield(T::FT, parameters) where {FT}

Computes the intrinsic quantum yield of photosystem II for both c3 and c4.
"""
function intrinsic_quantum_yield(T::FT, parameters) where {FT}
    return (
        c3_intrinsic_quantum_yield(T, parameters),
        c4_intrinsic_quantum_yield(T, parameters),
    )
end

function c3_intrinsic_quantum_yield(T::FT, parameters) where {FT}
    return parameters.temperature_dep_yield ?
           quadratic_intrinsic_quantum_yield(
        T,
        parameters.ϕa0_c3,
        parameters.ϕa1_c3,
        parameters.ϕa2_c3,
    ) : parameters.ϕ0_c3
end

function c4_intrinsic_quantum_yield(T::FT, parameters) where {FT}
    return parameters.temperature_dep_yield ?
           quadratic_intrinsic_quantum_yield(
        T,
        parameters.ϕa0_c4,
        parameters.ϕa1_c4,
        parameters.ϕa2_c4,
    ) : parameters.ϕ0_c4
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
    if T < FT(0)
        return FT(0)
    else
        ϕ = ϕa0 + ϕa1 * T + ϕa2 * T^FT(2)
        return min(max(ϕ, FT(0)), FT(1)) # Clip to [0,1]
    end
end


"""
    patek_dimensionless_viscosity(
        T::FT,
    ) where {FT}

Computes the viscosity of water in Pa s given temperature T (K),
as found by J. Pátek, J. Hrubý, J. Klomfar, M. Součková, and A. H. Harvey,
J. Phys. Chem. Ref. Data https://doi.org/10.1063/1.304357538, 21 (2009); as presented
in Equation 37 of Huber et al. (2009) [https://doi.org/10.1063/1.3088050].

This formula should not be used outside the range of 253.15 K ≤ T ≤ 383.15 K
"""
function patek_dimensionless_viscosity(T::FT) where {FT}
    a1 = FT(280.68)
    b1 = FT(-1.9)
    a2 = FT(511.45)
    b2 = FT(-7.7)
    a3 = FT(61.131)
    b3 = FT(-19.6)
    a4 = FT(0.45903)
    b4 = FT(-40)
    x = T/FT(300)
    μ = a1*x^b1 + a2*x^b2 + a3*x^b3 + a4*x^b4
    return μ
end

"""
    compute_viscosity_ratio(
        T::FT,
        To::FT,
    ) where {FT}

Computes η*, the ratio of the viscosity of water at temperature T to that at To = 25˚C.
"""
function compute_viscosity_ratio(T::FT, To::FT) where {FT}
    μ = patek_dimensionless_viscosity(max(min(T, FT(383.15)), FT(253.15)))
    μ25 = patek_dimensionless_viscosity(To)
    return μ/μ25
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
    intercellular_co2_pmodel(ξ::FT, ca_pp::FT, Γstar::FT, VPD::FT, vpd_ratio_min::FT,
    Γ_ratio_max::FT) where {FT}

Computes the intercellular co2 concentration (`ci`) as a function of the
optimal `ξ` (sensitivity to dryness), `ca_pp` (ambient CO2 partial pressure),
`Γstar` (CO2 compensation point), and `VPD` (vapor pressure deficit).

To prevent ci = ca, we limit both √VPD/ξ and Γ⋆/ca_pp. 
"""
function intercellular_co2_pmodel(
    ξ::FT,
    ca_pp::FT,
    Γstar::FT,
    VPD::FT,
    vpd_ratio_min::FT,
    Γ_ratio_max::FT,
) where {FT}
    # Guard ξ = 0 (e.g. before the optimal ξ is initialized); the ξ→0 limit is ci = Γstar
    x = max(sqrt(VPD) / max(ξ, eps(FT)), vpd_ratio_min)
    y = min(Γstar / ca_pp, Γ_ratio_max)
    return ca_pp * (1 + y * x) / (1 + x)
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
    compute_LUE(
        ϕ0::FT,
        βm::FT,
        mprime::FT,
        Mc::FT
    ) where {FT}

Computes light use efficiency (LUE) in kg C/mol from intrinsic quantum yield (`ϕ0`),
moisture stress factor (`β`), and a Jmax modified capacity (`mprime`); see Eqn 17 and 19
in Stocker et al. (2020). Mc is the molar mass of carbon (kg/mol) = 0.0120107 kg/mol.
"""
function compute_LUE(ϕ0::FT, βm::FT, mprime::FT, Mc::FT) where {FT}
    return ϕ0 * βm * mprime * Mc
end

"""
    compute_A0_and_χ(
        fractional_c3::FT,
        parameters::PModelParameters{FT},
        constants::PModelConstants{FT},
        T_air::FT,
        P_air::FT,
        VPD::FT,
        ca::FT,
        PPFD::FT,
        βm::FT,
        vpd_gs::FT,
    ) where {FT}

Compute potential GPP (A0) and ci/ca ratio used in the optimal LAI model (Zhou et al. 2025).

This function computes the potential GPP assuming fAPAR = 1 (full light absorption),
which represents the maximum carbon assimilation possible under given environmental
conditions.

`A0` responds to the instantaneous vapor pressure deficit, while `χ` is the
growing-season value the water-limitation term of `LAI_max` calls for, and so is
evaluated at the growing-season mean VPD `vpd_gs`. Both share the same optimal `ξ`,
which depends on temperature and pressure only.

# Arguments
- `fractional_c3::FT`: Photosynthesis mechanism (1 for C3, 0 for C4)
- `parameters`: PModelParameters containing cstar, beta, etc.
- `constants`: PModelConstants containing physical constants
- `earth_param_set`: Additional physical constants
- `T_air::FT`: Air temperature (K)
- `P_air::FT`: Atmospheric pressure (Pa)
- `q_air::FT`: Atmospheric specific humidity (kg/kg)
- `ca::FT`: Ambient CO2 concentration (mol/mol)
- `PPFD::FT`: Downwelling radiative flux in PAR band at the surface (moles photons/m^2/s)
- `βm::FT`: Soil moisture stress factor (dimensionless, 0-1)
- `vpd_gs::FT`: Growing-season mean vapor pressure deficit (Pa)

# Returns
- NamedTuple with `A0`, the potential GPP with fAPAR=1 (mol C m^-2 s^-1), and
  `χ = ci/ca` at the growing-season VPD (unitless)
"""
function compute_A0_and_χ(
    fractional_c3::FT,
    parameters::PModelParameters{FT},
    constants::PModelConstants{FT},
    earth_param_set,
    T_air::FT,
    P_air::FT,
    q_air::FT,
    ca::FT,
    PPFD::FT,
    βm::FT,
    vpd_gs::FT,
) where {FT}
    (; cstar, β_c3, β_c4) = parameters
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
        Mc,
        oi,
        ρ_water,
        vpd_ratio_min,
        Γ_ratio_max,
    ) = constants

    # Convert ca from mol/mol to partial pressure (Pa)
    ca_pp = ca * P_air
    VPD = Thermodynamics.vapor_pressure_deficit(
        LP.thermodynamic_parameters(earth_param_set),
        T_air,
        P_air,
        q_air,
    )
    # Compute P-model intermediate values
    ϕ0_c3, ϕ0_c4 = intrinsic_quantum_yield(T_air, parameters)
    Γstar = co2_compensation_pmodel(T_air, To, P_air, R, ΔHΓstar, Γstar25)
    ηstar = compute_viscosity_ratio(T_air, To)
    Kmm = compute_Kmm(T_air, P_air, Kc25, Ko25, ΔHkc, ΔHko, To, R, oi)

    # Compute optimal ξ and intercellular CO2
    ξ_opt_c3 = sqrt(β_c3 * (Kmm + Γstar) / (Drel * ηstar))
    ξ_opt_c4 = sqrt(β_c4 * (Kmm + Γstar) / (Drel * ηstar))
    ci_c3 = intercellular_co2_pmodel(
        ξ_opt_c3,
        ca_pp,
        Γstar,
        VPD,
        vpd_ratio_min,
        Γ_ratio_max,
    )
    ci_c4 = intercellular_co2_pmodel(
        ξ_opt_c4,
        ca_pp,
        Γstar,
        VPD,
        vpd_ratio_min,
        Γ_ratio_max,
    )

    mj_c3 = c3_compute_mj(Γstar, ci_c3)
    mj_c4 = c4_compute_mj(Γstar, ci_c4) # This is misleading, because we approximate mj = 1
    L_c3 = sqrt(1-min(cstar/max(mj_c3, eps(FT)), 1)^FT(2/3)) # If mj < cstar, L = 0
    L_c4 = sqrt(1-min(cstar/max(mj_c4, eps(FT)), 1)^FT(2/3)) # If mj < cstar, L = 0
    # Compute LUE with actual βm (soil moisture stress)
    LUE_daily_c3 = compute_LUE(ϕ0_c3, βm, mj_c3*L_c3, Mc)
    LUE_daily_c4 = compute_LUE(ϕ0_c4, βm, mj_c4*L_c4, Mc)

    # χ = ci/ca at the growing-season VPD, as the LAI_max water limitation requires
    ci_gs_c3 = intercellular_co2_pmodel(
        ξ_opt_c3,
        ca_pp,
        Γstar,
        vpd_gs,
        vpd_ratio_min,
        Γ_ratio_max,
    )
    ci_gs_c4 = intercellular_co2_pmodel(
        ξ_opt_c4,
        ca_pp,
        Γstar,
        vpd_gs,
        vpd_ratio_min,
        Γ_ratio_max,
    )

    # Potential GPP = PPFD * LUE (fAPAR = 1 is implicit in using full PPFD); LUE is
    # in kg C per mol photon, so dividing by Mc returns it in mol C m^-2 s^-1.
    return (;
        A0 = PPFD * blend(LUE_daily_c3, LUE_daily_c4, fractional_c3) / Mc,
        χ = blend(ci_gs_c3, ci_gs_c4, fractional_c3) / ca_pp,
    )
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
which expects T to be in Celsius:
lnR = constant + aT + bT^2, or
ln(R/R25) = a(T-To) + b(T^2-To^2)
"""
function inst_temp_scaling_rd(T_canopy::FT, To::FT, aRd::FT, bRd::FT) where {FT}
    return exp(
        aRd * (T_canopy - To) +
        bRd * ((T_canopy - FT(273.15))^2 - (To - FT(273.15))^2),
    )
end


"""
    c3_compute_mj(args...)

Computes the unitless factor `mj = (ci - Γstar)/(ci+2Γstar)` (for C3 plants);
"""
function c3_compute_mj(Γstar::FT, ci::FT) where {FT}
    mj = (ci - Γstar) / (ci + 2 * Γstar) # eqn 11 in Stocker et al. (2020)
    return mj
end
"""
    c4_compute_mj(args...)

Computes the unitless factor `mj = 1` for C4 plants
for both c3 and c4.
"""
function c4_compute_mj(::FT, ::FT) where {FT}
    return FT(1.0)
end

"""
    c3_compute_mc(args...)

Computes the unitless factor `mc = (ci - Γstar)/(ci+Kmm)` (for C3 plants).
"""
function c3_compute_mc(Γstar::FT, ci::FT, Kmm::FT) where {FT}
    mc = (ci - Γstar) / (ci + Kmm) # eqn 7 in Stocker et al. (2020)
    return mc
end

"""
    c4_compute_mc(args...)

Computes the unitless factor `mc = 1` (for C4 plants).
"""
function c4_compute_mc(::FT, ::FT, ::FT) where {FT}
    return FT(1)
end

"""
    blend(x1, x2, w)

Compute a weighted average of `x1` and `x2` weighted by
0 <= `w` <= 1.
"""
function blend(x1, x2, w)
    return w * x1 + (1 - w) * x2
end
