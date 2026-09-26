module Parameters
import Thermodynamics.Parameters.ThermodynamicsParameters
import Insolation.Parameters.InsolationParameters
import SurfaceFluxes.Parameters.SurfaceFluxesParameters
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP

abstract type AbstractLandParameters end
const ALP = AbstractLandParameters

Base.@kwdef struct LandParameters{FT, TP, SFP, IP} <: ALP
    K_therm::FT
    ρ_cloud_liq::FT
    ρ_cloud_ice::FT
    cp_l::FT
    cp_i::FT
    T_0::FT
    LH_v0::FT
    LH_s0::FT
    Stefan::FT
    T_freeze::FT
    grav::FT
    MSLP::FT
    D_vapor::FT
    gas_constant::FT
    molmass_water::FT
    h_Planck::FT
    light_speed::FT
    avogad::FT
    thermo_params::TP
    surf_flux_params::SFP
    insol_params::IP
end

Base.eltype(::LandParameters{FT}) where {FT} = FT
Base.broadcastable(ps::LandParameters) = tuple(ps)

# wrapper methods:
P_ref(ps::ALP) = ps.MSLP
K_therm(ps::ALP) = ps.K_therm
ρ_cloud_liq(ps::ALP) = ps.ρ_cloud_liq
ρ_cloud_ice(ps::ALP) = ps.ρ_cloud_ice
cp_l(ps::ALP) = ps.cp_l
cp_i(ps::ALP) = ps.cp_i
T_0(ps::ALP) = ps.T_0
LH_v0(ps::ALP) = ps.LH_v0
LH_s0(ps::ALP) = ps.LH_s0
Stefan(ps::ALP) = ps.Stefan
T_freeze(ps::ALP) = ps.T_freeze
grav(ps::ALP) = ps.grav
D_vapor(ps::ALP) = ps.D_vapor
gas_constant(ps::ALP) = ps.gas_constant
molar_mass_water(ps::ALP) = ps.molmass_water
planck_constant(ps::ALP) = ps.h_Planck
avogadro_constant(ps::ALP) = ps.avogad
light_speed(ps::ALP) = ps.light_speed
# Derived parameters
LH_f0(ps::ALP) = LH_s0(ps) - LH_v0(ps)
ρ_m_liq(ps::ALP) = ρ_cloud_liq(ps) / molar_mass_water(ps)
# Dependency parameter wrappers
thermodynamic_parameters(ps::ALP) = ps.thermo_params
surface_fluxes_parameters(ps::ALP) = ps.surf_flux_params
insolation_parameters(ps::ALP) = ps.insol_params

"""
    StableLimitedUniversalFunctionParams{FT, P} <: UF.AbstractUniversalFunctionParameters{FT}

Monin-Obukhov universal functions `base` (e.g. Gryanik et al. 2020) with
the stability parameter entering the flux-profile relations limited to
`ζ ≤ ζ_max_stable` under stable stratification, as in CLM (`zetamaxstable`;
2 in CLM4.5 and in CLM5/6 configurations with biomass heat storage or the
Meier et al. (2022) roughness, 0.5 in the CLM5.0 default). The limit
represents turbulent mixing that persists in very stable conditions
(intermittent turbulence, submeso motions, roughness-sublayer mixing over
tall canopies) and is not captured by surface-layer similarity. Without it,
the bulk exchange collapses (`u_* → 0`, `ζ → ζ_max = 100`) once the bulk
Richardson number becomes supercritical, which decouples tall canopies and
snow from the atmosphere at night ("runaway cooling"); over forests, where
the forcing height is typically only a few canopy roughness lengths above
the displacement height, this happens at wind speeds of 1–2 m/s and
surface–air temperature differences of a few kelvin.

Implementation: the dimensionless profiles `F_m`, `F_h` (and hence `u_*`
and the heat/moisture conductances) are evaluated at `min(ζ, ζ_max_stable)`,
which gives the same fluxes as CLM's `ζ = min(ζ, zetamaxstable)`. For
`ζ > ζ_max_stable`, the transfer coefficients become independent of
stability and the bulk Richardson number `Ri_b(ζ) = ζ F_h/F_m²` grows
linearly in `ζ`, so the Monin-Obukhov solve always has a root. Unlike CLM,
the `ζ` and `L_MO` reported by SurfaceFluxes are the unlimited values
consistent with the computed fluxes. Unstable conditions are unaffected;
`ζ_max_stable = Inf` recovers the unlimited universal functions.
"""
struct StableLimitedUniversalFunctionParams{
    FT,
    P <: UF.AbstractUniversalFunctionParameters{FT},
} <: UF.AbstractUniversalFunctionParameters{FT}
    base::P
    ζ_max_stable::FT
end

const SLUFP = StableLimitedUniversalFunctionParams
for f in (:phi, :psi, :Psi), tt in (:MomentumTransport, :HeatTransport)
    @eval @inline UF.$f(p::SLUFP, ζ, t::UF.$tt) = UF.$f(p.base, ζ, t)
end
for f in (:Pr_0, :a_m, :a_h, :b_m, :b_h, :c_h, :ζ_a, :γ)
    @eval UF.$f(p::SLUFP) = UF.$f(p.base)
end
@inline limit_stable_ζ(p::SLUFP, ζ) = min(ζ, oftype(ζ, p.ζ_max_stable))
@inline UF.dimensionless_profile(
    p::SLUFP,
    Δz_eff,
    ζ,
    z0,
    transport,
    scheme::UF.PointValueScheme,
) = UF.dimensionless_profile(
    p.base,
    Δz_eff,
    limit_stable_ζ(p, ζ),
    z0,
    transport,
    scheme,
)
@inline UF.dimensionless_profile(
    p::SLUFP,
    Δz_eff,
    ζ,
    z0,
    transport,
    scheme::UF.LayerAverageScheme,
) = UF.dimensionless_profile(
    p.base,
    Δz_eff,
    limit_stable_ζ(p, ζ),
    z0,
    transport,
    scheme,
)


# interfacing with ClimaParams
"""
    LandParameters(toml_dict::CP.ParamDict)

A constructor from `toml_dict` for the ClimaLand `earth_param_set`
(LandParameters) struct which contains the default values defined in ClimaParams
with type FT (Float32, Float64)

See [`ClimaLand.Parameters.create_toml_dict`](@ref).
"""
function LandParameters(toml_dict::CP.ParamDict)
    thermo_params = ThermodynamicsParameters(toml_dict)
    TP = typeof(thermo_params)

    insol_params = InsolationParameters(toml_dict)
    IP = typeof(insol_params)

    surf_flux_params_base =
        SurfaceFluxesParameters(toml_dict, UF.GryanikParams)
    (; zeta_max_stable) = CP.get_parameter_values(
        toml_dict,
        "zeta_max_stable",
        "Land",
    )
    surf_flux_params = SurfaceFluxesParameters(;
        von_karman_const = surf_flux_params_base.von_karman_const,
        ufp = StableLimitedUniversalFunctionParams(
            surf_flux_params_base.ufp,
            zeta_max_stable,
        ),
        thermo_params = surf_flux_params_base.thermo_params,
        z0m_fixed = surf_flux_params_base.z0m_fixed,
        z0s_fixed = surf_flux_params_base.z0s_fixed,
        gustiness_coeff = surf_flux_params_base.gustiness_coeff,
        gustiness_zi = surf_flux_params_base.gustiness_zi,
    )
    SFP = typeof(surf_flux_params)

    name_map = (;
        :light_speed => :light_speed,
        :planck_constant => :h_Planck,
        :density_ice_water => :ρ_cloud_ice,
        :avogadro_constant => :avogad,
        :thermodynamics_temperature_reference => :T_0,
        :temperature_water_freeze => :T_freeze,
        :density_liquid_water => :ρ_cloud_liq,
        :isobaric_specific_heat_ice => :cp_i,
        :latent_heat_sublimation_at_reference => :LH_s0,
        :molar_mass_water => :molmass_water,
        :mean_sea_level_pressure => :MSLP,
        :diffusivity_of_water_vapor => :D_vapor,
        :isobaric_specific_heat_liquid => :cp_l,
        :latent_heat_vaporization_at_reference => :LH_v0,
        :universal_gas_constant => :gas_constant,
        :thermal_conductivity_of_air => :K_therm,
        :gravitational_acceleration => :grav,
        :stefan_boltzmann_constant => :Stefan,
    )

    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    return LandParameters{FT, TP, SFP, IP}(;
        parameters...,
        thermo_params,
        surf_flux_params,
        insol_params,
    )
end

"""
    create_toml_dict(FT; override_files = [])

Construct a `ParamDict{FT}` struct from the default parameters with any
parameter overrides specified in `override_files`
"""
function create_toml_dict(FT; override_files = [])
    all(filepath -> endswith(filepath, ".toml"), override_files) ||
        error("File paths ($override_files) must be TOML files")
    toml_dict = CP.create_toml_dict(
        FT,
        override_file = CP.merge_toml_files(
            [DEFAULT_PARAMS_FILEPATH, override_files...],
            override = true,
        ),
    )
    return toml_dict
end

const DEFAULT_PARAMS_FILEPATH =
    joinpath(pkgdir(Parameters), "toml", "default_parameters.toml")

end # module
