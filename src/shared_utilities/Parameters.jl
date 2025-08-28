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

    surf_flux_params = SurfaceFluxesParameters(toml_dict, UF.BusingerParams)
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
    get_default_parameter(FT, climaparams_name)

Helper function for accessing and returning the default
parameter value (a float of type `FT`) from  the ClimaParams
dictionary. To look up the value, you must use the name
of the parameter from ClimaParams (`climaparams_name`).
For example, the following keys in the ClimaParams dictionary map
to the variables in the parameter structs as noted below:

:soil_momentum_roughness_length => :z_0m (EnergyHydrologyParameters)
:soil_scalar_roughness_length => :z_0b (EnergyHydrologyParameters)
:emissivity_bare_soil => :emissivity (EnergyHydrologyParameters)

:snow_momentum_roughness_length => :z_0m (SnowParameters)
:snow_scalar_roughness_length => :z_0b (SnowParameters)
:snow_emissivity => :ϵ_snow (SnowParameters)
:holding_capacity_of_water_in_snow => :θ_r (SnowParameters)
:wet_snow_hydraulic_conductivity => :Ksat (SnowParameters)

:canopy_emissivity => :ϵ_canopy (TwoStreamParameters, BeerLambertParameters)
:min_stomatal_conductance => :g0 (MedlynConductanceParameters)
:low_water_pressure_sensitivity => :sc (FarquharParameters)
:moisture_stress_ref_water_pressure => :pc (FarquharParameters)

:soil_conductivity => :κ_soil (BucketModelParameters)
:soil_heat_capacity => :ρc_soil (BucketModelParameters)
:critical_snow_water_equivalent => :σS_c (BucketModelParameters)
:land_bucket_capacity => :W_f (BucketModelParameters)
:critical_snow_fraction => :f_snow (BucketModelParameters)
:bucket_capacity_fraction => :f_bucket (BucketModelParameters)
:bucket_beta_decay_exponent => :p (BucketModelParameters)
"""
function get_default_parameter(FT, climaparams_name)
    toml_dict = CP.create_toml_dict(FT)
    return CP.get_parameter_values(toml_dict, string(climaparams_name))[climaparams_name]
end

"""
    create_toml_dict(FT, filepaths...; override = false)

Construct a `ParamDict{FT}` struct from `filepaths`.

If `override = false`, then non-unique TOML entries are not allowed. If
`override = true`, then parameters from later TOML files in `filepaths` will
overwrite the parameters from earlier TOML files.
"""
function create_toml_dict(FT, filepaths...; override = false)
    all(filepath -> endswith(filepath, ".toml"), filepaths) ||
        error("File paths ($filepaths) must be TOML files")
    toml_dict = CP.create_toml_dict(
        FT,
        override_file = CP.merge_toml_files(
            [filepaths...],
            override = override,
        ),
    )
    return toml_dict
end

const DEFAULT_PARAMS_FILEPATH =
    joinpath(pkgdir(Parameters), "toml", "default_parameters.toml")

end # module
