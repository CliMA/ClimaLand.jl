module CreateParametersExt

import ClimaCore
import Thermodynamics.Parameters.ThermodynamicsParameters
import Insolation.Parameters.InsolationParameters
import SurfaceFluxes.Parameters.SurfaceFluxesParameters
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP

import ClimaLand
import ClimaLand.Soil
# Parameter structs
import ClimaLand.Parameters as LP
import ClimaLand.Soil.EnergyHydrologyParameters
import ClimaLand.Canopy.AutotrophicRespirationParameters
import ClimaLand.Canopy.FarquharParameters
import ClimaLand.Canopy.OptimalityFarquharParameters
import ClimaLand.Canopy.MedlynConductanceParameters
import ClimaLand.Canopy.BeerLambertParameters
import ClimaLand.Canopy.TwoStreamParameters
import ClimaLand.Canopy.ConstantGFunction
import ClimaLand.Snow.SnowParameters
import ClimaLand.Bucket.BucketModelParameters
import ClimaLand.Soil.Biogeochemistry.SoilCO2ModelParameters

LP.LandParameters(::Type{FT}) where {FT <: AbstractFloat} =
    LP.LandParameters(CP.create_toml_dict(FT))

function LP.LandParameters(toml_dict::CP.AbstractTOMLDict)
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
        :gas_constant => :gas_constant,
        :thermal_conductivity_of_air => :K_therm,
        :gravitational_acceleration => :grav,
        :stefan_boltzmann_constant => :Stefan,
    )

    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    return LP.LandParameters{FT, TP, SFP, IP}(;
        parameters...,
        thermo_params,
        surf_flux_params,
        insol_params,
    )
end
Base.broadcastable(ps::LP.LandParameters) = tuple(ps)

"""
    AutotrophicRespirationParameters(FT; kwargs...)
    AutotrophicRespirationParameters(toml_dict; kwargs...)

Constructors for the AutotrophicRespirationParameters struct. Two variants:
1. Pass in the float-type and retrieve parameter values from the default TOML dict.
2. Pass in a TOML dictionary to retrieve parameter values.
With either constructor, you can manually override any parameter via kwargs:
```julia
AutotrophicRespirationParameters(FT; ne = 99999)
AutotrophicRespirationParameters(toml_dict; ne = 99999)
```
"""
AutotrophicRespirationParameters(
    ::Type{FT};
    kwargs...,
) where {FT <: AbstractFloat} =
    AutotrophicRespirationParameters(CP.create_toml_dict(FT); kwargs...)

function AutotrophicRespirationParameters(
    toml_dict::CP.AbstractTOMLDict;
    kwargs...,
)
    name_map = (;
        :N_factor_Vcmax25 => :ne,
        :live_stem_wood_coeff => :ηsl,
        :specific_leaf_density => :σl,
        :root_leaf_nitrogen_ratio => :μr,
        :relative_contribution_factor => :Rel,
        :stem_leaf_nitrogen_ratio => :μs,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    AutotrophicRespirationParameters{FT}(; parameters..., kwargs...)
end

"""
    function FarquharParameters(
        FT,
        is_c3::Union{FT, ClimaCore.Fields.Field};
        Vcmax25 = FT(5e-5),
        kwargs...  # For individual parameter overrides
    )

    function FarquharParameters(
        toml_dict::CP.AbstractTOMLDict,
        is_c3::Union{AbstractFloat, ClimaCore.Fields.Field};
        Vcmax25 = FT(5e-5),
        kwargs...  # For individual parameter overrides
    )

Constructors for the FarquharParameters struct. Two variants:
1. Pass in the float-type and retrieve parameter values from the default TOML dict.
2. Pass in a TOML dictionary to retrieve parameter values.Possible calls:
```julia
ClimaLand.Canopy.FarquharParameters(Float64, 1.0)
# Kwarg overrides
ClimaLand.Canopy.FarquharParameters(Float64, 1.0; Vcmax25 = 99999999, pc = 444444444)
# TOML Dictionary:
import ClimaParams as CP
toml_dict = CP.create_toml_dict(Float32);
ClimaLand.Canopy.FarquharParameters(toml_dict, 1.0f0; Vcmax25 = 99999999, pc = 444444444)
```
"""
FarquharParameters(
    ::Type{FT},
    is_c3::Union{FT, ClimaCore.Fields.Field};
    kwargs...,
) where {FT <: AbstractFloat} =
    FarquharParameters(CP.create_toml_dict(FT), is_c3; kwargs...)

function FarquharParameters(
    toml_dict::CP.AbstractTOMLDict,
    is_c3::Union{AbstractFloat, ClimaCore.Fields.Field};
    Vcmax25 = 5e-5,
    kwargs...,
)
    name_map = (;
        :Jmax_activation_energy => :ΔHJmax,
        :intercellular_O2_concentration => :oi,
        :CO2_compensation_point_25c => :Γstar25,
        :Farquhar_curvature_parameter => :θj,
        :kelvin_25C => :To,
        :photosystem_II_quantum_yield => :ϕ,
        :O2_michaelis_menten => :Ko25,
        :CO2_michaelis_menten => :Kc25,
        :dark_respiration_factor => :f,
        :O2_activation_energy => :ΔHko,
        :low_water_pressure_sensitivity => :sc,
        :Rd_activation_energy => :ΔHRd,
        :Vcmax_activation_energy => :ΔHVcmax,
        :Γstar_activation_energy => :ΔHΓstar,
        :CO2_activation_energy => :ΔHkc,
        :moisture_stress_ref_water_pressure => :pc,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    if maximum(is_c3) > 1
        error(
            "is_c3 has maximum of $(maximum(is_c3)). is_c3 should be between 0 and 1",
        )
    end
    if minimum(is_c3) < 0
        error(
            "is_c3 has minimum of $(minimum(is_c3)). is_c3 should be between 0 and 1",
        )
    end
    # if is_c3 is a field, is_c3 may contain values between 0.0 and 1.0 after regridding
    # this deals with that possibility by rounding to the closest int
    is_c3 = round.(is_c3)
    MECH = typeof(is_c3)
    Vcmax25 = FT.(Vcmax25)
    VC = typeof(Vcmax25)
    return FarquharParameters{FT, MECH, VC}(;
        is_c3,
        Vcmax25,
        parameters...,
        kwargs...,
    )
end

"""
    function OptimalityFarquharParameters(
        FT,
        kwargs...  # For individual parameter overrides
    )

    function OptimalityFarquharParameters(
        toml_dict::CP.AbstractTOMLDict,
        kwargs...  # For individual parameter overrides
    )

Constructors for the OptimalityFarquharParameters struct. Two variants:
1. Pass in the float-type and retrieve parameter values from the default TOML dict.
2. Pass in a TOML dictionary to retrieve parameter values.Possible calls:
```julia
ClimaLand.Canopy.OptimalityFarquharParameters(Float64)
# Kwarg overrides
ClimaLand.Canopy.OptimalityFarquharParameters(Float64; pc = 444444444)
# Toml Dictionary:
import ClimaParams as CP
toml_dict = CP.create_toml_dict(Float32);
ClimaLand.Canopy.OptimalityFarquharParameters(toml_dict; pc = 444444444)
```
"""
OptimalityFarquharParameters(
    ::Type{FT};
    kwargs...,
) where {FT <: AbstractFloat} =
    OptimalityFarquharParameters(CP.create_toml_dict(FT); kwargs...)

function OptimalityFarquharParameters(toml_dict; kwargs...)
    name_map = (;
        :Jmax_activation_energy => :ΔHJmax,
        :intercellular_O2_concentration => :oi,
        :CO2_compensation_point_25c => :Γstar25,
        :Farquhar_curvature_parameter => :θj,
        :kelvin_25C => :To,
        :photosystem_II_quantum_yield => :ϕ,
        :O2_michaelis_menten => :Ko25,
        :CO2_michaelis_menten => :Kc25,
        :dark_respiration_factor => :f,
        :O2_activation_energy => :ΔHko,
        :low_water_pressure_sensitivity => :sc,
        :Rd_activation_energy => :ΔHRd,
        :Vcmax_activation_energy => :ΔHVcmax,
        :electron_transport_maintenance => :c,
        :Γstar_activation_energy => :ΔHΓstar,
        :CO2_activation_energy => :ΔHkc,
        :moisture_stress_ref_water_pressure => :pc,
    )

    params = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    is_c3 = FT(1)
    return OptimalityFarquharParameters{FT, FT}(; params..., kwargs..., is_c3)
end


"""
    EnergyHydrologyParameters(
        ::Type{FT};
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        kwargs...,)

    EnergyHydrologyParameters(
        toml_dict;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        kwargs...,)

EnergyHydrologyParameters has two constructors: float-type and toml dict based.
Additional parameters must be added manually: ν, ν_ss_om, ν_ss_quartz, ν_ss_gravel, hydrology_cm, K_sat, S_s, and θ_r
All parameters can be manually overriden via keyword arguments. Note, however,
that certain parameters must have the same type (e.g, if a field is
supplied for porosity, it must be supplied for all other parameters
defined in the interior of the domain). Some parameters are defined only
on the surface of the domain (e.g albedo), while other are defined everywhere
(e.g. porosity). These are indicated with types `F` and `SF`.

Please see the EnergyHydrologyParameters documentation for a complete list.
"""
function EnergyHydrologyParameters(
    ::Type{FT};
    ν,
    ν_ss_om,
    ν_ss_quartz,
    ν_ss_gravel,
    hydrology_cm,
    K_sat,
    S_s,
    θ_r,
    kwargs...,
) where {FT <: AbstractFloat}
    return EnergyHydrologyParameters(
        CP.create_toml_dict(FT);
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        kwargs...,
    )
end

function EnergyHydrologyParameters(
    toml_dict::CP.AbstractTOMLDict;
    ν::F,
    ν_ss_om::F,
    ν_ss_quartz::F,
    ν_ss_gravel::F,
    hydrology_cm::C,
    K_sat::F,
    S_s::F,
    θ_r::F,
    PAR_albedo::SF = 0.2,
    NIR_albedo::SF = 0.4,
    kwargs...,
) where {
    F <: Union{<:AbstractFloat, ClimaCore.Fields.Field},
    SF <: Union{<:AbstractFloat, ClimaCore.Fields.Field},
    C,
}
    earth_param_set = LP.LandParameters(toml_dict)

    # Obtain parameters needed to calculate the derived parameters
    derived_param_name_map = (;
        :thermal_conductivity_of_quartz => :κ_quartz,
        :thermal_conductivity_of_soil_minerals => :κ_minerals,
        :thermal_conductivity_of_water_ice => :κ_ice,
        :thermal_conductivity_of_liquid_water => :κ_liq,
        :thermal_conductivity_of_organic_matter => :κ_om,
        :particle_density_quartz => :ρp_quartz,
        :particle_density_minerals => :ρp_minerals,
        :particle_density_organic_matter => :ρp_om,
        :vol_heat_capacity_quartz => :ρc_quartz,
        :vol_heat_capacity_organic_matter => :ρc_om,
        :vol_heat_capacity_minerals => :ρc_minerals,
    )
    p = CP.get_parameter_values(toml_dict, derived_param_name_map, "Land")
    # Particle density of the soil - per unit soil solids
    # Denoted ρ_ds in the Clima Design Docs (Equation 2.3)
    # where ν_ss_i = ν_i/(1-ν)
    ρp = @. (
        ν_ss_om * p.ρp_om +
        ν_ss_quartz * p.ρp_quartz +
        ν_ss_gravel * p.ρp_minerals +
        (1 - ν_ss_om - ν_ss_quartz - ν_ss_gravel) * p.ρp_minerals
    )
    # Volumetric heat capacity of soil solids - per unit volume soil solids
    # This is Equation 2.6a/2.10 in the Clima Design Docs
    # where ν_ss_i = ν_i/(1-ν)
    ρc_ss = @. (
        ν_ss_om * p.ρc_om +
        ν_ss_quartz * p.ρc_quartz +
        ν_ss_gravel * p.ρc_minerals +
        (1 - ν_ss_om - ν_ss_quartz - ν_ss_gravel) * p.ρc_minerals
    )
    # Volumetric heat capacity of dry soil - per unit volume of soil
    # This is denoted č_ds (Equation 2.11) and is used in equation 2.9
    ρc_ds = @. (1 - ν) * ρc_ss

    κ_solid =
        Soil.κ_solid.(ν_ss_om, ν_ss_quartz, p.κ_om, p.κ_quartz, p.κ_minerals)
    κ_dry = Soil.κ_dry.(ρp, ν, κ_solid, LP.K_therm(earth_param_set))
    κ_sat_frozen = Soil.κ_sat_frozen.(κ_solid, ν, p.κ_ice)
    κ_sat_unfrozen = Soil.κ_sat_unfrozen.(κ_solid, ν, p.κ_liq)

    name_map = (;
        :kersten_number_alpha => :α,
        :kersten_number_beta => :β,
        :ice_impedance_omega => :Ω,
        :temperature_factor_soil_hydraulic_conductivity => :γ,
        :temperature_reference_soil_hydraulic_conductivity => :γT_ref,
        :emissivity_bare_soil => :emissivity,
        :maximum_dry_soil_layer_depth => :d_ds,
        :soil_momentum_roughness_length => :z_0m,
        :soil_scalar_roughness_length => :z_0b,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    PSE = typeof(earth_param_set)
    FT = CP.float_type(toml_dict)
    EnergyHydrologyParameters{FT, F, SF, C, PSE}(;
        PAR_albedo,
        NIR_albedo,
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        κ_dry,
        κ_sat_frozen,
        κ_sat_unfrozen,
        ρc_ds,
        earth_param_set,
        parameters...,
        kwargs...,
    )
end


"""
    BucketModelParameters(
        ::Type{FT};
        albedo,
        z_0m,
        z_0b,
        τc,
        kwargs...,
    )

    BucketModelParameters(
        toml_dict::CP.AbstractTOMLDict;
        albedo,
        z_0m,
        z_0b,
        τc,
        kwargs...,
    )

BucketModelParameters has a float-type and a toml-dict based constructor.
Keyword arguments can be used to manually override any of the values in the struct.
```julia
BucketModelParameters(Float64; albedo, z_0m, z_0b, τc)
BucketModelParameters(toml_dict; albedo, z_0m, z_0b, τc)
```
"""
function BucketModelParameters(
    ::Type{FT};
    albedo,
    z_0m,
    z_0b,
    τc,
    kwargs...,
) where {FT <: AbstractFloat}
    return BucketModelParameters(
        CP.create_toml_dict(FT);
        albedo,
        z_0m,
        z_0b,
        τc,
        kwargs...,
    )
end

function BucketModelParameters(
    toml_dict::CP.AbstractTOMLDict;
    albedo,
    z_0m,
    z_0b,
    τc,
    kwargs...,
)

    name_map = (;
        :soil_conductivity => :κ_soil,
        :soil_heat_capacity => :ρc_soil,
        :critical_snow_water_equivalent => :σS_c,
        :land_bucket_capacity => :W_f,
        :critical_snow_fraction => :f_snow,
        :bucket_capacity_fraction => :f_bucket,
        :bucket_beta_decay_exponent => :p,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")

    AAM = typeof(albedo)
    earth_param_set = LP.LandParameters(toml_dict)
    PSE = typeof(earth_param_set)
    FT = CP.float_type(toml_dict)
    BucketModelParameters{FT, AAM, PSE}(;
        albedo,
        z_0m,
        z_0b,
        τc,
        earth_param_set,
        parameters...,
        kwargs...,
    )
end

"""
    SoilCO2ModelParameters(FT; kwargs...)
    SoilCO2ModelParameters(toml_dict; kwargs...)

SoilCO2ModelParameters has two constructors: float-type and toml dict based.
Keywords arguments can be used to directly override any parameters.
"""
SoilCO2ModelParameters(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    SoilCO2ModelParameters(CP.create_toml_dict(FT); kwargs...)

function SoilCO2ModelParameters(toml_dict::CP.AbstractTOMLDict; kwargs...)
    name_map = (;
        :CO2_diffusion_coefficient => :D_ref,
        :soil_C_substrate_diffusivity => :D_liq,
        :soilCO2_pre_expontential_factor => :α_sx,
        :soilCO2_activation_energy => :Ea_sx,
        :michaelis_constant => :kM_sx,
        :O2_michaelis_constant => :kM_o2,
        :O2_volume_fraction => :O2_a,
        :oxygen_diffusion_coefficient => :D_oa,
        :soluble_soil_carbon_fraction => :p_sx,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    earth_param_set = LP.LandParameters(toml_dict)
    return SoilCO2ModelParameters{FT, typeof(earth_param_set)}(;
        earth_param_set,
        parameters...,
        kwargs...,
    )
end

"""
    function MedlynConductanceParameters(FT::AbstractFloat;
        g1 = 790,
        kwargs...
    )
    function MedlynConductanceParameters(toml_dict;
        g1 = 790,
        kwargs...
    )

Floating-point and toml dict based constructor supplying default values
for the MedlynConductanceParameters struct.
Additional parameter values can be directly set via kwargs.
"""
MedlynConductanceParameters(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    MedlynConductanceParameters(CP.create_toml_dict(FT); kwargs...)

function MedlynConductanceParameters(
    toml_dict::CP.AbstractTOMLDict;
    g1 = 790,
    kwargs...,
)
    name_map = (;
        :relative_diffusivity_of_water_vapor => :Drel,
        :min_stomatal_conductance => :g0,
    )

    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    g1 = FT.(g1)
    G1 = typeof(g1)
    return MedlynConductanceParameters{FT, G1}(; g1, parameters..., kwargs...)
end

"""
    function TwoStreamParameters(FT::AbstractFloat;
        ld = (_) -> 0.5,
        α_PAR_leaf = 0.3,
        τ_PAR_leaf = 0.2,
        α_NIR_leaf = 0.4,
        τ_NIR_leaf = 0.25,
        Ω = 1,
        n_layers = UInt64(20),
        kwargs...
    )
    function TwoStreamParameters(toml_dict;
        ld = (_) -> 0.5,
        α_PAR_leaf = 0.3,
        τ_PAR_leaf = 0.2,
        α_NIR_leaf = 0.4,
        τ_NIR_leaf = 0.25,
        Ω = 1,
        n_layers = UInt64(20),
        kwargs...
    )

Floating-point and toml dict based constructor supplying default values
for the TwoStreamParameters struct. Additional parameter values can be directly set via kwargs.
"""
TwoStreamParameters(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    TwoStreamParameters(CP.create_toml_dict(FT); kwargs...)

function TwoStreamParameters(
    toml_dict::CP.AbstractTOMLDict;
    G_Function = ConstantGFunction(CP.float_type(toml_dict)(0.5)),
    α_PAR_leaf::F = 0.3,
    τ_PAR_leaf::F = 0.2,
    α_NIR_leaf::F = 0.4,
    τ_NIR_leaf::F = 0.25,
    Ω = 1,
    n_layers = UInt64(20),
    kwargs...,
) where {F}
    name_map = (;
        :wavelength_per_PAR_photon => :λ_γ_PAR,
        :canopy_emissivity => :ϵ_canopy,
    )

    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    # default value for keyword args must be converted manually
    # automatic conversion not possible to Union types
    α_PAR_leaf = FT.(α_PAR_leaf)
    τ_PAR_leaf = FT.(τ_PAR_leaf)
    α_NIR_leaf = FT.(α_NIR_leaf)
    τ_NIR_leaf = FT.(τ_NIR_leaf)
    return TwoStreamParameters{FT, typeof(G_Function), typeof(α_PAR_leaf)}(;
        G_Function,
        α_PAR_leaf,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        Ω,
        n_layers,
        parameters...,
        kwargs...,
    )
end

"""
    function BeerLambertParameters(FT::AbstractFloat;
        ld = (_) -> 0.5,
        α_PAR_leaf = 0.1,
        α_NIR_leaf = 0.4,
        Ω = 1,
        kwargs...
    )
    function BeerLambertParameters(toml_dict;
        ld = (_) -> 0.5,
        α_PAR_leaf = 0.1,
        α_NIR_leaf = 0.4,
        Ω = 1,
        kwargs...
    )

Floating-point and toml dict based constructor supplying default values
for the BeerLambertParameters struct. Additional parameter values can be directly set via kwargs.
"""
BeerLambertParameters(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    BeerLambertParameters(CP.create_toml_dict(FT); kwargs...)

function BeerLambertParameters(
    toml_dict::CP.AbstractTOMLDict;
    G_Function = ConstantGFunction(CP.float_type(toml_dict)(0.5)),
    α_PAR_leaf::F = 0.1,
    α_NIR_leaf::F = 0.4,
    Ω = 1,
    kwargs...,
) where {F}
    name_map = (;
        :wavelength_per_PAR_photon => :λ_γ_PAR,
        :canopy_emissivity => :ϵ_canopy,
    )

    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    # default value for keyword args must be converted manually
    # automatic conversion not possible to Union types
    α_PAR_leaf = FT.(α_PAR_leaf)
    α_NIR_leaf = FT.(α_NIR_leaf)
    return BeerLambertParameters{FT, typeof(G_Function), typeof(α_PAR_leaf)}(;
        G_Function,
        α_PAR_leaf,
        α_NIR_leaf,
        Ω,
        parameters...,
        kwargs...,
    )
end


SnowParameters(::Type{FT}, Δt; kwargs...) where {FT <: AbstractFloat} =
    SnowParameters(CP.create_toml_dict(FT), Δt; kwargs...)

function SnowParameters(toml_dict::CP.AbstractTOMLDict, Δt; kwargs...)
    name_map = (;
        :snow_momentum_roughness_length => :z_0m,
        :snow_scalar_roughness_length => :z_0b,
        :thermal_conductivity_of_water_ice => :κ_ice,
        :snow_density => :ρ_snow,
        :snow_albedo => :α_snow,
        :snow_emissivity => :ϵ_snow,
        :holding_capacity_of_water_in_snow => :θ_r,
        :wet_snow_hydraulic_conductivity => :Ksat,
        :snow_cover_fraction_crit_threshold => :fS_c,
    )

    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    earth_param_set = LP.LandParameters(toml_dict)
    PSE = typeof(earth_param_set)
    return SnowParameters{FT, PSE}(;
        Δt,
        earth_param_set,
        parameters...,
        kwargs...,
    )
end

end
