module CreateParametersExt

import Thermodynamics.Parameters.ThermodynamicsParameters
import Insolation.Parameters.InsolationParameters
import SurfaceFluxes.Parameters.SurfaceFluxesParameters
import SurfaceFluxes.UniversalFunctions as UF
import CLIMAParameters as CP

import ClimaLand
import ClimaLand.Parameters.LandParameters
import ClimaLand.Canopy.AutotrophicRespirationParameters
import ClimaLand.Canopy.FarquharParameters
import ClimaLand.Canopy.OptimalityFarquharParameters

LandParameters(::Type{FT}) where {FT <: AbstractFloat} =
    LandParameters(CP.create_toml_dict(FT))

function LandParameters(toml_dict::CP.AbstractTOMLDict)
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

    parameters = CP.get_parameter_values(toml_dict, name_map, "ClimaLand")
    FT = CP.float_type(toml_dict)
    return LandParameters{FT, TP, SFP, IP}(;
        parameters...,
        thermo_params,
        surf_flux_params,
        insol_params,
    )
end

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
        :mol_CO2_to_kg_C_factor => :f1,
        :relative_contribution_factor => :f2,
        :stem_leaf_nitrogen_ratio => :μs,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    AutotrophicRespirationParameters{FT}(; parameters..., kwargs...)
end

"""
    function FarquharParameters(
        FT, 
        mechanism::AbstractPhotosynthesisMechanism;
        Vcmax25 = FT(5e-5),
        kwargs...  # For individual parameter overrides
    )

    function FarquharParameters(
        toml_dict::CP.AbstractTOMLDict, 
        mechanism::AbstractPhotosynthesisMechanism;
        Vcmax25 = FT(5e-5),
        kwargs...  # For individual parameter overrides
    )
    
Constructors for the FarquharParameters struct. Two variants:
1. Pass in the float-type and retrieve parameter values from the default TOML dict.
2. Pass in a TOML dictionary to retrieve parameter values.Possible calls:
```julia
ClimaLand.Canopy.FarquharParameters(Float64, ClimaLand.Canopy.C3())
# Kwarg overrides
ClimaLand.Canopy.FarquharParameters(Float64, ClimaLand.Canopy.C3(); Vcmax25 = 99999999, pc = 444444444)
# TOML Dictionary:
import CLIMAParameters as CP
toml_dict = CP.create_toml_dict(Float32);
ClimaLand.Canopy.FarquharParameters(toml_dict, ClimaLand.Canopy.C3(); Vcmax25 = 99999999, pc = 444444444)
```
"""
FarquharParameters(
    ::Type{FT},
    mechanism;
    kwargs...,
) where {FT <: AbstractFloat} =
    FarquharParameters(CP.create_toml_dict(FT), mechanism; kwargs...)

function FarquharParameters(
    toml_dict::CP.AbstractTOMLDict,
    mechanism;
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
    MECH = typeof(mechanism)
    return FarquharParameters{FT, MECH}(;
        mechanism,
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
import CLIMAParameters as CP
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
    mechanism = ClimaLand.Canopy.C3()
    return OptimalityFarquharParameters{FT}(; params..., kwargs..., mechanism)
end

end
