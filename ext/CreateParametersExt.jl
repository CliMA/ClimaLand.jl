module CreateParametersExt

import Thermodynamics.Parameters.ThermodynamicsParameters
import Insolation.Parameters.InsolationParameters
import SurfaceFluxes.Parameters.SurfaceFluxesParameters
import SurfaceFluxes.UniversalFunctions as UF
import CLIMAParameters as CP

# Parameter structs to override
import ClimaLand.Parameters.LandParameters
import ClimaLand.Canopy.AutotrophicRespirationParameters

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

end
