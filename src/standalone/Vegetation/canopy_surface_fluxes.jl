export MoninObukhovCanopyFluxes, MoninObukhovHeightBased
abstract type AbstractCanopySFParameterization{FT <: AbstractFloat} end

"""
    MoninObukhovCanopyFluxes{FT, F <: Union{FT, ClimaCore.Fields.Field}} <: AbstractCanopySFParameterization{FT}

A parameterization specifying how to compute latent and sensible heat
fluxes between the atmosphere and the canopy based on Monin-Obukhov
Surface Theory.

You must specify
- a minimum roughness length (global constant)
- the roughness length for momentum (can be a constant or a field)
- the roughness length for scalars (can be a constant or a field)
- the displacement height (can be a constant or a field)
- the parameter `η` which specifies how to upscale from leaf level
  to canopy fluxes via the leaf boundary layer resistance.
"""
struct MoninObukhovCanopyFluxes{FT, F <: Union{FT, ClimaCore.Fields.Field}} <:
       AbstractCanopySFParameterization{FT}
    "Minimum roughness length (m)"
    z_0min::FT
    "Canopy roughness length for momentum (m)"
    z_0m::F
    "Canopy roughness length for scalars (m)"
    z_0b::F
    "Canopy displacement height (m)"
    displ::F
    "Constant used to compute the resistance for fluxes between the leaves and canopy airspace (unitless)"
    η::FT
end

"""
    MoninObukhovHeightBased(toml_dict, height)

A constructor for a MoninObukhovCanopyFluxes surface flux theory,
 specifying how to compute latent and sensible heat
fluxes between the atmosphere and the canopy based on Monin-Obukhov
Surface Theory.

In this specific constructor,
the roughness lengths and displacement heights are assumed to be
linear in canopy height:
z_0m = coeff1 * height + z_0min
z_0b = coeff2 * height + z_0min
displacement = coeff3*height

where the coefficients are read from the toml_dict.

Cowan 1968; Brutsaert 1982, pp. 113–116; Campbell and Norman 1998, p. 71; Shuttleworth 2012, p. 343; Monteith and Unsworth 2013, p. 304
"""
function MoninObukhovHeightBased(toml_dict, height)
    z_0min = toml_dict["canopy_z_0min"]
    z_0m = toml_dict["canopy_z_0m_coeff"] .* height .+ z_0min
    z_0b = toml_dict["canopy_z_0b_coeff"] .* height .+ z_0min
    displ = toml_dict["canopy_d_coeff"] .* height
    η = toml_dict["canopy_η"]
    FT = typeof(η)
    F = typeof(height)
    return MoninObukhovCanopyFluxes{FT, F}(z_0min, z_0m, z_0b, displ, η)
end
