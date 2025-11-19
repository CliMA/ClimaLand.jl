export MoninObukhovCanopyFluxes
abstract type AbstractCanopyFluxParameterization{FT <: AbstractFloat} end

"""
    MoninObukhovCanopyFluxes{FT, F <: Union{FT, ClimaCore.Fields.Field}} <: AbstractCanopyFluxParameterization{FT}

A parameterization specifying how to compute latent and sensible heat
fluxes between the atmosphere and the canopy based on Monin-Obukhov
Surface Theory.

You must specify
- a minimum roughness length (global constant)
- the roughness length for momentum (can be a constant or a field)
- the roughness length for scalars (can be a constant or a field)
- the displacement height (can be a constant or a field)
- the leaf-level drag coefficient (unitless)
"""
struct MoninObukhovCanopyFluxes{FT, F <: Union{FT, ClimaCore.Fields.Field}} <:
       AbstractCanopyFluxParameterization{FT}
    "Minimum roughness length (m)"
    z_0min::FT
    "Canopy roughness length for momentum (m)"
    z_0m::F
    "Canopy roughness length for scalars (m)"
    z_0b::F
    "Canopy displacement height (m)"
    displ::F
    "Leaf level drag coefficient (unitless)"
    Cd::FT
end

"""
    MoninObukhovCanopyFluxes(toml_dict, height)

A constructor for a MoninObukhovCanopyFluxes surface flux theory,
 specifying how to compute vapor fluxes, latent and sensible heat
fluxes, and momentum fluxes between the atmosphere and the canopy based on Monin-Obukhov Surface Theory, assuming that the roughness lengths
and displacment height are linear in the canopy height:
z_0m = coeff1 * height + z_0min
z_0b = coeff2 * height + z_0min
displacement = coeff3*height

where the coefficients are read from the toml_dict. The height can be 
either a float or a field.

Cowan 1968; Brutsaert 1982, pp. 113–116; Campbell and Norman 1998, p. 71; Shuttleworth 2012, p. 343; Monteith and Unsworth 2013, p. 304
"""
function MoninObukhovCanopyFluxes(toml_dict, height)
    z_0min = toml_dict["canopy_z_0min"]
    z_0m = toml_dict["canopy_z_0m_coeff"] .* height .+ z_0min
    z_0b = toml_dict["canopy_z_0b_coeff"] .* height .+ z_0min
    displ = toml_dict["canopy_d_coeff"] .* height
    Cd = toml_dict["leaf_Cd"]
    FT = typeof(Cd)
    F = typeof(height)
    return MoninObukhovCanopyFluxes{FT, F}(z_0min, z_0m, z_0b, displ, Cd)
end
