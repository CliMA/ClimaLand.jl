abstract type AbstractCanopySFParameterization{FT <: AbstractFloat} end

struct MoninObukhovHeightBased{FT, F <: Union{FT, ClimaCore.Fields.Field}} <:
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

function MoninObukhovHeightBased{FT}(toml_dict, height) where {FT}
    z_0min = toml_dict["canopy_z_0min"]
    z_0m = toml_dict["canopy_z_0m_coeff"] .* height .+ z_0min
    z_0b = toml_dict["canopy_z_0b_coeff"] .* height .+ z_0min
    displ = toml_dict["canopy_d_coeff"] .* height
    η = toml_dict["canopy_η"]
    F = typeof(height)
    return MoninObukhovHeightBased{FT, F}(z_0min, z_0m, z_0b, displ, η)
end
