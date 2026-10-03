export MoninObukhovCanopyFluxes, subcanopy_wind
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
- the extinction coefficient of the wind speed below the canopy (unitless)
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
    "Extinction coefficient of the wind speed below the canopy per unit plant area index (unitless)"
    subcanopy_wind_extinction::FT
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
either a float or a field. The leaf drag coefficient and the extinction
coefficient of the wind speed below the canopy are also read from the
toml_dict.

Cowan 1968; Brutsaert 1982, pp. 113–116; Campbell and Norman 1998, p. 71; Shuttleworth 2012, p. 343; Monteith and Unsworth 2013, p. 304
"""
function MoninObukhovCanopyFluxes(toml_dict, height)
    z_0min = toml_dict["canopy_z_0min"]
    z_0m = toml_dict["canopy_z_0m_coeff"] .* height .+ z_0min
    z_0b = toml_dict["canopy_z_0b_coeff"] .* height .+ z_0min
    displ = toml_dict["canopy_d_coeff"] .* height
    Cd = toml_dict["leaf_Cd"]
    subcanopy_wind_extinction =
        toml_dict["canopy_subcanopy_wind_extinction_coefficient"]
    FT = typeof(Cd)
    F = typeof(height)
    return MoninObukhovCanopyFluxes{FT, F}(
        z_0min,
        z_0m,
        z_0b,
        displ,
        Cd,
        subcanopy_wind_extinction,
    )
end

"""
    subcanopy_wind(u, gustiness, PAI, extinction)

Return the wind speed at the ground below a canopy of plant area index
`PAI` (leaf plus stem area index, unitless) given the wind `u` and
gustiness `gustiness` at the atmospheric reference height [m/s]:

    u_ground = exp(-extinction * PAI) * max(u, gustiness).

Momentum absorption by foliage and stems attenuates the wind
exponentially with plant area index (Brutsaert, 1982; Mahfouf and
Noilhan, 1991; Norman et al., 1995). The gustiness is folded into the
wind before attenuation so that it acts as a floor on the wind speed
above the canopy, as in the Monin-Obukhov solve of SurfaceFluxes.jl
(`SurfaceFluxes.windspeed`), rather than below it; the ground-level
solve must therefore be given zero gustiness. For `PAI = 0` the result
is the effective wind speed `max(u, gustiness)` of that solve.

When `u` is a two-component `SVector`, the attenuated wind keeps the
direction of `u`; a zero vector is mapped onto the first component.

Called from `subcanopy_wind(canopy::CanopyModel, p)` in
`canopy_boundary_fluxes.jl`.
"""
function subcanopy_wind(u, gustiness, PAI, extinction)
    return exp(-extinction * PAI) * max(u, gustiness)
end

function subcanopy_wind(u::StaticArrays.SVector{2}, gustiness, PAI, extinction)
    speed = hypot(u[1], u[2])
    ground_speed = subcanopy_wind(speed, gustiness, PAI, extinction)
    return ifelse(
        speed > 0,
        (ground_speed / speed) * u,
        typeof(u)(ground_speed, zero(ground_speed)),
    )
end
