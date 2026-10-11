export MoninObukhovCanopyFluxes,
    subcanopy_wind,
    subcanopy_reference_height,
    subcanopy_forcing,
    ground_gustiness
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
- the extinction coefficient of the wind speed within the canopy (unitless)
- the minimum height above the ground at which the fluxes of the ground
  below the canopy are evaluated (m)
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
    "Extinction coefficient of the wind speed within the canopy per unit plant area index (unitless)"
    subcanopy_wind_extinction::FT
    "Minimum height above the ground surface at which the fluxes of the ground below the canopy are evaluated (m)"
    subcanopy_min_reference_height::FT
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
either a float or a field. The leaf drag coefficient, the extinction
coefficient of the wind speed within the canopy, and the minimum reference
height of the ground below the canopy are also read from the toml_dict.

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
    subcanopy_min_reference_height =
        toml_dict["canopy_subcanopy_min_reference_height"]
    FT = typeof(Cd)
    F = typeof(height)
    return MoninObukhovCanopyFluxes{FT, F}(
        z_0min,
        z_0m,
        z_0b,
        displ,
        Cd,
        subcanopy_wind_extinction,
        subcanopy_min_reference_height,
    )
end

"""
    subcanopy_reference_height(h_sfc, displ, z_0m, z_min, h_atmos)

Return the height [m] at which the ground below a canopy exchanges heat and
water vapor with the air: the apparent sink height `displ + z_0m` of the
canopy, which is where the wind profile above the canopy places the canopy
air, but at least `z_min` above the ground surface at height `h_sfc` and no
higher than the forcing height `h_atmos`. The floor matters for canopies
shorter than about `z_min`.

The temperature and humidity at this height are taken to be the canopy-air
values (`p.canopy.turbulent_fluxes.T_sfc` and `q_sfc`) from the canopy
Monin-Obukhov solve, and the wind is attenuated across the full plant area
index of the canopy (see [`subcanopy_wind`](@ref) and
[`subcanopy_forcing`](@ref)).
"""
function subcanopy_reference_height(h_sfc, displ, z_0m, z_min, h_atmos)
    return min(max(displ + z_0m, h_sfc + z_min), h_atmos)
end

"""
    subcanopy_wind(u, gustiness, h_atmos, z, height, displ, z_0m, PAI, extinction)

Return the wind speed [m/s] driving the turbulent exchange of the ground below
a canopy of height `height`, displacement height `displ`, momentum roughness
length `z_0m`, and plant area index `PAI` (leaf plus stem area index,
unitless), evaluated with reference height `z` below the forcing height
`h_atmos`, given the wind `u` and gustiness `gustiness` at `h_atmos`:

    U = max(u, gustiness)
        * log((max(z, height) - displ) / z_0m) / log((h_atmos - displ) / z_0m)
        * exp(-extinction * PAI).

Above the canopy the wind follows the neutral logarithmic profile of the
canopy down to `max(z, height)`; across the canopy it is attenuated
exponentially by the full plant area index `PAI` between the canopy top and the
ground (Inoue, 1963; Cionco, 1965; Brutsaert, 1982; Shuttleworth and Wallace,
1985; Choudhury and Monteith, 1988), with the extinction coefficient
`extinction` per unit plant area index. Attenuating by the full `PAI` ensures
that dense short canopies (`height < z`, where `z` is floored at
`subcanopy_min_reference_height`) still shelter the ground beneath them. The
gustiness is folded into the wind at `h_atmos` so that it acts as a floor on
the wind speed above the canopy, as in the Monin-Obukhov solve of
SurfaceFluxes.jl (`SurfaceFluxes.windspeed`), rather than below it; the
ground-level solve is therefore given the gustiness model without its floor
(see [`ground_gustiness`](@ref)).

When `u` is a two-component `SVector`, the attenuated wind keeps the
direction of `u`; a zero vector is mapped onto the first component.

Called from [`subcanopy_forcing`](@ref).
"""
function subcanopy_wind(
    u,
    gustiness,
    h_atmos,
    z,
    height,
    displ,
    z_0m,
    PAI,
    extinction,
)
    log_profile =
        log((max(z, height) - displ) / z_0m) / log((h_atmos - displ) / z_0m)
    attenuation = exp(-extinction * PAI)
    return max(u, gustiness) * log_profile * attenuation
end

function subcanopy_wind(
    u::StaticArrays.SVector{2},
    gustiness,
    h_atmos,
    z,
    height,
    displ,
    z_0m,
    PAI,
    extinction,
)
    speed = hypot(u[1], u[2])
    ground_speed = subcanopy_wind(
        speed,
        gustiness,
        h_atmos,
        z,
        height,
        displ,
        z_0m,
        PAI,
        extinction,
    )
    return ifelse(
        speed > 0,
        (ground_speed / speed) * u,
        typeof(u)(ground_speed, zero(ground_speed)),
    )
end

"""
    ground_gustiness(gustiness, plants)

Return the gustiness model of the ground below a canopy as a lazy field over
the Boolean field `plants` (`plants_present`): where plants are present,
the model of the atmospheric forcing `gustiness` with its minimum wind speed
set to zero, since that floor is folded into the sub-canopy wind
([`subcanopy_wind`](@ref)), so that only its convective part remains; where
they are absent, the model of the forcing itself, so that bare ground is
forced as a standalone surface (see `SurfaceFluxes.without_floor`). Models
without a floor are returned as they are.
"""
function ground_gustiness(
    gustiness::Union{
        SurfaceFluxes.ConstantGustinessSpec,
        SurfaceFluxes.FlooredDeardorffGustinessSpec,
    },
    plants,
)
    unfloored = SurfaceFluxes.without_floor(gustiness)
    choose_gustiness(has_plants) = ifelse(has_plants, unfloored, gustiness)
    return @. lazy(choose_gustiness(plants))
end
ground_gustiness(gustiness::SurfaceFluxes.AbstractGustinessSpec, plants) =
    gustiness
