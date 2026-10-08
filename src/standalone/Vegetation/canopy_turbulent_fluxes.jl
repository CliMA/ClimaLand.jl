export MoninObukhovCanopyFluxes,
    canopy_z_0m,
    canopy_z_0b,
    canopy_displacement,
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
    "Excess resistance to scalar transfer kB⁻¹ = ln(z_0m / z_0b) (unitless)"
    kB_inv::FT
    "Frontal area index per unit plant area index of the Raupach (1994) roughness; 0 selects the static roughness `z_0m`, `displ` (unitless)"
    raupach_frontal_area::FT
    "Canopy height (m)"
    height::F
end

"""
    MoninObukhovCanopyFluxes(toml_dict, height)

A constructor for a MoninObukhovCanopyFluxes surface flux theory,
 specifying how to compute vapor fluxes, latent and sensible heat
fluxes, and momentum fluxes between the atmosphere and the canopy based on Monin-Obukhov Surface Theory, assuming that the roughness lengths
and displacment height are linear in the canopy height:
z_0m = coeff1 * height + z_0min
z_0b = z_0m * exp(-kB⁻¹)
displacement = coeff3*height

where the coefficients and kB⁻¹ = ln(z_0m / z_0b), the excess resistance
to scalar transfer of a rough vegetated surface (Garratt 1992, Ch. 4;
Brutsaert 1982, Ch. 5), are read from the toml_dict. The height can be
either a float or a field.

When `canopy_raupach_frontal_area` is positive, the momentum roughness
length and displacement height used in the flux solves are instead the
Raupach (1994) functions of canopy height and plant area index
([`canopy_z_0m`](@ref), [`canopy_displacement`](@ref)); the static values
above are then only used to check the forcing height at initialization. The leaf drag coefficient, the extinction
coefficient of the wind speed within the canopy, and the minimum reference
height of the ground below the canopy are also read from the toml_dict.

Cowan 1968; Brutsaert 1982, pp. 113–116; Garratt 1992, pp. 92–94; Campbell and Norman 1998, p. 71; Shuttleworth 2012, p. 343; Monteith and Unsworth 2013, p. 304
"""
function MoninObukhovCanopyFluxes(toml_dict, height)
    z_0min = toml_dict["canopy_z_0min"]
    z_0m = toml_dict["canopy_z_0m_coeff"] .* height .+ z_0min
    z_0b = z_0m .* exp(-toml_dict["canopy_kB_inv"])
    displ = toml_dict["canopy_d_coeff"] .* height
    Cd = toml_dict["leaf_Cd"]
    subcanopy_wind_extinction =
        toml_dict["canopy_subcanopy_wind_extinction_coefficient"]
    subcanopy_min_reference_height =
        toml_dict["canopy_subcanopy_min_reference_height"]
    kB_inv = toml_dict["canopy_kB_inv"]
    raupach_frontal_area = toml_dict["canopy_raupach_frontal_area"]
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
        kB_inv,
        raupach_frontal_area,
        height,
    )
end

# Raupach (1994, Boundary-Layer Meteorol. 71, 211-216) constants
const RAUPACH_CS = 0.003    # substrate drag coefficient
const RAUPACH_CR = 0.3      # roughness element drag coefficient
const RAUPACH_USTAR_MAX = 0.3   # maximum of u_* / U_h
const RAUPACH_CD1 = 7.5     # displacement height constant
const RAUPACH_PSI_H = 0.193 # roughness sublayer influence function
const RAUPACH_KAPPA = 0.4

"""
    raupach_displacement_fraction(Λ)
    raupach_roughness_fraction(Λ)

Displacement height and momentum roughness length of a canopy of frontal
area index `Λ`, as fractions of the canopy height (Raupach 1994, eqs. 8
and 9):

    d / h    = 1 - (1 - exp(-√(c_d1 Λ))) / √(c_d1 Λ)
    z_0m / h = (1 - d / h) exp(-κ U_h / u_* + ψ_h),  u_* / U_h = min(√(C_S + C_R Λ), 0.3)

`Λ` is floored at 0.05 so that leafless canopies keep the roughness of
their stems. Temporary workaround: this parameterization belongs in
SurfaceFluxes.jl and will migrate there.
"""
function raupach_displacement_fraction(Λ::FT) where {FT}
    x = sqrt(FT(RAUPACH_CD1) * max(Λ, FT(0.05)))
    return 1 - (1 - exp(-x)) / x
end

function raupach_roughness_fraction(Λ::FT) where {FT}
    Λ = max(Λ, FT(0.05))
    u_ratio =
        min(sqrt(FT(RAUPACH_CS) + FT(RAUPACH_CR) * Λ), FT(RAUPACH_USTAR_MAX))
    return (1 - raupach_displacement_fraction(Λ)) *
           exp(-FT(RAUPACH_KAPPA) / u_ratio + FT(RAUPACH_PSI_H))
end

"""
    canopy_z_0m(raupach_frontal_area, z_0min, z_0m, height, PAI)
    canopy_displacement(raupach_frontal_area, displ, height, PAI)
    canopy_z_0b(kB_inv, z_0m)

Momentum roughness length, displacement height and scalar roughness length
(m) of the canopy at a point, given the static values `z_0m`, `displ` of the
[`MoninObukhovCanopyFluxes`](@ref) parameterization, the canopy `height` and
the plant area index `PAI`. With `raupach_frontal_area > 0` the Raupach
(1994) values for the frontal area index `Λ = raupach_frontal_area * PAI`
are used, floored at `z_0min`; otherwise the static values.
`z_0b = z_0m exp(-kB⁻¹)`. The arguments are scalars so that the functions
broadcast over fields.
"""
function canopy_z_0m(
    raupach_frontal_area::FT,
    z_0min,
    z_0m,
    height,
    PAI,
) where {FT}
    return ifelse(
        raupach_frontal_area > 0,
        max(
            raupach_roughness_fraction(raupach_frontal_area * PAI) * height,
            z_0min,
        ),
        z_0m,
    )
end

function canopy_displacement(
    raupach_frontal_area::FT,
    displ,
    height,
    PAI,
) where {FT}
    return ifelse(
        raupach_frontal_area > 0,
        raupach_displacement_fraction(raupach_frontal_area * PAI) * height,
        displ,
    )
end

canopy_z_0b(kB_inv, z_0m) = z_0m * exp(-kB_inv)

"""
    subcanopy_reference_height(h_sfc, displ, z_0m, z_min, h_atmos)

Return the height [m] at which the ground below a canopy exchanges heat and
water vapor with the air: the apparent sink height `displ + z_0m` of the
canopy, which is where the wind profile above the canopy places the canopy
air, but at least `z_min` above the ground surface at height `h_sfc` and no
higher than the forcing height `h_atmos`. The floor matters for canopies
shorter than about `z_min` and for deep snow, which raises `h_sfc`.

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
the Boolean field `plants` (`LAI + SAI >= 0.05`): where plants are present,
the model of the atmospheric forcing `gustiness` with its minimum wind speed
set to zero, since that floor is folded into the sub-canopy wind
([`subcanopy_wind`](@ref)), so that only its convective part remains; where
they are absent, the model of the forcing itself, so that bare ground is
forced as a standalone surface. A number is a constant minimum wind speed.
Models without a floor are returned as they are.
"""
ground_gustiness(gustiness::Number, plants) =
    ground_gustiness(SurfaceFluxes.ConstantGustinessSpec(gustiness), plants)
function ground_gustiness(
    gustiness::SurfaceFluxes.ConstantGustinessSpec,
    plants,
)
    u_min = gustiness.value
    return @. lazy(
        SurfaceFluxes.ConstantGustinessSpec(ifelse(plants, zero(u_min), u_min)),
    )
end
function ground_gustiness(
    gustiness::SurfaceFluxes.FlooredDeardorffGustinessSpec,
    plants,
)
    u_min = gustiness.u_min
    return @. lazy(
        SurfaceFluxes.FlooredDeardorffGustinessSpec(
            ifelse(plants, zero(u_min), u_min),
        ),
    )
end
ground_gustiness(gustiness::SurfaceFluxes.AbstractGustinessSpec, plants) =
    gustiness
