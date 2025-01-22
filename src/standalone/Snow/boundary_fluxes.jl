"""
    AbstractSnowBC <: ClimaLand.AbstractBC

An abstract type for boundary conditions for the snow model.
"""
abstract type AbstractSnowBC <: ClimaLand.AbstractBC end

"""
    AtmosDrivenSnowBC{
        A <: AbstractAtmosphericDrivers,
        B <: AbstractRadiativeDrivers,
        C::Tuple
    } <: AbstractSnowBC

A struct used to specify the snow fluxes, referred
to as ``boundary conditions", at the surface and
bottom of the snowpack, for water and energy.

These fluxes include turbulent surface fluxes
computed with Monin-Obukhov theory, and radiative fluxes.
$(DocStringExtensions.FIELDS)
"""
struct AtmosDrivenSnowBC{
    A <: AbstractAtmosphericDrivers,
    B <: AbstractRadiativeDrivers,
    C <: Tuple,
} <: AbstractSnowBC
    "The atmospheric conditions driving the model"
    atmos::A
    "The radiative fluxes driving the model"
    radiation::B
    "Prognostic land components present"
    prognostic_land_components::C
end

function AtmosDrivenSnowBC(
    atmos,
    radiation;
    prognostic_land_components = (:snow,),
)
    args = (atmos, radiation, prognostic_land_components)
    return AtmosDrivenSnowBC(args...)
end


"""
    snow_boundary_fluxes!(bc::AtmosDrivenSnowBC, model::SnowModel, Y, p, t)

Updates in place various volumetric water flux (m/s) and energy
flux (W/m^2) terms for the snow model:

- p.snow.turbulent fluxes (latent, sensible, and evaporative fluxes)
- p.snow.R_n (radiative fluxes)
- p.snow.total_water_flux
- p.snow.total_energy_flux

The two latter fluxes also include contributions from fluxes due to melt and
precipitation, but note that precipitation and melt flux are not computed or
updated in `snow_boundary_fluxes` currently. Instead, they are updated in
`update_aux!`, which happens prior to the `snow_boundary_fluxes!` call, and used
in the `snow_boundary_fluxes!` call.


This function calls the `turbulent_fluxes!` and `net_radiation!` functions,
which use the snow surface conditions as well as the atmos and radiation
conditions in order to compute the surface fluxes using Monin Obukhov Surface
Theory. It also accounts for the presence of other components, if run as part of
an integrated land model, and their effect on boundary conditions.
"""
function snow_boundary_fluxes!(bc::AtmosDrivenSnowBC, model::SnowModel, Y, p, t)
    snow_boundary_fluxes!(
        bc,
        Val(bc.prognostic_land_components),
        model,
        Y,
        p,
        t,
    )
    return nothing
end

"""
    snow_boundary_fluxes!(
        bc::AtmosDrivenSnowBC,
        prognostic_land_components::Val{(:snow,)},
        model::SnowModel{FT},
        Y,
        p,
        t,
    ) where {FT}

Computes the boundary fluxes for the snow model in standalone mode.

The ground heat flux is assumed to be zero, and the snow surface is assumed to
be bare (no vegetation).
"""
function snow_boundary_fluxes!(
    bc::AtmosDrivenSnowBC,
    prognostic_land_components::Val{(:snow,)},
    model::SnowModel{FT},
    Y,
    p,
    t,
) where {FT}
    bc = model.boundary_conditions

    turbulent_fluxes!(p.snow.turbulent_fluxes, bc.atmos, model, Y, p, t)
    net_radiation!(p.snow.R_n, bc.radiation, model, Y, p, t)

    # How does rain affect the below?
    P_snow = p.drivers.P_snow

    @. p.snow.total_water_flux =
        P_snow +
        (p.snow.turbulent_fluxes.vapor_flux - p.snow.water_runoff) *
        p.snow.snow_cover_fraction

    # I think we want dU/dt to include energy of falling snow.
    # otherwise snow can fall but energy won't change
    # We are assuming that the sensible heat portion of snow is negligible.
    _LH_f0 = FT(LP.LH_f0(model.parameters.earth_param_set))
    _ρ_liq = FT(LP.ρ_cloud_liq(model.parameters.earth_param_set))
    ρe_falling_snow = -_LH_f0 * _ρ_liq # per unit vol of liquid water

    # positive fluxes are TOWARDS atmos
    @. p.snow.total_energy_flux =
        P_snow * ρe_falling_snow +
        (
            p.snow.turbulent_fluxes.lhf +
            p.snow.turbulent_fluxes.shf +
            p.snow.R_n - p.snow.energy_runoff
        ) * p.snow.snow_cover_fraction
    return nothing
end
