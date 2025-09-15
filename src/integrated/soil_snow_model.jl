export SoilSnowModel
using ClimaCore.Operators: column_integral_definite!


"""
    struct SoilSnowModel{
        FT,
        SnM <: Snow.SnowModel{FT},
        SoM <: Soil.EnergyHydrology{FT},
    } <: AbstractLandModel{FT}
        "The snow model to be used"
        snow::SnM
        "The soil model to be used"
        soil::SoM
    end

A concrete type of land model used for simulating systems with
snow and soil.

The inner constructor checks that the two models are consistent with
respect to the forcing (atmos, radiation), the parameters, the domain,
and the prognostic land components of the model.
$(DocStringExtensions.FIELDS)
"""
struct SoilSnowModel{
    FT,
    SnM <: Snow.SnowModel{FT},
    SoM <: Soil.EnergyHydrology{FT},
} <: AbstractLandModel{FT}
    "The snow model to be used"
    snow::SnM
    "The soil model to be used"
    soil::SoM
    function SoilSnowModel{FT}(; snow, soil) where {FT <: AbstractFloat}
        prognostic_land_components = (:snow, :soil)
        top_soil_bc = soil.boundary_conditions.top
        snow_bc = snow.boundary_conditions
        @assert top_soil_bc.prognostic_land_components ==
                prognostic_land_components
        @assert snow_bc.prognostic_land_components == prognostic_land_components

        @assert top_soil_bc.atmos == snow_bc.atmos
        @assert top_soil_bc.radiation == snow_bc.radiation

        @assert Domains.obtain_surface_domain(soil.domain) == snow.domain

        @assert soil.parameters.earth_param_set ==
                snow.parameters.earth_param_set
        new{FT, typeof(snow), typeof(soil)}(snow, soil)
    end

end

"""
    lsm_aux_vars(m::SoilSnowModel)

The names of the additional auxiliary variables that are
included in the integrated Soil-Snow model.
"""
lsm_aux_vars(m::SoilSnowModel) = (
    :ground_heat_flux,
    :effective_soil_sfc_T,
    :sfc_scratch,
    :subsfc_scratch,
    :effective_soil_sfc_depth,
)
"""
    lsm_aux_types(m::SoilSnowModel)

The types of the additional auxiliary variables that are
included in the integrated Soil-Snow model.
"""
lsm_aux_types(m::SoilSnowModel{FT}) where {FT} = (FT, FT, FT, FT, FT)

"""
    lsm_aux_domain_names(m::SoilSnowModel)

The domain names of the additional auxiliary variables that are
included in the integrated Soil-Snow model.
"""
lsm_aux_domain_names(m::SoilSnowModel) =
    (:surface, :surface, :surface, :subsurface, :surface)

"""
    make_update_boundary_fluxes(
        land::SoilSnowModel{FT, SnM, SoM},
    ) where {
        FT,
        SnM <: Snow.SnowModel{FT},
        SoM <: Soil.EnergyHydrology{FT},
        }

A method which makes a function; the returned function updates the additional
auxiliary variables for the integrated model, as well as updates the boundary
auxiliary variables for all component models.

This function is called each ode function evaluation, prior to the tendency function
evaluation.

In this method, we
1. Compute the ground heat flux between soil and snow. This is required to update the snow and soil boundary fluxes
2. Update the snow boundary fluxes, which also computes any excess flux of energy or water which occurs when the snow
completely melts in a step. In this case, that excess must go to the soil for conservation
3. Update the soil boundary fluxes use precomputed ground heat flux and excess fluxes from snow.
4. Compute the net flux for the atmosphere, which is useful for assessing conservation.
"""
function make_update_boundary_fluxes(
    land::SoilSnowModel{FT, SnM, SoM},
) where {FT, SnM <: Snow.SnowModel{FT}, SoM <: Soil.EnergyHydrology{FT}}
    update_soil_bf! = make_update_boundary_fluxes(land.soil)
    update_snow_bf! = make_update_boundary_fluxes(land.snow)
    function update_boundary_fluxes!(p, Y, t)
        earth_param_set = land.soil.parameters.earth_param_set
        # First compute the ground heat flux in place:
        update_soil_snow_ground_heat_flux!(
            p,
            Y,
            land.soil.parameters,
            land.snow.parameters,
            land.soil.domain,
            FT,
        )
        #Now update snow boundary conditions, which rely on the ground heat flux
        update_snow_bf!(p, Y, t)
        # Now we can update the soil BC, and use the excess fluxes there in order
        # to conserve energy and water
        update_soil_bf!(p, Y, t)
        return nothing
    end
    return update_boundary_fluxes!
end

"""
    update_soil_snow_ground_heat_flux!(p, Y, soil_params, snow_params, soil_domain, FT)

Computes and updates `p.ground_heat_flux` with the ground heat flux. We approximate this
as
    F_g =  - κ_eff (T_snow - T_soil)/Δz_eff,

where:
    κ_eff = κ_soil * κ_snow / (κ_snow * Δz_soil / 2 + κ_soil * Δz_snow / 2)* (Δz_soil + Δz_snow)/2
    Δz_eff =( Δz_soil + Δz_snow)/2

This is what JULES does to compute the diffusive heat flux between snow and soil, for example,
see equation 24 and 25, with k=N, of Best et al, Geosci. Model Dev., 4, 677–699, 2011

However, this is for a multi-layer snow model, with
T_snow and Δz_snow related to the bottom layer.
We only have a single layer snow model, and using Δz_snow = height of snow leads to a very small flux when the snow is deep.
Due to this, we cap Δz_snow at 10 cm.
"""
function update_soil_snow_ground_heat_flux!(
    p,
    Y,
    soil_params,
    snow_params,
    soil_domain,
    FT,
)
    # Thermal conductivities of soil and snow
    κ_snow = p.snow.κ
    κ_soil = ClimaLand.Domains.top_center_to_surface(p.soil.κ)

    # Depths
    Δz_snow = @. lazy(max(p.snow.z_snow, FT(0.1)))
    Δz_soil = soil_domain.fields.Δz_top

    # Temperatures
    T_snow = p.snow.T
    T_soil = ClimaLand.Domains.top_center_to_surface(p.soil.T)

    # compute the flux
    @. p.ground_heat_flux =
        -κ_soil * κ_snow / (κ_snow * Δz_soil / 2 + κ_soil * Δz_snow / 2) *
        (T_snow - T_soil)
    return nothing
end


### Extensions of existing functions to account for prognostic soil/snow
"""
    snow_boundary_fluxes!(
        bc::AtmosDrivenSnowBC,
        prognostic_land_components::Val{(:snow, :soil)},
        model::SnowModel{FT},
        Y,
        p,
        t,
    ) where {FT}

A method of `snow_boundary_fluxes!` which computes the boundary fluxes for the
snow model accounting for a heat flux between the soil and snow.

The snow surface is assumed to be bare (no vegetation).

Currently this is almost identical to the method for snow alone, except for the
inclusion of the ground heat flux (precomputed by the integrated land model).
However, this will change more if e.g. we allow for transmission of radiation
through the snowpack.
"""
function snow_boundary_fluxes!(
    bc::Snow.AtmosDrivenSnowBC,
    prognostic_land_components::Val{(:snow, :soil)},
    model::SnowModel{FT},
    Y,
    p,
    t,
) where {FT}
    turbulent_fluxes!(p.snow.turbulent_fluxes, bc.atmos, model, Y, p, t)
    net_radiation!(p.snow.R_n, bc.radiation, model, Y, p, t)
    P_snow = p.drivers.P_snow
    P_liq = p.drivers.P_liq

    @. p.snow.total_water_flux =
        P_snow +
        (P_liq + p.snow.turbulent_fluxes.vapor_flux - p.snow.water_runoff) *
        p.snow.snow_cover_fraction

    @. p.snow.liquid_water_flux =
        (
            P_liq + p.snow.turbulent_fluxes.vapor_flux * p.snow.q_l -
            p.snow.water_runoff
        ) * p.snow.snow_cover_fraction

    e_flux_falling_snow =
        Snow.energy_flux_falling_snow(bc.atmos, p, model.parameters)
    e_flux_falling_rain =
        Snow.energy_flux_falling_rain(bc.atmos, p, model.parameters)

    # positive fluxes are TOWARDS atmos
    @. p.snow.total_energy_flux =
        e_flux_falling_snow +
        (
            p.snow.turbulent_fluxes.lhf +
            p.snow.turbulent_fluxes.shf +
            p.snow.R_n - p.snow.energy_runoff - p.ground_heat_flux +
            e_flux_falling_rain
        ) * p.snow.snow_cover_fraction
end

"""
    soil_boundary_fluxes!(
        bc::AtmosDrivenFluxBC{<:PrescribedAtmosphere, <:PrescribedRadiativeFluxes},
         prognostic_land_components::Val{(:snow, :soil)},
        soil::EnergyHydrology,
        Y,
        p,
        t,
    )

A method of `ClimaLand.Soil.soil_boundary_fluxes!` which is used for
integrated land surface models; this computes and returns the net
energy and water flux at the surface of the soil for use as boundary
conditions, taking into account the presence of snow on the surface.
"""
function soil_boundary_fluxes!(
    bc::AtmosDrivenFluxBC,
    prognostic_land_components::Val{(:snow, :soil)},
    soil::EnergyHydrology,
    Y,
    p,
    t,
)
    turbulent_fluxes!(p.soil.turbulent_fluxes, bc.atmos, soil, Y, p, t)
    net_radiation!(p.soil.R_n, bc.radiation, soil, Y, p, t)
    # Liquid influx is a combination of precipitation and snowmelt in general
    liquid_influx =
        Soil.compute_liquid_influx(p, soil, prognostic_land_components)
    # This partitions the influx into runoff and infiltration
    Soil.update_infiltration_water_flux!(
        p,
        bc.runoff,
        liquid_influx,
        Y,
        t,
        soil,
    )
    # This computes the energy of the infiltrating water
    infiltration_energy_flux = Soil.compute_infiltration_energy_flux(
        p,
        bc.runoff,
        bc.atmos,
        prognostic_land_components,
        liquid_influx,
        soil,
        Y,
        t,
    )
    # The actual boundary condition is a mix of liquid water infiltration and
    # evaporation. The infiltration already has accounted for snow cover fraction,
    # because the influx it is computed from has accounted for that.
    # The last term, `excess water flux`, arises when snow melts in a timestep but
    # has a nonzero sublimation which was applied for the entire step.
    @. p.soil.top_bc.water =
        p.soil.infiltration +
        (1 - p.snow.snow_cover_fraction) *
        p.soil.turbulent_fluxes.vapor_flux_liq
    @. p.soil.top_bc.heat =
        (1 - p.snow.snow_cover_fraction) * (
            p.soil.R_n +
            p.soil.turbulent_fluxes.lhf +
            p.soil.turbulent_fluxes.shf
        ) +
        p.snow.snow_cover_fraction * p.ground_heat_flux +
        infiltration_energy_flux

    return nothing
end

"""
   compute_liquid_influx(p,
                         model,
                         prognostic_land_components::Val{(:snow, :soil,)},
    )

Returns the liquid water volume flux at the surface of the soil,
accounting for snowmelt and rainfall.
"""
function Soil.compute_liquid_influx(
    p,
    model,
    prognostic_land_components::Val{(:snow, :soil)},
)
    return @. lazy(
        p.snow.water_runoff * p.snow.snow_cover_fraction +
        (1 - p.snow.snow_cover_fraction) * p.drivers.P_liq,
    )
end

"""
    compute_infiltration_energy_flux(
        p,
        runoff,
        atmos,
        prognostic_land_components::Val{(::snow, :soil,)},
        liquid_influx,
        model::EnergyHydrology,
        Y,
        t,
    )

Computes the energy associated with infiltration of
liquid water into the soil.

For a mix of snowmelt and rainfall, we compute the energy flux fraction due to each like:
Infiltration * (P_liq*(1-σ)/Influx * ρe(T_air)+snowmelt_water_flux *σ/Influx * ρe(T_snow))
i.e.,
Infiltration * (f_rain * ρe(T_air) + f_snow * ρe(T_snow)), where f_rain  = 1- f_snow, since
Influx = P_liq * (1-σ) + snowmelt * σ. Here, σ is the snow cover fraction, and ρe(T) refers
to the volumetric internal energy of liquid water at temperature T.
"""
function Soil.compute_infiltration_energy_flux(
    p,
    runoff,
    atmos,
    prognostic_land_components::Val{(:snow, :soil)},
    liquid_influx,
    model::EnergyHydrology,
    Y,
    t,
)
    FT = eltype(Y.soil.ϑ_l)
    earth_param_set = model.parameters.earth_param_set
    infiltration_fraction = @. lazy(
        Soil.compute_infiltration_fraction(p.soil.infiltration, liquid_influx),
    )
    return @. lazy(
        infiltration_fraction * (
            p.drivers.P_liq *
            (1 - p.snow.snow_cover_fraction) *
            Soil.volumetric_internal_energy_liq(p.drivers.T, earth_param_set) +
            p.snow.energy_runoff * p.snow.snow_cover_fraction
        ),
    )
end

function ClimaLand.Soil.sublimation_source(::Val{(:snow, :soil)}, FT)
    return SoilSublimationwithSnow{FT}()
end

"""
    SoilSublimationwithSnow{FT} <: AbstractSoilSource{FT}

Soil Sublimation source type. Used to defined a method
of `ClimaLand.source!` for soil sublimation with snow present;
treated implicitly
in ϑ_l, ρe_int but explicitly in θ_i.

"""
@kwdef struct SoilSublimationwithSnow{FT} <:
              ClimaLand.Soil.AbstractSoilSource{FT}
    explicit::Bool = false
end

"""
     source!(dY::ClimaCore.Fields.FieldVector,
             src::SoilSublimationwithSnow{FT},
             Y::ClimaCore.Fields.FieldVector,
             p::NamedTuple,
             model
             )

Updates `dY.soil.θ_i` in place with a term due to sublimation; this only affects
the surface layer of soil.
"""
function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::SoilSublimationwithSnow{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model,
) where {FT}
    _ρ_i = FT(LP.ρ_cloud_ice(model.parameters.earth_param_set))
    _ρ_l = FT(LP.ρ_cloud_liq(model.parameters.earth_param_set))
    z = model.domain.fields.z
    Δz_top = model.domain.fields.Δz_top # this returns the center-face distance, not layer thickness
    @. dY.soil.θ_i +=
        -p.soil.turbulent_fluxes.vapor_flux_ice *
        (1 - p.snow.snow_cover_fraction) *
        _ρ_l / _ρ_i * heaviside(z + 2 * Δz_top) / (2 * Δz_top) # only apply to top layer, recall that z is negative
    @. dY.soil.∫F_vol_liq_water_dt +=
        -p.soil.turbulent_fluxes.vapor_flux_ice *
        (1 - p.snow.snow_cover_fraction) # The integral of the source is designed to be this
    return nothing
end

function ClimaLand.get_drivers(model::SoilSnowModel)
    return (
        model.snow.boundary_conditions.atmos,
        model.snow.boundary_conditions.radiation,
    )
end
