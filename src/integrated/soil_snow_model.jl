export LandHydrologyModel

"""
    AtmosDrivenFluxBCwithSnow{
        A <: AbstractAtmosphericDrivers,
        B <: AbstractRadiativeDrivers,
        R <: ClimaLand.Soil.Runoff.AbstractRunoffModel,
    } <: ClimaLand.Soil.AbstractAtmosDrivenFluxBC

A type of Soil.AbstractAtmosDrivenFluxBC for use when snow
changes the boundary conditions of the soil, by modifying
the fraction of the surface exposed to the atmosphere
and by providing snowmelt as additional surface infiltration.
"""
struct AtmosDrivenFluxBCwithSnow{
    A <: AbstractAtmosphericDrivers,
    B <: AbstractRadiativeDrivers,
    R <: ClimaLand.Soil.Runoff.AbstractRunoffModel,
} <: ClimaLand.Soil.AbstractAtmosDrivenFluxBC
    "The atmospheric conditions driving the model"
    atmos::A
    "The radiative fluxes driving the model"
    radiation::B
    "The runoff model. The default is no runoff."
    runoff::R
end

function AtmosDrivenFluxBCwithSnow(atmos, radiation)
    return AtmosDrivenFluxBCwithSnow{
        typeof(atmos),
        typeof(radiation),
        ClimaLand.Soil.Runoff.NoRunoff,
    }(
        atmos,
        radiation,
        ClimaLand.Soil.Runoff.NoRunoff(),
    )
end

"""
    struct LandHydrologyModel{
        FT,
        SnM <: Snow.SnowModel{FT},
        SoM <: Soil.EnergyHydrology{FT},
    } <: AbstractLandModel{FT}
        "The snow model to be used"
        snow::SoM
        "The soil model to be used"
        soil::SoM
    end

A concrete type of land model used for simulating systems with
snow and soil (and eventually rivers).
$(DocStringExtensions.FIELDS)
"""
struct LandHydrologyModel{
    FT,
    SoM <: Snow.SnowModel{FT},
    SnM <: Soil.EnergyHydrology{FT},
} <: AbstractLandModel{FT}
    "The snow model to be used"
    snow::SoM
    "The soil model to be used"
    soil::SnM
end

"""
    LandHydrologyModel{FT}(;
        land_args::NamedTuple = (;),
        snow_model_type::Type{SnM},
        snow_args::NamedTuple = (;),
        soil_model_type::Type{SoM},
        soil_args::NamedTuple = (;),
        ) where {
            FT,
            SnM <: Snow.SnowModel{FT},
            SoM <: Soil.EnergyHydrology{FT},
            }

A constructor for the `LandHydrology`, which takes in the concrete model
type and required arguments for each component, constructs those models,
and constructs the `LandHydrologyModel` from them.

Each component model is constructed with everything it needs to be stepped
forward in time, including boundary conditions, source terms, and interaction
terms.
"""
function LandHydrologyModel{FT}(;
    land_args::NamedTuple = (;),
    snow_model_type::Type{SnM},
    snow_args::NamedTuple = (;),
    soil_model_type::Type{SoM},
    soil_args::NamedTuple = (;),
) where {FT, SnM <: Snow.SnowModel, SoM <: Soil.EnergyHydrology{FT}}
    (; atmos, radiation) = land_args
    if :runoff ∈ propertynames(land_args)
        top_bc = ClimaLand.AtmosDrivenFluxBCwithSnow(
            atmos,
            radiation,
            land_args.runoff,
        )
    else #no runoff model
        top_bc = ClimaLand.AtmosDrivenFluxBCwithSnow(atmos, radiation)
    end
    sources = (Soil.PhaseChange{FT}(),)
    zero_flux = Soil.HeatFluxBC((p, t) -> 0.0)
    boundary_conditions = (;
        top = top_bc,
        bottom = Soil.WaterHeatBC(;
            water = Soil.FreeDrainage(),
            heat = zero_flux,
        ),
    )
    soil = soil_model_type(;
        boundary_conditions = boundary_conditions,
        sources = sources,
        soil_args...,
    )
    snow = snow_model_type(; atmos = atmos, radiation = radiation, snow_args...)


    return LandHydrologyModel{FT, typeof(snow), typeof(soil)}(snow, soil)
end

function ClimaLand.Soil.sublimation_source(
    bc::AtmosDrivenFluxBCwithSnow{AbstractAtmosphericDrivers{FT}},
) where {FT}
    return SoilSublimationwithSnow{FT}()
end

"""
    SoilSublimationwithSnow{FT} <: AbstractSoilSource{FT}

Soil Sublimation source type. Used to defined a method
of `ClimaLand.source!` for soil sublimation with snow present.
"""
struct SoilSublimationwithSnow{FT} <: ClimaLand.Soil.AbstractSoilSource{FT} end

"""
     source!(dY::ClimaCore.Fields.FieldVector,
             src::SoilSublimationwithSnow{FT},
             Y::ClimaCore.Fields.FieldVector,
             p::NamedTuple,
             model
             )

Updates dY.soil.θ_i in place with a term due to sublimation; this only affects
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
        _ρ_l / _ρ_i * heaviside(z + 2 * Δz_top) # only apply to top layer, recall that z is negative
end


"""
    lsm_aux_vars(m::LandHydrologyModel)

The names of the additional auxiliary variables that are
included in the integrated Soil-Snow model.
"""
lsm_aux_vars(m::LandHydrologyModel) = (
    :excess_water_flux,
    :excess_heat_flux,
    :atmos_energy_flux,
    :atmos_water_flux,
)
"""
    lsm_aux_types(m::LandHydrologyModel)

The types of the additional auxiliary variables that are
included in the integrated Soil-Snow model.
"""
lsm_aux_types(m::LandHydrologyModel{FT}) where {FT} = (FT, FT, FT, FT)

"""
    lsm_aux_domain_names(m::LandHydrologyModel)

The domain names of the additional auxiliary variables that are
included in the integrated Soil-Snow model.
"""
lsm_aux_domain_names(m::LandHydrologyModel) =
    (:surface, :surface, :surface, :surface)

"""
    make_update_boundary_fluxes(
        land::LandHydrologyModel{FT, SnM, SoM},
    ) where {
        FT,
        SnM <: Snow.SnowModel{FT},
        SoM <: Soil.EnergyHydrology{FT},
        }

A method which makes a function; the returned function
updates the additional auxiliary variables for the integrated model,
as well as updates the boundary auxiliary variables for all component
models. 

This function is called each ode function evaluation, prior to the tendency function
evaluation.
"""
function make_update_boundary_fluxes(
    land::LandHydrologyModel{FT, SnM, SoM},
) where {FT, SnM <: Snow.SnowModel{FT}, SoM <: Soil.EnergyHydrology{FT}}
    update_soil_bf! = make_update_boundary_fluxes(land.soil)
    update_snow_bf! = make_update_boundary_fluxes(land.snow)
    function update_boundary_fluxes!(p, Y, t)
        update_snow_bf!(p, Y, t)
        # Now we have access to the actual applied and initially computed fluxes for snow
        @. p.excess_water_flux =
            (p.snow.total_water_flux - p.snow.applied_water_flux)
        @. p.excess_heat_flux =
            (p.snow.total_energy_flux - p.snow.applied_energy_flux)
        update_soil_bf!(p, Y, t)
        _LH_f0 = FT(LP.LH_f0(land.soil.parameters.earth_param_set))
        _ρ_liq = FT(LP.ρ_cloud_liq(land.soil.parameters.earth_param_set))
        ρe_falling_snow = -_LH_f0 * _ρ_liq # per unit vol of liquid water
        @. p.atmos_energy_flux =
            (1 - p.snow.snow_cover_fraction) * (
                p.soil.turbulent_fluxes.lhf +
                p.soil.turbulent_fluxes.shf +
                p.soil.R_n
            ) +
            p.snow.snow_cover_fraction * (
                p.snow.turbulent_fluxes.lhf +
                p.snow.turbulent_fluxes.shf +
                p.snow.R_n
            ) +
            p.drivers.P_snow * ρe_falling_snow
        @. p.atmos_water_flux =
            p.drivers.P_snow +
            p.drivers.P_liq +
            (1 - p.snow.snow_cover_fraction) * (
                p.soil.turbulent_fluxes.vapor_flux_liq +
                p.soil.turbulent_fluxes.vapor_flux_ice
            ) +
            p.snow.snow_cover_fraction * p.snow.turbulent_fluxes.vapor_flux

    end
    return update_boundary_fluxes!
end


### Extensions of existing functions to account for prognostic soil/snow
"""
    soil_boundary_fluxes!(
        bc::AtmosDrivenFluxBCwithSnow,
        boundary::ClimaLand.TopBoundary,
        soil::EnergyHydrology{FT},
        Δz,
        Y,
        p,
        t,
    ) where {FT}

A method of `ClimaLand.Soil.soil_boundary_fluxes!` which is used for
integrated land surface models; this computes and returns the net
energy and water flux at the surface of the soil for use as boundary
conditions.
"""
function soil_boundary_fluxes!(
    bc::AtmosDrivenFluxBCwithSnow,
    boundary::ClimaLand.TopBoundary,
    soil::EnergyHydrology{FT},
    Δz,
    Y,
    p,
    t,
) where {FT}
    bc = soil.boundary_conditions.top
    p.soil.turbulent_fluxes .= turbulent_fluxes(bc.atmos, soil, Y, p, t)
    p.soil.R_n .= net_radiation(bc.radiation, soil, Y, p, t)
    Soil.Runoff.update_runoff!(
        p,
        bc.runoff,
        p.drivers.P_liq .+ p.snow.water_runoff .* p.snow.snow_cover_fraction .+
        p.excess_water_flux,
        Y,
        t,
        soil,
    )
    @. p.soil.top_bc.water =
        p.soil.infiltration +
        (1 - p.snow.snow_cover_fraction) *
        p.soil.turbulent_fluxes.vapor_flux_liq
    T_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.T)
    @. p.soil.top_bc.heat =
        (1 - p.snow.snow_cover_fraction) * (
            p.soil.R_n +
            p.soil.turbulent_fluxes.lhf +
            p.soil.turbulent_fluxes.shf
        ) + p.excess_heat_flux
end

function ClimaLand.get_drivers(model::LandHydrologyModel)
    return (model.snow.atmos, model.snow.radiation)
end
