export LandModel
"""
    struct LandModel{
        FT,
        MM <: Soil.Biogeochemistry.SoilCO2Model{FT},
        SM <: Soil.EnergyHydrology{FT},
        VM <: Canopy.CanopyModel{FT},
        SnM <: Snow.SnowModel{FT},
    } <: AbstractLandModel{FT}
        "The soil microbe model to be used"
        soilco2::MM
        "The soil model to be used"
        soil::SM
        "The canopy model to be used"
        canopy::VM
        "The snow model to be used"
        snow::SnM
    end

A concrete type of land model used for simulating systems with
soil, canopy, snow, soilco2.

ClimaLand v1: SoilCO2 is still under testing, but errors in soilco2
do not propagate into the other components.

$(DocStringExtensions.FIELDS)
"""
struct LandModel{
    FT,
    MM <: Soil.Biogeochemistry.SoilCO2Model{FT},
    SM <: Soil.EnergyHydrology{FT},
    VM <: Canopy.CanopyModel{FT},
    SnM <: Snow.SnowModel{FT},
} <: AbstractLandModel{FT}
    "The soil microbe model to be used"
    soilco2::MM
    "The soil model to be used"
    soil::SM
    "The canopy model to be used"
    canopy::VM
    "The snow model to be used"
    snow::SnM
    function LandModel{FT}(
        canopy::VM,
        snow::SnM,
        soil::SM,
        soilco2::MM,
    ) where {
        FT,
        MM <: Soil.Biogeochemistry.SoilCO2Model{FT},
        SM <: Soil.EnergyHydrology{FT},
        VM <: Canopy.CanopyModel{FT},
        SnM <: Snow.SnowModel{FT},
    }
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
        top_soil_bc = soil.boundary_conditions.top
        canopy_bc = canopy.boundary_conditions
        snow_bc = snow.boundary_conditions

        # Integrated model checks
        @assert top_soil_bc.prognostic_land_components ==
                prognostic_land_components
        @assert canopy_bc.prognostic_land_components ==
                prognostic_land_components
        @assert snow_bc.prognostic_land_components == prognostic_land_components

        @assert top_soil_bc.atmos == soilco2.drivers.atmos
        @assert top_soil_bc.atmos == canopy_bc.atmos
        @assert top_soil_bc.radiation == canopy_bc.radiation
        @assert top_soil_bc.atmos == snow_bc.atmos
        @assert top_soil_bc.radiation == snow_bc.radiation

        @assert Domains.obtain_surface_domain(soil.domain) == canopy.domain
        @assert Domains.obtain_surface_domain(soilco2.domain) == canopy.domain
        @assert snow.domain == canopy.domain

        @assert soil.parameters.earth_param_set ==
                soilco2.parameters.earth_param_set
        @assert soil.parameters.earth_param_set == canopy.earth_param_set
        @assert soil.parameters.earth_param_set ==
                snow.parameters.earth_param_set

        # Check that soil moisture stress parameters are consistent between canopy and soil
        if canopy.soil_moisture_stress isa PiecewiseMoistureStressModel
            # Note that these functions allocate. These checks should not occur except on initialization.
            check_land_equality(
                canopy.soil_moisture_stress.θ_high,
                soil.parameters.ν,
            )
            check_land_equality(
                canopy.soil_moisture_stress.θ_low,
                soil.parameters.θ_r,
            )
        end

        # LandModel-specific checks
        # Runoff and sublimation are also automatically included in the soil model
        @assert RootExtraction{FT}() in soil.sources
        @assert Soil.PhaseChange{FT}() in soil.sources
        @assert canopy_bc.ground isa PrognosticGroundConditions{FT}
        @assert soilco2.drivers.met isa PrognosticMet
        comparison = PrognosticMet(soil.parameters)
        # check_land_equality allocates, and should only be used in initialization
        for property in propertynames(soilco2.drivers.met)
            check_land_equality(
                getproperty(soilco2.drivers.met, property),
                getproperty(comparison, property),
            )
        end

        return new{FT, MM, SM, VM, SnM}(soilco2, soil, canopy, snow)
    end
end

"""
    LandModel{FT}(
        forcing,
        LAI,
        toml_dict::CP.ParamDict,
        domain::Union{ClimaLand.Domains.Column, ClimaLand.Domains.SphericalShell},
        Δt;
        soil = Soil.EnergyHydrology{FT}(
            domain,
            forcing,
            toml_dict;
            prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
            additional_sources = (ClimaLand.RootExtraction{FT}(),),
        ),
        soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(
            domain,
            Soil.Biogeochemistry.SoilDrivers(
               PrognosticMet(soil.parameters),
                forcing.atmos,
            ),
            toml_dict,
        ),
        canopy = Canopy.CanopyModel{FT}(
            Domains.obtain_surface_domain(domain),
            (;
                atmos = forcing.atmos,
                radiation = forcing.radiation,
                ground = ClimaLand.PrognosticGroundConditions{FT}(),
            ),
            LAI,
            toml_dict;
            prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        ),
       snow = Snow.SnowModel(
            FT,
            ClimaLand.Domains.obtain_surface_domain(domain),
            forcing,
            toml_dict,
            Δt;
            prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        ),
    ) where {FT}

A convenience constructor for setting up the default LandModel,
where all the parameterizations and parameter values are set to default values
or passed in via the `toml_dict`. The boundary conditions of all models
correspond to `forcing` with the atmosphere, as specified by `forcing`, a NamedTuple
of the form (;atmos, radiation), with `atmos` an AbstractAtmosphericDriver and `radiation`
and AbstractRadiativeDriver. The leaf area index `LAI` must be provided (prescribed)
as a TimeVaryingInput, and the domain must be a ClimaLand domain with a vertical extent.
Finally, since the snow model requires the timestep, that is a required argument as well.
"""
function LandModel{FT}(
    forcing,
    LAI,
    toml_dict::CP.ParamDict,
    domain::Union{
        ClimaLand.Domains.Column,
        ClimaLand.Domains.SphericalShell,
        ClimaLand.Domains.HybridBox,
    },
    Δt;
    soil = Soil.EnergyHydrology{FT}(
        domain,
        forcing,
        toml_dict;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
    ),
    soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(
        domain,
        Soil.Biogeochemistry.SoilDrivers(
            PrognosticMet(soil.parameters),
            forcing.atmos,
        ),
        toml_dict,
    ),
    canopy = Canopy.CanopyModel{FT}(
        Domains.obtain_surface_domain(domain),
        (;
            atmos = forcing.atmos,
            radiation = forcing.radiation,
            ground = ClimaLand.PrognosticGroundConditions{FT}(),
        ),
        LAI,
        toml_dict;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
    ),
    snow = Snow.SnowModel(
        FT,
        ClimaLand.Domains.obtain_surface_domain(domain),
        forcing,
        toml_dict,
        Δt;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
    ),
) where {FT}
    return LandModel{FT}(canopy, snow, soil, soilco2)
end

"""
    ClimaLand.land_components(land::LandModel)

Returns the components of the `LandModel`.

Currently, this method is required in order to preserve an ordering in how we
update the component models' auxiliary states. The canopy `update_aux!` step
depends on snow and soil albedo, but those are only updated in the snow and soil
`update_aux!` steps. So those must occur first (as controlled by the order of
the components returned by `land_components!`.

This needs to be fixed.
"""
ClimaLand.land_components(land::LandModel) = (:soil, :snow, :soilco2, :canopy)
"""
    lsm_aux_vars(m::LandModel)

The names of the additional auxiliary variables that are
included in the land model.
"""
lsm_aux_vars(m::LandModel) = (
    :root_extraction,
    :root_energy_extraction,
    :LW_u,
    :SW_u,
    :scratch1,
    :scratch2,
    :scratch3,
    :excess_water_flux,
    :excess_heat_flux,
    :ground_heat_flux,
    :effective_soil_sfc_T,
    :sfc_scratch,
    :subsfc_scratch,
    :effective_soil_sfc_depth,
    :T_sfc,
    :ϵ_sfc,
    :α_sfc,
    :α_ground,
)

"""
    lsm_aux_types(m::LandModel)

The types of the additional auxiliary variables that are
included in the land model.
"""
lsm_aux_types(m::LandModel{FT}) where {FT} = (
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    NamedTuple{(:PAR, :NIR), Tuple{FT, FT}},
)

"""
    lsm_aux_domain_names(m::LandModel)

The domain names of the additional auxiliary variables that are
included in the land model.
"""
lsm_aux_domain_names(m::LandModel) = (
    :subsurface,
    :subsurface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :subsurface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
)

"""
    make_update_boundary_fluxes(
        land::LandModel{FT, MM, SM, RM, SnM},
    ) where {
        FT,
        MM <: Soil.Biogeochemistry.SoilCO2Model{FT},
        SM <: Soil.RichardsModel{FT},
        RM <: Canopy.CanopyModel{FT}
        SnM <: Snow.SnowModel{FT}
        }

A method which makes a function; the returned function
updates the additional auxiliary variables for the integrated model,
as well as updates the boundary auxiliary variables for all component
models.

This function is called each ode function evaluation, prior to the tendency function
evaluation.
"""
function make_update_boundary_fluxes(
    land::LandModel{FT, MM, SM, RM, SnM},
) where {
    FT,
    MM <: Soil.Biogeochemistry.SoilCO2Model{FT},
    SM <: Soil.EnergyHydrology{FT},
    RM <: Canopy.CanopyModel{FT},
    SnM <: Snow.SnowModel{FT},
}
    update_soil_bf! = make_update_boundary_fluxes(land.soil)
    update_soilco2_bf! = make_update_boundary_fluxes(land.soilco2)
    update_canopy_bf! = make_update_boundary_fluxes(land.canopy)
    update_snow_bf! = make_update_boundary_fluxes(land.snow)

    function update_boundary_fluxes!(p, Y, t)
        earth_param_set = land.soil.parameters.earth_param_set
        # update root extraction
        update_root_extraction!(p, Y, t, land) # defined in src/integrated/soil_canopy_root_interactions.jl

        # Radiation - updates Rn for soil and snow also
        lsm_radiant_energy_fluxes!(
            p,
            land,
            land.canopy.radiative_transfer,
            Y,
            t,
        )

        # Effective (radiative) land properties
        set_eff_land_radiation_properties!(
            p,
            land.soil.parameters.earth_param_set,
        )

        # Compute the ground heat flux in place:
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

        # Now we have access to the actual applied and initially computed fluxes for snow
        @. p.excess_water_flux =
            (p.snow.total_water_flux - p.snow.applied_water_flux)
        @. p.excess_heat_flux =
            (p.snow.total_energy_flux - p.snow.applied_energy_flux)

        # Now we can update the soil BC, and use the precomputed excess
        # fluxes from snow in that function in order to ensure conservation
        update_soil_bf!(p, Y, t)
        # Update canopy
        update_canopy_bf!(p, Y, t)
        # Update soil CO2
        update_soilco2_bf!(p, Y, t)
    end
    return update_boundary_fluxes!
end

"""
    lsm_radiant_energy_fluxes!(p,land::LandModel{FT},
                                canopy_radiation::Canopy.AbstractRadiationModel{FT},
                                Y,
                                t,
                                ) where {FT}


A function which computes the net radiation at the ground surface
give the canopy radiation model, as well as the outgoing radiation,
and the net canopy radiation.

Returns the correct radiative fluxes for bare ground in the case
where the canopy LAI is zero. Note also that this serves the role of
`canopy_radiant_energy_fluxes!`, which computes the net canopy radiation
when the Canopy is run in standalone mode.
"""
function lsm_radiant_energy_fluxes!(
    p,
    land::LandModel{FT},
    canopy_radiation::Canopy.AbstractRadiationModel{FT},
    Y,
    t,
) where {FT}
    canopy = land.canopy
    canopy_bc = canopy.boundary_conditions
    radiation = canopy_bc.radiation
    earth_param_set = canopy.earth_param_set
    _σ = LP.Stefan(earth_param_set)
    LW_d = p.drivers.LW_d
    SW_d = p.drivers.SW_d

    T_canopy = ClimaLand.Canopy.canopy_temperature(canopy.energy, canopy, Y, p)

    α_soil_PAR = p.soil.PAR_albedo
    α_soil_NIR = p.soil.NIR_albedo
    ϵ_soil = land.soil.parameters.emissivity
    T_soil = ClimaLand.Domains.top_center_to_surface(p.soil.T)

    α_snow_NIR = p.snow.α_snow
    α_snow_PAR = p.snow.α_snow
    ϵ_snow = land.snow.parameters.ϵ_snow
    T_snow = p.snow.T_sfc

    # in W/m^2
    LW_d_canopy = p.scratch1
    LW_u_soil = p.scratch2
    LW_u_snow = p.scratch3
    LW_net_canopy = p.canopy.radiative_transfer.LW_n
    SW_net_canopy = p.canopy.radiative_transfer.SW_n
    R_net_soil = p.soil.R_n
    R_net_snow = p.snow.R_n
    LW_u = p.LW_u
    SW_u = p.SW_u
    par_d = p.canopy.radiative_transfer.par_d
    nir_d = p.canopy.radiative_transfer.nir_d
    f_abs_par = p.canopy.radiative_transfer.par.abs
    f_abs_nir = p.canopy.radiative_transfer.nir.abs
    f_refl_par = p.canopy.radiative_transfer.par.refl
    f_refl_nir = p.canopy.radiative_transfer.nir.refl
    f_trans_par = p.canopy.radiative_transfer.par.trans
    f_trans_nir = p.canopy.radiative_transfer.nir.trans
    # in total: d - u = CANOPY_ABS + (1-α_ground)*CANOPY_TRANS
    # SW_u  = reflected par + reflected nir
    @. SW_u = par_d * f_refl_par + f_refl_nir * nir_d

    # net canopy
    @. SW_net_canopy = f_abs_par * par_d + f_abs_nir * nir_d

    # net radiative flux for soil = -((1-α)*trans for par and nir)
    @. R_net_soil .= -(
        f_trans_nir * nir_d * (1 - α_soil_NIR) +
        f_trans_par * par_d * (1 - α_soil_PAR)
    )

    @. R_net_snow .= -(
        f_trans_nir * nir_d * (1 - α_snow_NIR) +
        f_trans_par * par_d * (1 - α_snow_PAR)
    )

    ϵ_canopy = p.canopy.radiative_transfer.ϵ # this takes into account LAI/SAI

    # Working through the math, this satisfies: LW_d - LW_u = LW_c + LW_soil + LW_snow
    @. LW_d_canopy = ((1 - ϵ_canopy) * LW_d + ϵ_canopy * _σ * T_canopy^4) # double checked
    @. LW_u_soil = ϵ_soil * _σ * T_soil^4 + (1 - ϵ_soil) * LW_d_canopy # double checked
    @. LW_u_snow = ϵ_snow * _σ * T_snow^4 + (1 - ϵ_snow) * LW_d_canopy # identical to soil, checked
    @. R_net_soil -= ϵ_soil * LW_d_canopy - ϵ_soil * _σ * T_soil^4 # double checked
    @. R_net_snow -= ϵ_snow * LW_d_canopy - ϵ_snow * _σ * T_snow^4 # identical to soil, checked
    @. LW_net_canopy =
        ϵ_canopy * LW_d - 2 * ϵ_canopy * _σ * T_canopy^4 +
        ϵ_canopy * LW_u_soil * (1 - p.snow.snow_cover_fraction) +
        ϵ_canopy * LW_u_snow * p.snow.snow_cover_fraction # area weighted by snow cover fraction, OK
    @. LW_u =
        (1 - ϵ_canopy) * LW_u_soil * (1 - p.snow.snow_cover_fraction) +
        (1 - ϵ_canopy) * LW_u_snow * p.snow.snow_cover_fraction +
        ϵ_canopy * _σ * T_canopy^4 # area weighed by snow cover fraction, OK
end


### Extensions of existing functions to account for prognostic soil/canopy
"""
    soil_boundary_fluxes!(
        bc::AtmosDrivenFluxBC{<:PrescribedAtmosphere, <:PrescribedRadiativeFluxes},
        prognostic_land_components::Val{(:canopy, :snow, :soil,:soilco2,)},
        soil::EnergyHydrology{FT},
        Y,
        p,
        t,
    ) where {FT}

A method of `ClimaLand.Soil.soil_boundary_fluxes!` which is used for
integrated land surface models; this computes and returns the net
energy and water flux at the surface of the soil for use as boundary
conditions when a canopy and Soil CO2  model is also included, though only
the presence of the canopy modifies the soil BC.
"""
function soil_boundary_fluxes!(
    bc::AtmosDrivenFluxBC,
    prognostic_land_components::Val{(:canopy, :snow, :soil, :soilco2)},
    soil::EnergyHydrology,
    Y,
    p,
    t,
)
    turbulent_fluxes!(p.soil.turbulent_fluxes, bc.atmos, soil, Y, p, t)
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
    @. p.soil.top_bc.water =
        p.soil.infiltration +
        p.excess_water_flux +
        (1 - p.snow.snow_cover_fraction) *
        p.soil.turbulent_fluxes.vapor_flux_liq
    # The actual boundary condition is a mix of liquid water infiltration and
    # evaporation. The infiltration already has accounted for snow cover fraction,
    # because the influx it is computed from has accounted for that.
    # The last term, `excess water flux`, arises when snow melts in a timestep but
    # has a nonzero sublimation which was applied for the entire step.
    @. p.soil.top_bc.heat =
        (1 - p.snow.snow_cover_fraction) * (
            p.soil.R_n +
            p.soil.turbulent_fluxes.lhf +
            p.soil.turbulent_fluxes.shf
        ) +
        p.excess_heat_flux +
        p.snow.snow_cover_fraction * p.ground_heat_flux +
        infiltration_energy_flux

    return nothing
end

"""
   compute_liquid_influx(p,
                         model,
                         prognostic_land_components::Val{(:canopy, :snow, :soil, :soilco2)},
    )

Returns the liquid water volume flux at the surface of the soil; uses
the same method as the soil+snow integrated model.
"""
function Soil.compute_liquid_influx(
    p,
    model,
    prognostic_land_components::Val{(:canopy, :snow, :soil, :soilco2)},
)
    Soil.compute_liquid_influx(p, model, Val((:snow, :soil)))
end

"""
    compute_infiltration_energy_flux(
        p,
        runoff,
        atmos,
        prognostic_land_components:::Val{(:canopy, :snow, :soil, :soilco2)},
        liquid_influx,
        model::EnergyHydrology,
        Y,
        t,
    )

Computes the energy associated with infiltration of
liquid water into the soil; uses the same method as
the soil+snow integrated model.
"""
function Soil.compute_infiltration_energy_flux(
    p,
    runoff,
    atmos,
    prognostic_land_components::Val{(:canopy, :snow, :soil, :soilco2)},
    liquid_influx,
    model::EnergyHydrology,
    Y,
    t,
)
    Soil.compute_infiltration_energy_flux(
        p,
        runoff,
        atmos,
        Val((:snow, :soil)),
        liquid_influx,
        model,
        Y,
        t,
    )
end

function snow_boundary_fluxes!(
    bc::Snow.AtmosDrivenSnowBC,
    prognostic_land_components::Val{(:canopy, :snow, :soil, :soilco2)},
    model::SnowModel{FT},
    Y,
    p,
    t,
) where {FT}
    turbulent_fluxes!(p.snow.turbulent_fluxes, bc.atmos, model, Y, p, t)
    # How does rain affect the below?
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

    # positive fluxes are TOWARDS atmos, but R_n positive if snow absorbs energy
    @. p.snow.total_energy_flux =
        e_flux_falling_snow +
        (
            p.snow.turbulent_fluxes.lhf +
            p.snow.turbulent_fluxes.shf +
            p.snow.R_n - p.snow.energy_runoff - p.ground_heat_flux +
            e_flux_falling_rain
        ) * p.snow.snow_cover_fraction
    return nothing
end


function ClimaLand.Soil.sublimation_source(
    prognostic_land_components::Val{(:canopy, :snow, :soil, :soilco2)},
    FT,
)
    return SoilSublimationwithSnow{FT}()
end

"""
    ClimaLand.get_drivers(model::LandModel)

Returns the "drivers", or forcing variables, for the LandModel.

These consist of atmospheric and radiative forcing.
"""
function ClimaLand.get_drivers(model::LandModel)
    return (
        model.canopy.boundary_conditions.atmos,
        model.canopy.boundary_conditions.radiation,
    )
end

"""
    ClimaLand.surface_albedo(::LandModel, Y, p)

Returns the surface albedo for the land model, which is computed
from the ratio of the upwelling and downwelling shortwave radiation.
"""
ClimaLand.surface_albedo(::LandModel, Y, p) = p.α_sfc

"""
    make_set_initial_cache(model::Union{LandModel, SoilCanopyModel})

Creates the function with arguments (p,Y0,t0) that updates the cache
`p` with initial values corresponding to Y0 and t0.

We require a different method from the default for a model
with a canopy, so this method is for any model with
type ∈ Union{LandModel, SoilCanopyModel}. This is a close copy of
the method for the CanopyModel, except unpacking `model.canopy` rather
than using `model` directly.
"""
function make_set_initial_cache(model::Union{LandModel, SoilCanopyModel})
    update_cache! = make_update_cache(model)
    canopy = model.canopy
    function set_initial_cache!(p, Y0, t0)
        update_cache!(p, Y0, t0)
        Canopy.set_historical_cache!(p, Y0, canopy.photosynthesis, canopy)
        # Make sure that the hydraulics scheme and the biomass scheme are compatible
        hydraulics = canopy.hydraulics
        n_stem = hydraulics.n_stem
        n_leaf = hydraulics.n_leaf
        Canopy.lai_consistency_check.(
            n_stem,
            n_leaf,
            p.canopy.biomass.area_index,
        )
    end
    return set_initial_cache!
end
