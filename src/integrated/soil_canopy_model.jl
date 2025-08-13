export SoilCanopyModel
"""
    struct SoilCanopyModel{
        FT,
        MM <: Soil.Biogeochemistry.SoilCO2Model{FT},
        SM <: Soil.EnergyHydrology{FT},
        VM <: Canopy.CanopyModel{FT},
    } <: AbstractLandModel{FT}
        "The soil microbe model to be used"
        soilco2::MM
        "The soil model to be used"
        soil::SM
        "The canopy model to be used"
        canopy::VM
    end

A concrete type of land model used for simulating systems with a
canopy and a soil component.
$(DocStringExtensions.FIELDS)
"""
struct SoilCanopyModel{
    FT,
    MM <: Soil.Biogeochemistry.SoilCO2Model{FT},
    SM <: Soil.EnergyHydrology{FT},
    VM <: Canopy.CanopyModel{FT},
} <: AbstractLandModel{FT}
    "The soil microbe model to be used"
    soilco2::MM
    "The soil model to be used"
    soil::SM
    "The canopy model to be used"
    canopy::VM
    function SoilCanopyModel{FT}(
        soilco2::MM,
        soil::SM,
        canopy::VM,
    ) where {
        FT,
        MM <: Soil.Biogeochemistry.SoilCO2Model{FT},
        SM <: Soil.EnergyHydrology{FT},
        VM <: Canopy.CanopyModel{FT},
    }
        prognostic_land_components = (:canopy, :soil, :soilco2)
        top_soil_bc = soil.boundary_conditions.top
        canopy_bc = canopy.boundary_conditions

        # Integrated model checks
        @assert top_soil_bc.prognostic_land_components ==
                prognostic_land_components
        @assert canopy_bc.prognostic_land_components ==
                prognostic_land_components

        @assert top_soil_bc.atmos == soilco2.drivers.atmos
        @assert top_soil_bc.atmos == canopy_bc.atmos
        @assert top_soil_bc.radiation == canopy_bc.radiation

        @assert Domains.obtain_surface_domain(soil.domain) == canopy.domain
        @assert Domains.obtain_surface_domain(soilco2.domain) == canopy.domain

        @assert soil.parameters.earth_param_set ==
                soilco2.parameters.earth_param_set
        @assert soil.parameters.earth_param_set ==
                canopy.parameters.earth_param_set

        # SoilCanopyModel-specific checks
        # Runoff and sublimation are also automatically included in the soil model
        @assert RootExtraction{FT}() in soil.sources
        @assert Soil.PhaseChange{FT}() in soil.sources
        @assert canopy.hydraulics.transpiration isa
                Canopy.PlantHydraulics.DiagnosticTranspiration{FT}
        @assert canopy_bc.ground isa PrognosticSoilConditions{FT}
        @assert soilco2.drivers.met isa Soil.Biogeochemistry.PrognosticMet

        return new{FT, MM, SM, VM}(soilco2, soil, canopy)
    end
end

"""
    SoilCanopyModel{FT}(
        forcing,
        LAI,
        earth_param_set,
        domain::Union{ClimaLand.Domains.Column, ClimaLand.Domains.SphericalShell};
        soil = Soil.EnergyHydrology{FT}(
            domain,
            forcing,
            earth_param_set;
            prognostic_land_components = (:canopy, :soil, :soilco2),
            additional_sources = (ClimaLand.RootExtraction{FT}(),),
        ),
        soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(
            domain,
            Soil.Biogeochemistry.SoilDrivers(
               Soil.Biogeochemistry.PrognosticMet(soil.parameters),
                PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5)),
                forcing.atmos,
            ),
        ),
        canopy = Canopy.CanopyModel{FT}(
            Domains.obtain_surface_domain(domain),
            (;
                atmos = forcing.atmos,
                radiation = forcing.radiation,
                ground = ClimaLand.PrognosticSoilConditions{FT}(),
            ),
            LAI,
            earth_param_set;
            prognostic_land_components = (:canopy, :soil, :soilco2),
        ),
   ) where {FT}

A convenience constructor for setting up the default SoilCanpyModel,
where all the parameterizations and parameter values are set to default values
or passed in via the `earth_param_set`. The boundary conditions of all models
correspond to `forcing` with the atmosphere, as specified by `forcing`, a NamedTuple
of the form (;atmos, radiation), with `atmos` an AbstractAtmosphericDriver and `radiation`
and AbstractRadiativeDriver. The leaf area index `LAI` must be provided (prescribed)
as a TimeVaryingInput, and the domain must be a ClimaLand domain with a vertical extent.
"""
function SoilCanopyModel{FT}(
    forcing,
    LAI,
    earth_param_set,
    domain::Union{ClimaLand.Domains.Column, ClimaLand.Domains.SphericalShell};
    soil = Soil.EnergyHydrology{FT}(
        domain,
        forcing,
        earth_param_set;
        prognostic_land_components = (:canopy, :soil, :soilco2),
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
    ),
    soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(
        domain,
        Soil.Biogeochemistry.SoilDrivers(
            Soil.Biogeochemistry.PrognosticMet(soil.parameters),
            PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5)),
            forcing.atmos,
        ),
    ),
    canopy = Canopy.CanopyModel{FT}(
        Domains.obtain_surface_domain(domain),
        (;
            atmos = forcing.atmos,
            radiation = forcing.radiation,
            ground = ClimaLand.PrognosticSoilConditions{FT}(),
        ),
        LAI,
        earth_param_set;
        prognostic_land_components = (:canopy, :soil, :soilco2),
    ),
) where {FT}
    return SoilCanopyModel{FT}(soilco2, soil, canopy)
end



"""
    SoilCanopyModel{FT}(;
        soilco2_type::Type{MM},
        soilco2_args::NamedTuple = (;),
        land_args::NamedTuple = (;),
        soil_model_type::Type{SM},
        soil_args::NamedTuple = (;),
        canopy_component_types::NamedTuple = (;),
        canopy_component_args::NamedTuple = (;),
        canopy_model_args::NamedTuple = (;),
        ) where {
            FT,
            SM <: Soil.EnergyHydrology{FT},
            MM <: Soil.Biogeochemistry.SoilCO2Model{FT},
            }

A constructor for the `SoilCanopyModel`, which takes in the concrete model
type and required arguments for each component, constructs those models,
and constructs the `SoilCanopyModel` from them.

Each component model is constructed with everything it needs to be stepped
forward in time, including boundary conditions, source terms, and interaction
terms.
"""
function SoilCanopyModel{FT}(;
    soilco2_type::Type{MM},
    soilco2_args::NamedTuple = (;),
    land_args::NamedTuple = (;),
    soil_model_type::Type{SM},
    soil_args::NamedTuple = (;),
    canopy_component_types::NamedTuple = (;),
    canopy_component_args::NamedTuple = (;),
    canopy_model_args::NamedTuple = (;),
) where {
    FT,
    SM <: Soil.EnergyHydrology{FT},
    MM <: Soil.Biogeochemistry.SoilCO2Model{FT},
}

    (; atmos, radiation, soil_organic_carbon) = land_args
    # These should always be set by the constructor.
    sources = (RootExtraction{FT}(), Soil.PhaseChange{FT}())
    prognostic_land_components = (:canopy, :soil, :soilco2)
    if :runoff ∈ propertynames(land_args)
        runoff_model = land_args.runoff
    else
        runoff_model = ClimaLand.Soil.Runoff.NoRunoff()
    end
    top_bc = ClimaLand.AtmosDrivenFluxBC(
        atmos,
        radiation,
        runoff_model,
        prognostic_land_components,
    )
    zero_flux = Soil.HeatFluxBC((p, t) -> 0.0)
    boundary_conditions =
        (; top = top_bc, bottom = Soil.EnergyWaterFreeDrainage())
    soil = soil_model_type(;
        boundary_conditions = boundary_conditions,
        sources = sources,
        soil_args...,
    )

    transpiration = Canopy.PlantHydraulics.DiagnosticTranspiration{FT}()
    ground_conditions = PrognosticSoilConditions{FT}()
    if :energy in propertynames(canopy_component_args)
        energy_model = canopy_component_types.energy(
            canopy_component_args.energy.parameters,
        )
    else
        energy_model = PrescribedCanopyTempModel{FT}()
    end

    canopy = Canopy.CanopyModel{FT}(;
        autotrophic_respiration = canopy_component_types.autotrophic_respiration(
            canopy_component_args.autotrophic_respiration...,
        ),
        radiative_transfer = canopy_component_types.radiative_transfer(
            canopy_component_args.radiative_transfer...,
        ),
        photosynthesis = canopy_component_types.photosynthesis(
            canopy_component_args.photosynthesis...,
        ),
        conductance = canopy_component_types.conductance(
            canopy_component_args.conductance...,
        ),
        hydraulics = canopy_component_types.hydraulics(;
            transpiration = transpiration,
            canopy_component_args.hydraulics...,
        ),
        energy = energy_model,
        boundary_conditions = Canopy.AtmosDrivenCanopyBC(
            atmos,
            radiation,
            ground_conditions,
            prognostic_land_components,
        ),
        canopy_model_args...,
    )


    co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
    soilco2_drivers = Soil.Biogeochemistry.SoilDrivers(
        co2_prognostic_soil,
        soil_organic_carbon,
        atmos,
    )
    # Set the soil CO2 BC to being atmospheric CO2
    soilco2_top_bc = Soil.Biogeochemistry.AtmosCO2StateBC()
    soilco2_bot_bc = Soil.Biogeochemistry.SoilCO2FluxBC((p, t) -> 0.0) # no flux
    soilco2_sources = (Soil.Biogeochemistry.MicrobeProduction{FT}(),)

    soilco2_boundary_conditions =
        (; top = soilco2_top_bc, bottom = soilco2_bot_bc)
    soilco2 = soilco2_type(
        soilco2_args.domain,
        soilco2_drivers;
        boundary_conditions = soilco2_boundary_conditions,
        sources = soilco2_sources,
        parameters = soilco2_args.parameters,
    )

    return SoilCanopyModel{FT}(soilco2, soil, canopy)
end

"""
    lsm_aux_vars(m::SoilCanopyModel)

The names of the additional auxiliary variables that are
included in the integrated Soil-Canopy model.

These include the broadband albedo of the land surface
`α_sfc`, defined as the ratio of SW_u/SW_d,
and `T_sfc`, defined as the temperature a blackbody with emissivity
`ϵ_sfc` would have
in order to emit the same LW_u as the land surface does. This is called the
[effective temperature](https://en.wikipedia.org/wiki/Effective_temperature) in some fields,
and is not the same as the skin temperature (defined e.g. Equation 7.13 of  Bonan, 2019, Climate Change and Terrestrial Ecosystem Modeling.  DOI: 10.1017/9781107339217).
"""
lsm_aux_vars(m::SoilCanopyModel) = (
    :root_extraction,
    :root_energy_extraction,
    :LW_u,
    :SW_u,
    :T_sfc,
    :ϵ_sfc,
    :α_sfc,
    :scratch1,
    :scratch2,
)

"""
    lsm_aux_types(m::SoilCanopyModel)

The types of the additional auxiliary variables that are
included in the integrated Soil-Canopy model.
"""
lsm_aux_types(m::SoilCanopyModel{FT}) where {FT} =
    (FT, FT, FT, FT, FT, FT, FT, FT, FT)

"""
    lsm_aux_domain_names(m::SoilCanopyModel)

The domain names of the additional auxiliary variables that are
included in the integrated Soil-Canopy model.
"""
lsm_aux_domain_names(m::SoilCanopyModel) = (
    :subsurface,
    :subsurface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
)

"""
    make_update_boundary_fluxes(
        land::SoilCanopyModel{FT, MM, SM, RM},
    ) where {
        FT,
        MM <: Soil.Biogeochemistry.SoilCO2Model{FT},
        SM <: Soil.RichardsModel{FT},
        RM <: Canopy.CanopyModel{FT}
        }

A method which makes a function; the returned function
updates the additional auxiliary variables for the integrated model,
as well as updates the boundary auxiliary variables for all component
models.

This function is called each ode function evaluation, prior to the tendency function
evaluation.
"""
function make_update_boundary_fluxes(
    land::SoilCanopyModel{FT, MM, SM, RM},
) where {
    FT,
    MM <: Soil.Biogeochemistry.SoilCO2Model{FT},
    SM <: Soil.EnergyHydrology{FT},
    RM <: Canopy.CanopyModel{FT},
}
    update_soil_bf! = make_update_boundary_fluxes(land.soil)
    update_soilco2_bf! = make_update_boundary_fluxes(land.soilco2)
    update_canopy_bf! = make_update_boundary_fluxes(land.canopy)
    function update_boundary_fluxes!(p, Y, t)
        # update root extraction
        update_root_extraction!(p, Y, t, land)
        # Radiation
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

        update_soil_bf!(p, Y, t)
        update_canopy_bf!(p, Y, t)
        update_soilco2_bf!(p, Y, t)
    end
    return update_boundary_fluxes!
end


"""
    lsm_radiant_energy_fluxes!(p, land::SoilCanopyModel{FT},
                                canopy_radiation::Canopy.AbstractRadiationModel{FT},
                                Y,
                                t,
                                ) where {FT}


A function which computes the net radiation at the ground surface
given the canopy radiation model, as well as the upwelling radiation - and hence
effective albedo, emissivity, and temperature -  and the net canopy radiation.

Returns the correct radiative fluxes for bare ground in the case
where the canopy LAI is zero. Note also that this serves the role of
`canopy_radiant_energy_fluxes!`, which computes the net canopy radiation
when the Canopy is run in standalone mode.
"""
function lsm_radiant_energy_fluxes!(
    p,
    land::SoilCanopyModel{FT},
    canopy_radiation::Canopy.AbstractRadiationModel{FT},
    Y,
    t,
) where {FT}
    canopy = land.canopy
    earth_param_set = canopy.parameters.earth_param_set
    _σ = LP.Stefan(earth_param_set)
    LW_d = p.drivers.LW_d
    SW_d = p.drivers.SW_d

    T_canopy = ClimaLand.Canopy.canopy_temperature(canopy.energy, canopy, Y, p)

    α_soil_PAR = p.soil.PAR_albedo
    α_soil_NIR = p.soil.NIR_albedo
    ϵ_soil = land.soil.parameters.emissivity
    T_soil = ClimaLand.Domains.top_center_to_surface(p.soil.T)

    # in W/m^2
    LW_d_canopy = p.scratch1
    LW_u_soil = p.scratch2
    LW_net_canopy = p.canopy.radiative_transfer.LW_n
    SW_net_canopy = p.canopy.radiative_transfer.SW_n
    R_net_soil = p.soil.R_n
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
    # in total: d - u = CANOPY_ABS + (1-α_soil)*CANOPY_TRANS
    # SW upwelling  = reflected par + reflected nir
    @. SW_u = par_d * f_refl_par + f_refl_nir * nir_d

    # net canopy
    @. SW_net_canopy = f_abs_par * par_d + f_abs_nir * nir_d

    # net soil = (1-α)*trans for par and nir
    @. R_net_soil .=
        f_trans_nir * nir_d * (1 - α_soil_NIR) +
        f_trans_par * par_d * (1 - α_soil_PAR)

    # Working through the math, this satisfies: LW_d - LW_u = LW_c + LW_soil
    ϵ_canopy = p.canopy.radiative_transfer.ϵ # this takes into account LAI/SAI
    @. LW_d_canopy = ((1 - ϵ_canopy) * LW_d + ϵ_canopy * _σ * T_canopy^4) # double checked
    @. LW_u_soil = ϵ_soil * _σ * T_soil^4 + (1 - ϵ_soil) * LW_d_canopy # double checked
    # This is a sign inconsistency. Here Rn is positive if towards soil. X_X
    @. R_net_soil += ϵ_soil * LW_d_canopy - ϵ_soil * _σ * T_soil^4 # double checked
    @. LW_net_canopy =
        ϵ_canopy * LW_d - 2 * ϵ_canopy * _σ * T_canopy^4 + ϵ_canopy * LW_u_soil

    @. LW_u = (1 - ϵ_canopy) * LW_u_soil + ϵ_canopy * _σ * T_canopy^4 # double checked
end


### Extensions of existing functions to account for prognostic soil/canopy
"""
    soil_boundary_fluxes!(
        bc::AtmosDrivenFluxBC{<:PrescribedAtmosphere, <:PrescribedRadiativeFluxes},
        prognostic_land_components::Val{(:canopy, :soil,:soilco2,)},
        soil::EnergyHydrology,
        Y,
        p,
        t,
    )

A method of `ClimaLand.Soil.soil_boundary_fluxes!` which is used for
integrated land surface models; this computes and returns the net
energy and water flux at the surface of the soil for use as boundary
conditions when a canopy and Soil CO2  model is also included, though only
the presence of the canopy modifies the soil BC.
"""
function soil_boundary_fluxes!(
    bc::AtmosDrivenFluxBC,
    prognostic_land_components::Val{(:canopy, :soil, :soilco2)},
    model::EnergyHydrology,
    Y,
    p,
    t,
)
    turbulent_fluxes!(p.soil.turbulent_fluxes, bc.atmos, model, Y, p, t)
    # Liquid influx is a combination of precipitation and snowmelt in general
    liquid_influx =
        Soil.compute_liquid_influx(p, model, prognostic_land_components)
    # This partitions the influx into runoff and infiltration
    Soil.update_infiltration_water_flux!(
        p,
        bc.runoff,
        liquid_influx,
        Y,
        t,
        model,
    )
    # This computes the energy of the infiltrating water
    infiltration_energy_flux = Soil.compute_infiltration_energy_flux(
        p,
        bc.runoff,
        bc.atmos,
        prognostic_land_components,
        liquid_influx,
        model,
        Y,
        t,
    )
    # The actual boundary condition is a mix of liquid water infiltration and
    # evaporation.
    @. p.soil.top_bc.water =
        p.soil.infiltration + p.soil.turbulent_fluxes.vapor_flux_liq
    @. p.soil.top_bc.heat =
        -p.soil.R_n +
        p.soil.turbulent_fluxes.lhf +
        p.soil.turbulent_fluxes.shf +
        infiltration_energy_flux
    return nothing
end

"""
    compute_infiltration_energy_flux(
        p,
        runoff,
        atmos,
        prognostic_land_components::Val{(:canopy, :soil, :soilco2)},
        liquid_influx,
        model::EnergyHydrology,
        Y,
        t,
    )

Computes the energy associated with infiltration of
liquid water into the soil; uses the same method
as the standalone soil model.
"""
function Soil.compute_infiltration_energy_flux(
    p,
    runoff,
    atmos,
    prognostic_land_components::Val{(:canopy, :soil, :soilco2)},
    liquid_influx,
    model::EnergyHydrology,
    Y,
    t,
)
    return Soil.compute_infiltration_energy_flux(
        p,
        runoff,
        atmos,
        Val((:soil,)),
        liquid_influx,
        model,
        Y,
        t,
    )
end

"""
   compute_liquid_influx(p,
                         model,
                         prognostic_land_components::Val{(:canopy, :soil, :soilco2)},
    )

Returns the liquid water volume flux at the surface of the soil; in
 a model without snow as a prognostic variable, the influx is
the liquid precipitation as a volume flux.
"""
function Soil.compute_liquid_influx(
    p,
    model,
    prognostic_land_components::Val{(:canopy, :soil, :soilco2)},
)
    return p.drivers.P_liq
end

function ClimaLand.Soil.sublimation_source(
    prognostic_land_components::Val{(:canopy, :soil, :soilco2)},
    FT,
)
    return ClimaLand.Soil.SoilSublimation{FT}()
end

function ClimaLand.get_drivers(model::SoilCanopyModel)
    return (
        model.canopy.boundary_conditions.atmos,
        model.canopy.boundary_conditions.radiation,
        model.soilco2.drivers.soc,
    )
end
