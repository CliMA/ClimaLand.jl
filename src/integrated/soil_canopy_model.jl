export SoilCanopyModel
"""
    struct SoilCanopyModel{
        FT,
        SM <: Soil.EnergyHydrology{FT},
        VM <: Canopy.CanopyModel{FT},
    } <: AbstractLandModel{FT}
        soil::SM
        canopy::VM
    end

A concrete type of land model used for simulating systems with a 
canopy and a soil component.
$(DocStringExtensions.FIELDS)
"""
struct SoilCanopyModel{
    FT,
    SM <: Soil.EnergyHydrology{FT},
    VM <: Canopy.CanopyModel{FT},
} <: AbstractLandModel{FT}
    "The soil model to be used"
    soil::SM
    "The canopy model to be used"
    canopy::VM
end

"""
    CanopyRadiativeFluxes{FT} <: AbstractRadiativeDrivers{FT}

A struct used to compute radiative fluxes in land surface models,
indicating that 
canopy absorption and emission is taken into account when computing
radiation at the surface of the soil or snow.

The only other alternative at this stage is
ClimaLSM.PrescribedRadiativeFluxes, where the prescribed downwelling
short and longwave radiative fluxes are used directly,
without accounting for the canopy. There is a different method
of the function `soil_boundary_fluxes` in this case. 
"""
struct CanopyRadiativeFluxes{FT} <: AbstractRadiativeDrivers{FT} end

"""
    soil_boundary_fluxes(
        bc::AtmosDrivenFluxBC{<:PrescribedAtmosphere, <:CanopyRadiativeFluxes},
        boundary::ClimaLSM.TopBoundary,
        model::EnergyHydrology{FT},
        Y,
        Δz,
        p,
        t,
    ) where {FT}

A method of `ClimaLSM.Soil.soil_boundary_fluxes` which is used for
integrated land surface models; this computes and returns the net
energy and water flux at the surface of the soil for use as boundary
conditions.

The net radiative, sensible heat, latent heat, and evaporative fluxes 
are computed and stored in the auxiliary state of the integrated land 
surface models, and this function simply returns those. They are updated
each time step in `update_aux`.
"""
function soil_boundary_fluxes(
    bc::AtmosDrivenFluxBC{<:PrescribedAtmosphere, <:CanopyRadiativeFluxes},
    boundary::ClimaLSM.TopBoundary,
    model::EnergyHydrology{FT},
    Y,
    Δz,
    p,
    t,
) where {FT}
    infiltration = soil_surface_infiltration(
        bc.runoff,
        bc.atmos.liquid_precip(t) .+ p.soil_evap,
        Y,
        p,
        t,
        model.parameters,
    )
    net_energy_flux = @. p.soil_Rn + p.soil_lhf + p.soil_shf
    return infiltration, net_energy_flux
end

"""
    SoilCanopyModel{FT}(;
                         land_args::NamedTuple = (;),
                         soil_model_type::Type{SM},
                         soil_args::NamedTuple = (;),
                         canopy_component_types::NamedTuple = (;),
                         canopy_component_args::NamedTuple = (;),
                         canopy_model_args::NamedTuple = (;),
                         ) where {FT, SM <: Soil.EnergyHydrology{FT}}
A constructor for the `SoilCanopyModel`, which takes in the concrete model
type and required arguments for each component, constructs those models,
and constructs the `SoilCanopyModel` from them.

Each component model is constructed with everything it needs to be stepped
forward in time, including boundary conditions, source terms, and interaction
terms.
"""
function SoilCanopyModel{FT}(;
    land_args::NamedTuple = (;),
    soil_model_type::Type{SM},
    soil_args::NamedTuple = (;),
    canopy_component_types::NamedTuple = (;),
    canopy_component_args::NamedTuple = (;),
    canopy_model_args::NamedTuple = (;),
) where {FT, SM <: Soil.EnergyHydrology{FT}}

    # These may be passed in, or set, depending on use scenario.
    (; atmos, radiation) = land_args
    # These should always be set by the constructor.
    Δz = minimum(
        ClimaCore.Fields.Δz_field(ClimaLSM.coordinates(soil_args.domain)),
    )
    sources = (RootExtraction{FT}(), Soil.PhaseChange{FT}(Δz))
    # add heat BC
    top_bc = ClimaLSM.Soil.AtmosDrivenFluxBC(atmos, CanopyRadiativeFluxes{FT}())
    zero_flux = FluxBC((p, t) -> eltype(t)(0.0))
    boundary_conditions = (;
        top = top_bc,
        bottom = (water = Soil.FreeDrainage(), heat = zero_flux),
    )
    soil = soil_model_type(;
        boundary_conditions = boundary_conditions,
        sources = sources,
        soil_args...,
    )

    transpiration = Canopy.PlantHydraulics.DiagnosticTranspiration{FT}()

    soil_driver = PrognosticSoil(;
        soil_α_PAR = soil.parameters.PAR_albedo,
        soil_α_NIR = soil.parameters.NIR_albedo,
    )

    canopy = Canopy.CanopyModel{FT}(;
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
        soil_driver = soil_driver,
        atmos = atmos,
        radiation = radiation,
        canopy_model_args...,
    )

    return SoilCanopyModel{FT, typeof(soil), typeof(canopy)}(soil, canopy)
end

"""
    interaction_vars(m::SoilCanopyModel)

The names of the additional auxiliary variables that are 
included in the integrated Soil-Canopy model.
"""
interaction_vars(m::SoilCanopyModel) = (
    :root_extraction,
    :soil_Rn,
    :soil_evap,
    :soil_shf,
    :soil_lhf,
    :T_soil,
    :α_soil,
    :ϵ_soil,
)

"""
    interaction_types(m::SoilCanopyModel)

The types of the additional auxiliary variables that are 
included in the integrated Soil-Canopy model.
"""
interaction_types(m::SoilCanopyModel{FT}) where {FT} =
    (FT, FT, FT, FT, FT, FT, FT, FT)

"""
    interaction_domain_names(m::SoilCanopyModel)

The domain names of the additional auxiliary variables that are 
included in the integrated Soil-Canopy model.
"""
interaction_domain_names(m::SoilCanopyModel) = (
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
    make_interactions_update_aux(
        land::SoilCanopyModel{FT, SM, RM},
    ) where {FT, SM <: Soil.RichardsModel{FT}, RM <: Canopy.CanopyModel{FT}}

A method which makes a function; the returned function 
updates the auxiliary variable `p.root_extraction`, which
is needed for a sink term for the soil model and to create the
lower water boundary condition for the canopy model.
It also updates
the soil surface fluxes, which are affected by the presence of a canopy.

This function is called each ode function evaluation.
"""
function make_interactions_update_aux(
    land::SoilCanopyModel{FT, SM, RM},
) where {FT, SM <: Soil.EnergyHydrology{FT}, RM <: Canopy.CanopyModel{FT}}
    function update_aux!(p, Y, t)
        z = ClimaCore.Fields.coordinate_field(land.soil.domain.space).z
        (; conductivity_model) = land.canopy.hydraulics.parameters
        area_index = p.canopy.hydraulics.area_index
        @. p.root_extraction =
            (
                area_index.root + getproperty(
                    area_index,
                    land.canopy.hydraulics.compartment_labels[1],
                )
            ) / 2 *
            PlantHydraulics.flux(
                z,
                land.canopy.hydraulics.compartment_midpoints[1],
                p.soil.ψ,
                p.canopy.hydraulics.ψ[1],
                PlantHydraulics.hydraulic_conductivity(
                    conductivity_model,
                    p.soil.ψ,
                ),
                PlantHydraulics.hydraulic_conductivity(
                    conductivity_model,
                    p.canopy.hydraulics.ψ[1],
                ),
            ) *
            (land.canopy.hydraulics.parameters.root_distribution(z))

        # Soil boundary fluxes under canopy or for bare soil
        bc = land.soil.boundary_conditions.top
        soil_conditions = surface_fluxes(bc.atmos, land.soil, Y, p, t)

        r_soil = ClimaLSM.Domains.top_center_to_surface(
            ClimaLSM.Soil.soil_resistance.(
                p.soil.θ_l,
                Y.soil.ϑ_l,
                Y.soil.θ_i,
                land.soil.parameters,
            ),
        )
        r_ae = soil_conditions.r_ae
        @. p.soil_evap = soil_conditions.vapor_flux * r_ae / (r_soil + r_ae)
        @. p.soil_lhf = soil_conditions.lhf * r_ae / (r_soil + r_ae)
        @. p.soil_shf = soil_conditions.shf
        p.T_soil .= surface_temperature(land.soil, Y, p, t)
        p.α_soil .= surface_albedo(land.soil, Y, p)
        p.ϵ_soil .= surface_emissivity(land.soil, Y, p)
        p.soil_Rn .= net_radiation_at_ground(
            land.canopy.radiative_transfer,
            land.canopy,
            land.soil,
            Y,
            p,
            t,
        )

    end
    return update_aux!
end

"""
    net_radiation_at_ground(
        canopy_radiation::Canopy.AbstractRadiationModel{FT},
        canopy,
        ground_model::Soil.EnergyHydrology,
        Y,
        p,
        t,
    ) where {FT}


A function which computes the net radiation at the ground surface
give the canopy radiation model.

Returns the correct radiative fluxes for bare ground in the case
where the canopy LAI is zero.
"""
function net_radiation_at_ground(
    canopy_radiation::Canopy.AbstractRadiationModel{FT},
    canopy,
    ground_model::Soil.EnergyHydrology,
    Y,
    p,
    t,
) where {FT}
    radiation = canopy.radiation
    earth_param_set = canopy.parameters.earth_param_set
    _σ = LSMP.Stefan(earth_param_set)
    LW_d::FT = radiation.LW_d(t)
    SW_d::FT = radiation.SW_d(t)
    c = FT(LSMP.light_speed(earth_param_set))
    h = FT(LSMP.planck_constant(earth_param_set))
    N_a = FT(LSMP.avogadro_constant(earth_param_set))
    (; λ_γ_PAR, λ_γ_NIR) = canopy_radiation.parameters
    energy_per_photon_PAR = h * c / λ_γ_PAR
    energy_per_photon_NIR = h * c / λ_γ_NIR
    SW_d_beneath_canopy = @. SW_d - (
        (energy_per_photon_PAR * N_a * p.canopy.radiative_transfer.apar) +
        (energy_per_photon_NIR * N_a * p.canopy.radiative_transfer.anir)
    )
    LW_d_beneath_canopy = LW_d # Assumes T_canopy = T_air
    return @. (
        -(1 - p.α_soil) * SW_d_beneath_canopy -
        p.ϵ_soil * (LW_d_beneath_canopy - _σ * p.T_soil^4)
    )
end

"""
    ClimaLSM.source!(dY::ClimaCore.Fields.FieldVector,
                     src::RootExtraction,
                     Y::ClimaCore.Fields.FieldVector,
                     p::NamedTuple
                     model::EnergyHydrology)

An extension of the `ClimaLSM.source!` function,
 which computes source terms for the 
soil model; this method returns the water and energy loss/gain due
to root extraction.
"""
function ClimaLSM.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::RootExtraction,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model::EnergyHydrology,
)
    @. dY.soil.ϑ_l += -1 * p.root_extraction
    @. dY.soil.ρe_int +=
        -1 *
        p.root_extraction *
        volumetric_internal_energy_liq(p.soil.T, model.parameters)
    # if flow is negative, towards soil -> soil water increases, add in sign here.
end
