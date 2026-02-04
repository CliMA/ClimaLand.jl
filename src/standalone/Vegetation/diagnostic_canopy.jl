
"""
     DiagnosticCanopyModel{FT, AR, RM, PM, SM, SMSM, SIFM, B, PS, D} <: ClimaLand.AbstractDiagnosticModel{FT}

$(DocStringExtensions.FIELDS)
"""
struct CanopyModel{FT, AR, RM, PM, SM, SMSM, SIFM, BM, B, PSE, D} <:
       ClimaLand.AbstractDiagnosticModel{FT}
    "Autotrophic respiration model, a canopy component model"
    autotrophic_respiration::AR
    "Radiative transfer model, a canopy component model"
    radiative_transfer::RM
    "Photosynthesis model, a canopy component model"
    photosynthesis::PM
    "Stomatal conductance model, a canopy component model"
    conductance::SM
    "Soil moisture stress parameterization, a canopy component model"
    soil_moisture_stress::SMSM
    "SIF model, a canopy component model"
    sif::SIFM
    "Biomass parameterization, a canopy component model"
    biomass::BM
    "Boundary Conditions"
    boundary_conditions::B
    "Shared parameters between component models"
    earth_param_set::PSE
    "Canopy model domain"
    domain::D
end

function DiagnosticCanopyModel{FT}(;
    autotrophic_respiration::AbstractAutotrophicRespirationModel{FT},
    radiative_transfer::AbstractRadiationModel{FT},
    photosynthesis::AbstractPhotosynthesisModel{FT},
    conductance::AbstractStomatalConductanceModel{FT},
    soil_moisture_stress::AbstractSoilMoistureStressModel{FT},
    sif::AbstractSIFModel{FT},
    biomass::PrescribedBiomassModel{FT},
    boundary_conditions::B,
    earth_param_set::PSE,
    domain::Union{
        ClimaLand.Domains.Point,
        ClimaLand.Domains.Plane,
        ClimaLand.Domains.SphericalSurface,
    },
) where {FT, B, PSE}

    if typeof(photosynthesis) <: PModel{FT}
        @assert typeof(conductance) <: PModelConductance{FT} "When using PModel for photosynthesis, you must also use PModelConductance for stomatal conductance"
    end

    if typeof(conductance) <: PModelConductance{FT}
        @assert typeof(photosynthesis) <: PModel{FT} "When using PModelConductance for stomatal conductance, you must also use PModel for photosynthesis"
    end

    args = (
        autotrophic_respiration,
        radiative_transfer,
        photosynthesis,
        conductance,
        soil_moisture_stress,
        sif,
        biomass,
        boundary_conditions,
        earth_param_set,
        domain,
    )
    return DiagnosticCanopyModel{FT, typeof.(args)...}(args...)
end

"""
    function DiagnosticCanopyModel{FT}(
        domain::Union{
            ClimaLand.Domains.Point,
            ClimaLand.Domains.Plane,
            ClimaLand.Domains.SphericalSurface,
        },
        forcing::NamedTuple,
        LAI::AbstractTimeVaryingInput,
        toml_dict::CP.ParamDict;
        prognostic_land_components = (:canopy,),
        autotrophic_respiration = AutotrophicRespirationModel{FT}(toml_dict),
        radiative_transfer = TwoStreamModel{FT}(domain, toml_dict),
        photosynthesis = FarquharModel{FT}(domain, toml_dict),
        conductance = MedlynConductanceModel{FT}(domain, toml_dict),
        soil_moisture_stress = TuzetMoistureStressModel{FT}(toml_dict),
        biomass= PrescribedBiomassModel{FT}(domain, LAI, toml_dict),
        sif = Lee2015SIFModel{FT}(toml_dict),
        turbulent_flux_parameterization = MoninObukhovCanopyFluxes(toml_dict, biomass.height),
    ) where {FT, PSE}

"""
function DiagnosticCanopyModel{FT}(
    domain::Union{
        ClimaLand.Domains.Point,
        ClimaLand.Domains.Plane,
        ClimaLand.Domains.SphericalSurface,
    },
    forcing::NamedTuple,
    LAI::AbstractTimeVaryingInput,
    toml_dict::CP.ParamDict;
    prognostic_land_components = (:canopy,),
    autotrophic_respiration = AutotrophicRespirationModel{FT}(toml_dict),
    radiative_transfer = TwoStreamModel{FT}(domain, toml_dict),
    photosynthesis = FarquharModel{FT}(domain, toml_dict),
    conductance = MedlynConductanceModel{FT}(domain, toml_dict),
    soil_moisture_stress = TuzetMoistureStressModel{FT}(toml_dict),
    biomass = PrescribedBiomassModel{FT}(domain, LAI, toml_dict),
    turbulent_flux_parameterization = MoninObukhovCanopyFluxes(
        toml_dict,
        biomass.height,
    ),
    sif = Lee2015SIFModel{FT}(toml_dict),
) where {FT}
    (; atmos, radiation, ground) = forcing

    # Confirm that each spatially-varying parameter is on the correct domain
    for component in [
        autotrophic_respiration,
        radiative_transfer,
        photosynthesis,
        conductance,
        soil_moisture_stress,
        biomass,
        sif,
    ]
        # For component models without parameters, skip the check
        !hasproperty(component, :parameters) && continue

        @assert !(component.parameters isa ClimaCore.Fields.Field) ||
                axes(component.parameters) == domain.space.surface
    end

    # Confirm that the LAI passed agrees with the LAI of the biomass model
    @assert biomass.plant_area_index.LAI == LAI
    boundary_conditions = AtmosDrivenCanopyBC(
        atmos,
        radiation,
        ground,
        turbulent_flux_parameterization,
        prognostic_land_components,
    )

    earth_param_set = LP.LandParameters(toml_dict)
    args = (
        autotrophic_respiration,
        radiative_transfer,
        photosynthesis,
        conductance,
        soil_moisture_stress,
        sif,
        biomass,
        boundary_conditions,
        earth_param_set,
        domain,
    )
    return DiagnosticCanopyModel{FT, typeof.(args)...}(args...)
end

ClimaLand.name(::DiagnosticCanopyModel) = :canopy

"""
    canopy_components(::DiagnosticCanopyModel)

Returns the names of the components of the CanopyModel.

These names are used for storing prognostic and auxiliary variables
in a hierarchical manner within the state vectors.

These names must match the field names of the CanopyModel struct.
"""
canopy_components(::DiagnosticCanopyModel) = (
    :conductance,
    :photosynthesis,
    :radiative_transfer,
    :autotrophic_respiration,
    :sif,
    :soil_moisture_stress,
    :biomass,
)



"""
     ClimaLand.make_update_aux(canopy::CanopyModel)

Creates the `update_aux!` function for the `CanopyModel`

Please note that the plant hydraulics model has auxiliary variables
that are updated in its prognostic `compute_exp_tendency!` function.
While confusing, this is better for performance as it saves looping
over the state vector multiple times.

The other sub-components rely heavily on each other,
so the version of the `CanopyModel` with these subcomponents
has a single update_aux! function, given here.
"""
function ClimaLand.make_update_aux(canopy::DiagnosticCanopyModel)
    function update_aux!(p, Y, t)
        # all of these use the current temperature
        update_biomass!(p, Y, t, canopy.biomass, canopy)
        update_radiative_transfer!(p, Y, t, canopy.radiative_transfer, canopy)
        update_soil_moisture_stress!(p, Y, canopy.soil_moisture_stress, canopy)
        update_photosynthesis!(p, Y, canopy.photosynthesis, canopy)
        update_SIF!(p, Y, canopy.sif, canopy)
        update_canopy_conductance!(p, Y, canopy.conductance, canopy)
        update_autotrophic_respiration!(
            p,
            Y,
            canopy.autotrophic_respiration,
            canopy,
        )
    end
    return update_aux!
end

"""
    ClimaLand.component_temperature(model::DiagnosticCanopyModel, Y, p)

a helper function which returns the component temperature for the canopy
model, which is stored in the aux state.
"""
function ClimaLand.component_temperature(model::DiagnosticCanopyModel, Y, p)
    return p.canopy.energy.T # guess using previous vlaue
end

"""
    ClimaLand.component_specific_humidity(model::CanopyModel, Y, p)

a helper function which returns the surface specific humidity for the canopy
model.
"""
function ClimaLand.component_specific_humidity(model::DiagnosticCanopyModel, Y, p)
    earth_param_set = get_earth_param_set(model)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    surface_flux_params = LP.surface_fluxes_parameters(earth_param_set)
    T_sfc = component_temperature(model, Y, p)
    T_air = p.drivers.T
    P_air = p.drivers.P
    q_air = p.drivers.q
    h_sfc = ClimaLand.surface_height(model, Y, p)
    # Below we approximate the surface air density with the
    # atmospheric density to make it independent of T
    # This makes our estimate of the derivatives more exact later on
    atmos = model.boundary_conditions.atmos
    q_sfc = @. lazy(
        Thermodynamics.q_vap_saturation(
            thermo_params,
            T_sfc,
            Thermodynamics.air_density(thermo_params, T_air, P_air, q_air),
            Thermodynamics.Liquid(),
        ),
    )
    return q_sfc
end

"""
    ClimaLand.surface_displacement_height(model::CanopyModel, Y, p)

a helper function which returns the displacement height for the canopy
model.
"""
function ClimaLand.surface_displacement_height(
    model::DiagnosticCanopyModel{FT},
    Y,
    p,
) where {FT}
    sfp = model.boundary_conditions.turbulent_flux_parameterization
    return sfp.displ
end

"""
    ClimaLand.surface_roughness_model(model::CanopyModel, Y, p)

a helper function which returns the surface roughness model for the canopy
model.
"""
function ClimaLand.surface_roughness_model(
    model::DiagnosticCanopyModel{FT},
    Y,
    p,
) where {FT}
    sfp = model.boundary_conditions.turbulent_flux_parameterization
    return @. lazy(
        SurfaceFluxes.ConstantRoughnessParams{FT}(sfp.z_0m, sfp.z_0b),
    )
end

"""
    ClimaLand.get_update_surface_humidity_function(model::CanopyModel, Y, p)

a helper function which computes and returns the function which updates the guess 
for surface specific humidity to the actual value, for the canopy model.
"""
function ClimaLand.get_update_surface_humidity_function(
    model::DiagnosticCanopyModel,
    Y,
    p,
)
    sfp = model.boundary_conditions.turbulent_flux_parameterization
    Cd = sfp.Cd
    LAI = p.canopy.biomass.area_index.leaf
    r_stomata_canopy = p.canopy.conductance.r_stomata_canopy
    function update_q_vap_sfc_at_a_point(
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        T_sfc,
        u_star,
        z_0m,
        z_0b,
        leaf_Cd,
        LAI,
        r_stomata_canopy,
    )
        FT = eltype(param_set)
        g_leaf = leaf_Cd * max(u_star, FT(1)) * LAI # TODO - change clipping Issue 1600
        g_stomata = 1 / r_stomata_canopy
        g_land = g_stomata * g_leaf / (g_leaf + g_stomata)
        g_h = SurfaceFluxes.heat_conductance(
            param_set,
            ζ,
            u_star,
            inputs,
            z_0m,
            z_0b,
            scheme,
        )

        q_vap_int = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
        q_canopy = inputs.q_vap_sfc_guess

        # Solve for q_sfc analytically to satisfy balance of fluxes:
        # Flux_aero = ρ * g_h * (q_sfc - q_atm)
        # Flux_stom = ρ * (q_canopy - q_sfc) / r_land
        # Equating fluxes: g_h * (q_sfc - q_atm) = (q_canopy - q_sfc) / r_land
        # q_sfc * (g_h + 1/r_land) = q_canopy/r_land + g_h * q_atm
        # q_sfc = (q_canopy + g_h * r_land * q_atm) / (1 + g_h * r_land)

        q_new = (g_land / g_h * q_canopy + q_vap_int) / (1 + g_land / g_h)
        return q_new
    end
    # Closure
    update_q_vap_sfc_field(LAI_val, r_val, leaf_Cd) =
        (args...) ->
            update_q_vap_sfc_at_a_point(args..., leaf_Cd, LAI_val, r_val)
    return @. lazy(update_q_vap_sfc_field(LAI, r_stomata_canopy, Cd))
end

"""
    ClimaLand.get_update_surface_temperature_function(model::CanopyModel, Y, p)

a helper function which computes and returns the function which updates the guess 
for surface temperature to the actual value, for the canopy model.
"""
function ClimaLand.get_update_surface_temperature_function(
    model::DiagnosticCanopyModel,
    Y,
    p,
)
    sfp = model.boundary_conditions.turbulent_flux_parameterization
    Cd = sfp.Cd
    AI = @. lazy(
        p.canopy.biomass.area_index.leaf + p.canopy.biomass.area_index.stem,
    )
    function update_T_sfc_at_a_point(
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        u_star,
        z_0m,
        z_0b,
        leaf_Cd,
        AI,
        SW_d,
        LW_d,
        σϵ_ground,
        σϵ_canopy
        
    )
        FT = eltype(param_set)
        Φ_sfc = SurfaceFluxes.surface_geopotential(inputs)
        Φ_int = SurfaceFluxes.interior_geopotential(param_set, inputs)
        T_int = inputs.T_int
        T_canopy = inputs.T_sfc_guess
        g_h = SurfaceFluxes.heat_conductance(
            param_set,
            ζ,
            u_star,
            inputs,
            z_0m,
            z_0b,
            scheme,
        )
        ws = SurfaceFluxes.windspeed(param_set, ζ, u_star, inputs)
        g_land = leaf_Cd * max(u_star, FT(1)) * AI # TODO - change clipping Issue 1600

        ΔΦ = Φ_int - Φ_sfc
        cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
        T_sfc =
            (T_int + T_canopy * g_land / g_h + ΔΦ / cp_d) / (1 + g_land / g_h)
        return T_sfc
    end
    # Closure
    update_T_sfc_field(AI_val, leaf_Cd) =
        (args...) -> update_T_sfc_at_a_point(args..., leaf_Cd, AI_val)
    return @. lazy(update_T_sfc_field(AI, Cd))
end

function make_update_boundary_fluxes(canopy::DiagnosticCanopyModel)
    function update_implicit_boundary_fluxes!(p, Y, t)
        bc = canopy.boundary_conditions
        FT = eltype(earth_param_set)
        par_d = p.canopy.radiative_transfer.par_d
        nir_d = p.canopy.radiative_transfer.nir_d
        f_abs_par = p.canopy.radiative_transfer.par.abs
        f_abs_nir = p.canopy.radiative_transfer.nir.abs
        @. p.canopy.radiative_transfer.SW_n = f_abs_par * par_d + f_abs_nir * nir_d
        ϵ_canopy = p.canopy.radiative_transfer.ϵ # this takes into account LAI/SAI
        # Long wave: use ground conditions from the ground driver
        T_ground = p.drivers.T_ground
        ϵ_ground = ground.ϵ
        _σ = FT(LP.Stefan(earth_param_set))
        LW_d = p.drivers.LW_d
        function LW(T_canopy; T_ground, ϵ_canopy, LW_d, _σ, ϵ_ground)
            LW_d_canopy = (1 - ϵ_canopy) * LW_d + ϵ_canopy * _σ * T_canopy^4
            LW_u_ground = ϵ_ground * _σ * T_ground^4 + (1 - ϵ_ground) * LW_d_canopy
            return ϵ_canopy * LW_d - 2 * ϵ_canopy * _σ * T_canopy^4 + ϵ_canopy * LW_u_ground
        end
        
    # Compute transpiration, SHF, LHF
    ClimaLand.turbulent_fluxes!(
        p.canopy.turbulent_fluxes,
        bc.atmos,
        canopy,
        Y,
        p,
        t,
    )
    # Due to roundoff problem when multiplying and dividing by cp_d, set
    # SHF to zero if LAI < 0.01
    zero_on_lai(X::FT, lai::FT) where {FT} = lai < FT(0.05) ? FT(0) : X
    @. p.canopy.turbulent_fluxes.shf = zero_on_lai(
        p.canopy.turbulent_fluxes.shf,
        p.canopy.biomass.area_index.leaf + p.canopy.biomass.area_index.stem,
    )
    @. p.canopy.turbulent_fluxes.∂shf∂T = zero_on_lai(
        p.canopy.turbulent_fluxes.∂shf∂T,
        p.canopy.biomass.area_index.leaf + p.canopy.biomass.area_index.stem,
    )

    end
    return update_implicit_boundary_fluxes!
end

function ClimaLand.total_energy_per_area!(
    surface_field,
    model::DiagnosticCanopyModel,
    Y,
    p,
    t,
)
    ClimaLand.total_energy_per_area!(surface_field, model.energy, Y, p, t)
end

function ClimaLand.total_liq_water_vol_per_area!(
    surface_field,
    model::DiagnosticCanopyModel,
    Y,
    p,
    t,
)
    ClimaLand.total_liq_water_vol_per_area!(
        surface_field,
        model.hydraulics,
        Y,
        p,
        t,
    )
end
"""
    ClimaLand.make_set_initial_cache(model::CanopyModel)

Set the initial cache `p` for the canopy model. Note that if the photosynthesis model
is the P-model, then `set_initial_cache!` will also run `set_historical_cache!` which
sets the (t-1) values for Vcmax25_opt, Jmax25_opt, and ξ_opt.
"""
function ClimaLand.make_set_initial_cache(model::DiagnosticCanopyModel)
    drivers = get_drivers(model)
    update_drivers! = make_update_drivers(drivers)
    update_cache! = make_update_cache(model)
    function set_initial_cache!(p, Y0, t0)
        update_drivers!(p, t0)
        p.canopy.energy.T .= p.drivers.T # for the first step's update aux
        update_cache!(p, Y0, t0)
        set_historical_cache!(p, Y0, model.photosynthesis, model)
        
    end
    return set_initial_cache!
end

end
