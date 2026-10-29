nonzero_norm(norm::FT) where {FT} = norm == 0 ? eps(FT) : norm

"""
    update_root_extraction!(p, Y, t, land)

Updates p.root_extraction and p.root_energy_extraction in place to account
for the flux of water and energy between the soil and the canopy via
root extraction.
"""
NVTX.@annotate function update_root_extraction!(p, Y, t, land)
    update_root_extraction!(
        p,
        Y,
        t,
        land,
        land.canopy.hydraulics,
        land.canopy.soil_moisture_stress,
    )
end

"""
    update_root_extraction!(p, Y, t, land, hydraulics::Canopy.PlantHydraulicsModel, sms)

Computes the root extraction at each level in the soil by computing a Darcy flux between the canopy and soil layers. This is only suitable
when using a PlantHydraulicsModel (not the SteadyStateModel), but it works with any soil moisture stress model.
"""
function update_root_extraction!(
    p,
    Y,
    t,
    land,
    hydraulics::Canopy.PlantHydraulicsModel,
    sms,
)
    z = land.soil.domain.fields.z
    (; conductivity_model) = land.canopy.hydraulics.parameters
    LAI = p.canopy.biomass.area_index.leaf
    RAI = p.canopy.biomass.area_index.root
    # Note that we model the flux between each soil layer and the canopy as:
    # Flux = -K_eff x [(ψ_canopy - ψ_soil)/(z_canopy - z_soil) + 1], where
    # K_eff = K_soil K_canopy /(K_canopy + K_soil)

    # Note that in `PrescribedSoil` mode, we compute the flux using K_soil = K_plant(ψ_soil)
    # and K_canopy = K_plant(ψ_canopy). In `PrognosticSoil` mode here, we compute the flux using
    # K_soil = K_soil(ψ_soil) and K_canopy = K_plant(ψ_canopy).
    effective_height =
        Canopy.hydraulic_height(p, Y, land.canopy.hydraulics, land.canopy)
    @. p.root_extraction =
        Canopy.harmonic_mean(LAI, RAI) *
        Canopy.water_flux(
            z,
            effective_height / 2,
            p.soil.ψ,
            p.canopy.hydraulics.ψ,
            p.soil.K,
            Canopy.hydraulic_conductivity(
                conductivity_model,
                p.canopy.hydraulics.ψ,
            ),
        ) *
        Canopy.root_distribution(z, land.canopy.biomass.rooting_depth)
    @. p.root_energy_extraction =
        p.root_extraction * ClimaLand.Soil.volumetric_internal_energy_liq(
            p.soil.T,
            land.soil.parameters.earth_param_set,
        )
end

"""
    update_root_extraction!(p, Y, t, land, hydraulics::Canopy.SteadyStateModel, sms::Canopy.PiecewiseMoistureStressModel)

Computes the root extraction at each level in the soil by assuming a steady state between the canopy and soil layers. This is only suitable
when using a SteadyStateModel and a PiecewiseMoistureStressModel, because we weight the flux from each layer using the soil moisture
stress factor. This ensures that layers give water to the canopy in proportion to the root density and the soil moisture available.
"""
function update_root_extraction!(
    p,
    Y,
    t,
    land,
    hydraulics::Canopy.SteadyStateModel,
    sms::Canopy.PiecewiseMoistureStressModel,
)
    z = land.soil.domain.fields.z
    (; θ_high, θ_low, c) = land.canopy.soil_moisture_stress
    βm =
        @. lazy(compute_piecewise_moisture_stress(θ_high, θ_low, c, p.soil.θ_l))
    root_dryness_weight = @. lazy(
        Canopy.root_distribution(z, land.canopy.biomass.rooting_depth)*βm,
    )
    norm = p.root_dryness_norm # this was computed in Canopy update aux, = ∫ p(z) β(z) dz = ∫ root dryness weight dz
    FT = eltype(z)
    @. p.root_extraction =
        root_dryness_weight * p.canopy.turbulent_fluxes.vapor_flux /
        nonzero_norm(norm) # Integrates to equal transpiration.
    @. p.root_energy_extraction =
        p.root_extraction * ClimaLand.Soil.volumetric_internal_energy_liq(
            p.soil.T,
            land.soil.parameters.earth_param_set,
        )
end

"""
    update_root_extraction!(p, Y, t, land, hydraulics::Canopy.SteadyStateModel, sms)

Computes the root extraction at each level in the soil by assuming a steady state between the canopy and soil layers. This is suitable
when using a SteadyStateModel and arbitrary moisture stress models that are independent of soil moisture (e.g. the Tuzet model, which depends
on leaf water potential, or the no stress model). Here, we weight the flux from each layer by the root density and the hydraulic conductivity, so that layers at the residual water content give no water to the canopy. This may not be stable
for all situations since we treat root extraction explicitly.
"""
function update_root_extraction!(
    p,
    Y,
    t,
    land,
    hydraulics::Canopy.SteadyStateModel,
    sms,
)
    z = land.soil.domain.fields.z
    FT = eltype(z)
    norm = p.scratch1 # A surface scratch field
    root_water_weight = @. lazy(
        Canopy.root_distribution(z, land.canopy.biomass.rooting_depth) *
        p.soil.K,
    )
    ClimaCore.Operators.column_integral_definite!(norm, root_water_weight)
    @. p.root_extraction =
        root_water_weight * p.canopy.turbulent_fluxes.vapor_flux /
        nonzero_norm(norm) # Integrates to equal transpiration.
    @. p.root_energy_extraction =
        p.root_extraction * ClimaLand.Soil.volumetric_internal_energy_liq(
            p.soil.T,
            land.soil.parameters.earth_param_set,
        )
end

"""
    root_water_flux_per_ground_area!(
        fa::ClimaCore.Fields.Field,
        s::PrognosticGroundConditions,
        model::Canopy.PlantHydraulicsModel,
        canopy,
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t,
    )

An extension of the `root_water_flux_per_ground_area!` function,
 which returns the
net flux of water between the
roots and the soil, per unit ground area,
when both soil and plant
hydraulics are modeled prognostically. This is for use in an LSM.

It is computed by summing the flux of water per ground area between
roots and soil at each soil layer.
"""
function Canopy.root_water_flux_per_ground_area!(
    fa::ClimaCore.Fields.Field,
    s::PrognosticGroundConditions,
    model::Union{Canopy.PlantHydraulicsModel, Canopy.SteadyStateModel},
    canopy,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    ClimaCore.Operators.column_integral_definite!(fa, p.root_extraction)
end

"""
    root_energy_flux_per_ground_area!(
        fa_energy::ClimaCore.Fields.Field,
        s::PrognosticGroundConditions,
        model::Canopy.AbstractCanopyEnergyModel,
        canopy,
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t,
    )


A method computing the energy flux associated with the root-soil
water flux, which returns 0 in cases where we do not need to track
this quantity: in this case, when the canopy energy is tracked,
but we are using a `PrescribedSoil` model (non-prognostic soil model).

Note that this energy flux is not typically included in land surface
models. We account for it when the soil model is prognostic because
the soil model includes the energy in the soil water in its energy
balance; therefore, in order to conserve energy, the canopy model
must account for it as well.
"""
function Canopy.root_energy_flux_per_ground_area!(
    fa_energy::ClimaCore.Fields.Field,
    s::PrognosticGroundConditions,
    model::Canopy.AbstractCanopyEnergyModel,
    canopy,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    ClimaCore.Operators.column_integral_definite!(
        fa_energy,
        p.root_energy_extraction,
    )
end

"""
    RootExtraction{FT} <: Soil.AbstractSoilSource{FT}

Concrete type of Soil.AbstractSoilSource, used for dispatch
in an LSM with both soil and plant hydraulic components.

This is paired with the source term `Canopy.PrognosticSoil`:both
are used at the same time,
ensuring that the water flux into the roots is extracted correctly
from the soil; treated explicitly in all prognostic variables.
"""
@kwdef struct RootExtraction{FT} <: ClimaLand.Soil.AbstractSoilSource{FT}
    explicit::Bool = true
end
"""
    ClimaLand.source!(dY::ClimaCore.Fields.FieldVector,
                     src::RootExtraction,
                     Y::ClimaCore.Fields.FieldVector,
                     p::NamedTuple
                     model::EnergyHydrology)

An extension of the `ClimaLand.source!` function,
 which computes source terms for the
soil model; this method returns the water and energy loss/gain due
to root extraction.
"""
NVTX.@annotate function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::RootExtraction,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model::EnergyHydrology,
)
    @. dY.soil.ϑ_l += -1 * p.root_extraction
    @. dY.soil.ρe_int += -1 * p.root_energy_extraction

    # if flow is negative, towards soil -> soil water increases, add in sign here.
    # Currently, these column integrals are done twice rather than add space to the cache.
    # We can revisit this as needed.
    ClimaCore.Operators.column_integral_definite!(
        p.scratch1,
        p.root_energy_extraction,
    )
    @. dY.soil.∫F_e_dt += p.scratch1

    ClimaCore.Operators.column_integral_definite!(p.scratch1, p.root_extraction)
    @. dY.soil.∫F_vol_liq_water_dt += p.scratch1
end

"""
    update_piecewise_soil_moisture_stress!(ground::PrognosticGroundConditions, p, Y, model, canopy)

Updates the soil moisture stress using the piecewise model for a
prognostic soil model, where θ is vertically resolved and prognostic.
"""
NVTX.@annotate function update_piecewise_soil_moisture_stress!(
    ground::PrognosticGroundConditions,
    p,
    Y,
    model,
    canopy,
)
    θ_l = p.soil.θ_l
    (; θ_high, θ_low, c) = model
    z = ClimaCore.Fields.coordinate_field(axes(θ_l)).z
    # normalized distribution for root density
    root_norm = p.scratch1 # A surface scratch field; included in lsm_aux_vars
    # for soilcanopy and land models.
    root_dist =
        @. lazy(Canopy.root_distribution(z, canopy.biomass.rooting_depth))
    # compute the root zone-averaged βm
    ClimaCore.Operators.column_integral_definite!(root_norm, root_dist)
    # root dryness norm was precomputed in the Canopy update_hydraulics! call
    # It is equal to ∫ β(z) p(z) dz. Therefore, the below is the average β.
    @. p.canopy.soil_moisture_stress.βm = p.root_dryness_norm / root_norm

end

"""
     update_steady_state_leaf_water_potential!(p, hydraulics::SteadyStateModel, sms::Canopy.PiecewiseMoistureStressModel, ground::PrognosticGroundConditions{FT}, canopy) where {FT}

Updates the leaf water potential with the estimated steady state value in the case of a prognostic ground condition,
when the PiecewiseMoistureStress function is chosen. We approximate ψ_leaf as the root distribution and
soil moisture stress weighted soil ψ. When every layer is at or below the wilting point, ψ_leaf is the
wilting potential `ψ_wilt`, the limit of the weighted soil ψ as the last layer dries.
"""
function update_steady_state_leaf_water_potential!(
    p,
    hydraulics::Canopy.SteadyStateModel,
    sms::Canopy.PiecewiseMoistureStressModel,
    ground::PrognosticGroundConditions{FT},
    canopy,
) where {FT}
    θ_l = p.soil.θ_l
    cds = ClimaCore.Fields.coordinate_field(axes(θ_l))
    norm = p.root_dryness_norm # A surface scratch field
    (; θ_high, θ_low, c, ψ_wilt) = sms
    root_dryness_weight = @. lazy(
        Canopy.root_distribution(cds.z, canopy.biomass.rooting_depth) *
        compute_piecewise_moisture_stress(θ_high, θ_low, c, θ_l),
    )
    Operators.column_integral_definite!(norm, root_dryness_weight)
    ψ_leaf = p.canopy.hydraulics.ψ
    Operators.column_integral_definite!(
        ψ_leaf,
        @. lazy(p.soil.ψ * root_dryness_weight)
    )
    @. ψ_leaf = ifelse(norm > 0, ψ_leaf / norm, ψ_wilt)
end

"""
     update_steady_state_leaf_water_potential!(p, hydraulics::SteadyStateModel, sms, ground::PrognosticGroundConditions{FT}, canopy) where {FT}

Updates the leaf water potential with the estimated steady state value in the case of a prognostic ground condition for
arbitrary soil moisture stress parameterizations. Here we weight the potential by the hydraulic conductivity in each layer
"""
function update_steady_state_leaf_water_potential!(
    p,
    hydraulics::Canopy.SteadyStateModel,
    sms,
    ground::PrognosticGroundConditions{FT},
    canopy,
) where {FT}
    K = p.soil.K
    cds = ClimaCore.Fields.coordinate_field(axes(K))
    norm = p.scratch1 # A surface scratch field
    root_weight = @. lazy(
        Canopy.root_distribution(cds.z, canopy.biomass.rooting_depth) * K,
    )
    Operators.column_integral_definite!(norm, root_weight)
    field = @. lazy(p.soil.ψ * root_weight / nonzero_norm(norm))
    Operators.column_integral_definite!(p.canopy.hydraulics.ψ, field)
end

"""
    check_moisture_stress_thresholds(sms, soil_parameters)

Check that the thresholds of a `PiecewiseMoistureStressModel` are those that
`Canopy.piecewise_moisture_stress_thresholds` computes from the soil
parameters; other soil moisture stress models need no check.

This allocates, and should only be called at initialization.
"""
check_moisture_stress_thresholds(sms, soil_parameters) = nothing

function check_moisture_stress_thresholds(
    sms::Canopy.PiecewiseMoistureStressModel,
    soil_parameters,
)
    (; ν, θ_r, hydrology_cm) = soil_parameters
    θ_high, θ_low = Canopy.piecewise_moisture_stress_thresholds(
        ν,
        θ_r,
        hydrology_cm,
        sms.ψ_wilt,
    )
    msg = "The thresholds of the PiecewiseMoistureStressModel differ from those of the soil model; construct it with `soil_params = (; ν, θ_r, hydrology_cm)` taken from the soil parameters."
    check_land_equality(sms.θ_high, θ_high, msg)
    check_land_equality(sms.θ_low, θ_low, msg)
end
