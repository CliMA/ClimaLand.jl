nonzero_norm(norm::FT) where {FT} = norm == 0 ? eps(FT) : norm

"""
    has_root_dryness_norm(hydraulics, sms)
    has_root_dryness_norm(canopy::Canopy.CanopyModel)

Whether an integrated model with these plant hydraulics and soil moisture stress models
caches `p.root_dryness_norm`, the column integral of `root_water_weight`. The steady state
hydraulics read it to distribute root extraction and to compute the leaf water potential,
and the piecewise moisture stress reads it to average βm over the root zone.
"""
has_root_dryness_norm(hydraulics, sms) = false
has_root_dryness_norm(hydraulics::Canopy.SteadyStateModel, sms) = true
has_root_dryness_norm(hydraulics, sms::Canopy.PiecewiseMoistureStressModel) =
    true
has_root_dryness_norm(
    hydraulics::Canopy.SteadyStateModel,
    sms::Canopy.PiecewiseMoistureStressModel,
) = true
has_root_dryness_norm(canopy::Canopy.CanopyModel) =
    has_root_dryness_norm(canopy.hydraulics, canopy.soil_moisture_stress)

"""
    root_water_weight(p, sms, canopy)

Returns, lazily, the weight w(z) of each soil layer in root water uptake: the root
distribution p(z) times the layer stress factor β(z) for the piecewise moisture stress
model, so that layers give water in proportion to the soil moisture available, and times
the soil hydraulic conductivity K(z) for moisture stress models that do not depend on the
soil moisture of each layer, so that layers at the residual water content give no water.
"""
function root_water_weight(p, sms::Canopy.PiecewiseMoistureStressModel, canopy)
    θ_l = p.soil.θ_l
    z = ClimaCore.Fields.coordinate_field(axes(θ_l)).z
    (; θ_high, θ_low, c) = sms
    return @. lazy(
        Canopy.root_distribution(z, canopy.biomass.rooting_depth) *
        compute_piecewise_moisture_stress(θ_high, θ_low, c, θ_l),
    )
end

function root_water_weight(p, sms, canopy)
    K = p.soil.K
    z = ClimaCore.Fields.coordinate_field(axes(K)).z
    return @. lazy(
        Canopy.root_distribution(z, canopy.biomass.rooting_depth) * K,
    )
end

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
    update_root_extraction!(p, Y, t, land, hydraulics::Canopy.SteadyStateModel, sms)

Computes the root extraction at each level in the soil by assuming a steady state between the canopy and soil layers,
so that root extraction integrates to transpiration. The flux from each layer is proportional to its `root_water_weight`:
the root density times the soil moisture stress factor for the piecewise moisture stress model, or times the hydraulic
conductivity for moisture stress models that do not depend on soil moisture (e.g. the Tuzet model, which depends on leaf
water potential, or the no stress model). This may not be stable for all situations since we treat root extraction explicitly.
"""
function update_root_extraction!(
    p,
    Y,
    t,
    land,
    hydraulics::Canopy.SteadyStateModel,
    sms,
)
    weight = root_water_weight(p, sms, land.canopy)
    # p.root_dryness_norm = ∫ weight dz, cached in update_hydraulics!
    @. p.root_extraction =
        weight * p.canopy.turbulent_fluxes.vapor_flux /
        nonzero_norm(p.root_dryness_norm)
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
    # p.root_dryness_norm = ∫ β(z) p(z) dz, cached in update_hydraulics!
    @. p.canopy.soil_moisture_stress.βm = p.root_dryness_norm / root_norm
end

"""
     update_steady_state_leaf_water_potential!(p, hydraulics::SteadyStateModel, sms::Canopy.PiecewiseMoistureStressModel, ground::PrognosticGroundConditions, canopy)

Updates the leaf water potential with the estimated steady state value in the case of a prognostic ground condition,
when the PiecewiseMoistureStress function is chosen. We approximate ψ_leaf as the root distribution and
soil moisture stress weighted soil ψ. When every layer is at or below the wilting point, ψ_leaf is the
wilting potential `ψ_wilt`, the limit of the weighted soil ψ as the last layer dries.
"""
function update_steady_state_leaf_water_potential!(
    p,
    hydraulics::Canopy.SteadyStateModel,
    sms::Canopy.PiecewiseMoistureStressModel,
    ground::PrognosticGroundConditions,
    canopy,
)
    norm = p.root_dryness_norm
    (; ψ_wilt) = sms
    weight = root_water_weight(p, sms, canopy)
    ψ_leaf = p.canopy.hydraulics.ψ
    Operators.column_integral_definite!(ψ_leaf, @. lazy(p.soil.ψ * weight))
    @. ψ_leaf = ifelse(norm > 0, ψ_leaf / norm, ψ_wilt)
end

"""
     update_steady_state_leaf_water_potential!(p, hydraulics::SteadyStateModel, sms, ground::PrognosticGroundConditions, canopy)

Updates the leaf water potential with the estimated steady state value in the case of a prognostic ground condition for
arbitrary soil moisture stress parameterizations. Here we weight the potential by the hydraulic conductivity in each layer.
"""
function update_steady_state_leaf_water_potential!(
    p,
    hydraulics::Canopy.SteadyStateModel,
    sms,
    ground::PrognosticGroundConditions,
    canopy,
)
    norm = p.root_dryness_norm
    weight = root_water_weight(p, sms, canopy)
    Operators.column_integral_definite!(
        p.canopy.hydraulics.ψ,
        @. lazy(p.soil.ψ * weight / nonzero_norm(norm))
    )
end

"""
     update_root_dryness_norm!(p, hydraulics, sms, ground::PrognosticGroundConditions, canopy)

Updates `p.root_dryness_norm` with ∫ w(z) dz, where w(z) is the `root_water_weight` of each
soil layer, in the models that cache it (see `has_root_dryness_norm`); otherwise does nothing.
"""
function update_root_dryness_norm!(
    p,
    hydraulics,
    sms,
    ground::PrognosticGroundConditions,
    canopy,
)
    has_root_dryness_norm(hydraulics, sms) || return nothing
    Operators.column_integral_definite!(
        p.root_dryness_norm,
        root_water_weight(p, sms, canopy),
    )
    return nothing
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
