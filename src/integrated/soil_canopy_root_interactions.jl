"""
    update_root_extraction!(p, Y, t, land)

Updates p.root_extraction and p.root_energy_extraction in place to account
for the flux of water and energy between the soil and the canopy via
root extraction.
"""
function update_root_extraction!(p, Y, t, land)
    z = land.soil.domain.fields.z
    (; conductivity_model) = land.canopy.hydraulics.parameters
    area_index = p.canopy.hydraulics.area_index
    above_ground_area_index = p.scratch1
    above_ground_area_index .=
        PlantHydraulics.harmonic_mean.(
            getproperty(
                area_index,
                land.canopy.hydraulics.compartment_labels[1],
            ),
            getproperty(area_index, :root),
        )
    # Note that we model the flux between each soil layer and the canopy as:
    # Flux = -K_eff x [(ψ_canopy - ψ_soil)/(z_canopy - z_soil) + 1], where
    # K_eff = K_soil K_canopy /(K_canopy + K_soil)

    # Note that in `PrescribedSoil` mode, we compute the flux using K_soil = K_plant(ψ_soil)
    # and K_canopy = K_plant(ψ_canopy). In `PrognosticSoil` mode here, we compute the flux using
    # K_soil = K_soil(ψ_soil) and K_canopy = K_plant(ψ_canopy).
    @. p.root_extraction =
        above_ground_area_index *
        PlantHydraulics.water_flux(
            z,
            land.canopy.hydraulics.compartment_midpoints[1],
            p.soil.ψ,
            p.canopy.hydraulics.ψ.:1,
            p.soil.K,
            PlantHydraulics.hydraulic_conductivity(
                conductivity_model,
                p.canopy.hydraulics.ψ.:1,
            ),
        ) *
        Canopy.PlantHydraulics.root_distribution(
            z,
            land.canopy.hydraulics.parameters.rooting_depth,
        )

    @. p.root_energy_extraction =
        p.root_extraction * ClimaLand.Soil.volumetric_internal_energy_liq(
            p.soil.T,
            land.soil.parameters.earth_param_set,
        )
end

"""
    PlantHydraulics.root_water_flux_per_ground_area!(
        fa::ClimaCore.Fields.Field,
        s::Union{PrognosticSoilConditions, PrognosticGroundConditions},
        model::Canopy.PlantHydraulics.PlantHydraulicsModel,
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t,
    )

An extension of the `PlantHydraulics.root_water_flux_per_ground_area!` function,
 which returns the
net flux of water between the
roots and the soil, per unit ground area,
when both soil and plant
hydraulics are modeled prognostically. This is for use in an LSM.

It is computed by summing the flux of water per ground area between
roots and soil at each soil layer.
"""
function PlantHydraulics.root_water_flux_per_ground_area!(
    fa::ClimaCore.Fields.Field,
    s::Union{PrognosticSoilConditions, PrognosticGroundConditions},
    model::Canopy.PlantHydraulics.PlantHydraulicsModel,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    ClimaCore.Operators.column_integral_definite!(fa, p.root_extraction)
end

"""
    root_energy_flux_per_ground_area!(
        fa_energy::ClimaCore.Fields.Field,
        s::Union{PrognosticSoilConditions, PrognosticGroundConditions},
        model::Canopy.AbstractCanopyEnergyModel,
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
    s::Union{PrognosticSoilConditions, PrognosticGroundConditions},
    model::Canopy.AbstractCanopyEnergyModel,
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
function ClimaLand.source!(
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
    # We don't update ∫Sdz here because it cancels out by design with the bottom flux into the canopy (within the land model)
end
