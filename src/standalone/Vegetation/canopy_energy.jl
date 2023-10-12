export PrescribedCanopyTempModel,
    BigLeafEnergyModel,
    AbstractCanopyEnergyModel,
    canopy_temperature,
    root_energy_flux_per_ground_area!,
    BigLeafEnergyParameters

abstract type AbstractCanopyEnergyModel{FT} <: AbstractCanopyComponent{FT} end

ClimaLSM.name(model::AbstractCanopyEnergyModel) = :energy


ClimaLSM.auxiliary_vars(model::AbstractCanopyEnergyModel) =
    (:shf, :lhf, :fa_energy_roots, :r_ae)
ClimaLSM.auxiliary_types(model::AbstractCanopyEnergyModel{FT}) where {FT} =
    (FT, FT, FT, FT)
ClimaLSM.auxiliary_domain_names(model::AbstractCanopyEnergyModel) =
    (:surface, :surface, :surface, :surface)

"""
    PrescribedCanopyTempModel{FT} <: AbstractCanopyEnergyModel{FT}

A model for the energy of the canopy which assumes the canopy temperature
is the same as the atmosphere temperature prescribed in the
`PrescribedAtmos` struct.

No equation for the energy of the canopy is solved.
"""
struct PrescribedCanopyTempModel{FT} <: AbstractCanopyEnergyModel{FT} end

"""
    canopy_temperature(model::PrescribedCanopyTempModel, canopy, Y, p, t)

Returns the canopy temperature under the `PrescribedCanopyTemp` model,
where the canopy temperature is assumed to be the same as the atmosphere
temperature.
"""
function canopy_temperature(
    model::PrescribedCanopyTempModel{FT},
    canopy,
    Y,
    p,
    t,
) where {FT}
    FT.(canopy.atmos.T(t))
end

## Prognostic Canopy Temperature
"""
    BigLeafEnergyParameters{FT <: AbstractFloat}

$(DocStringExtensions.FIELDS)

"""
struct BigLeafEnergyParameters{FT <: AbstractFloat}
    "Specific heat per emitting area [J/m^2/K]"
    ac_canopy::FT
end

function BigLeafEnergyParameters{FT}(; ac_canopy = FT(2e3)) where {FT}
    return BigLeafEnergyParameters{FT}(ac_canopy)
end


"""
    BigLeafEnergyModel{FT} <: AbstractCanopyEnergyModel{FT}
"""
struct BigLeafEnergyModel{FT} <: AbstractCanopyEnergyModel{FT}
    parameters::BigLeafEnergyParameters{FT}
end

ClimaLSM.prognostic_vars(model::BigLeafEnergyModel) = (:T,)
ClimaLSM.prognostic_types(model::BigLeafEnergyModel{FT}) where {FT} = (FT,)
ClimaLSM.prognostic_domain_names(model::BigLeafEnergyModel) = (:surface,)

"""
    canopy_temperature(model::BigLeafEnergyModel, canopy, Y, p, t)

Returns the canopy temperature under the `BigLeafEnergyModel` model,
where the canopy temperature is modeled prognostically.
"""
canopy_temperature(model::BigLeafEnergyModel, canopy, Y, p, t) =
    Y.canopy.energy.T

function make_compute_exp_tendency(
    model::BigLeafEnergyModel{FT},
    canopy,
) where {FT}
    function compute_exp_tendency!(dY, Y, p, t)
        area_index = p.canopy.hydraulics.area_index
        ac_canopy = model.parameters.ac_canopy
        # Energy Equation:
        # (ρc_canopy h_canopy AI) ∂T∂t = -∑F
        # or( ac_canopy AI)∂T∂t = -∑F
        # where ∑F = F_sfc - F_bot, and both F_sfc and F_bot are per
        # unit area ground [W/m^2].
        # Because they are per unit area ground, we need the factor of
        # area index on the LHF, as ac_canopy [J/m^2/K]
        # is per unit area plant.

        net_energy_flux = @. -p.canopy.radiative_transfer.LW_n -
           p.canopy.radiative_transfer.SW_n +
           p.canopy.energy.shf +
           p.canopy.energy.lhf - p.canopy.energy.fa_energy_roots

        # To prevent dividing by zero, change AI" to
        # "max(AI, eps(FT))"
        c_per_ground_area =
            @. ac_canopy * max(area_index.leaf + area_index.stem, eps(FT))
        @. dY.canopy.energy.T = -net_energy_flux / c_per_ground_area
    end
    return compute_exp_tendency!
end

"""
    root_energy_flux_per_ground_area!(
        fa_energy::ClimaCore.Fields.Field,
        s::PrescribedSoil{FT},
        model::AbstractCanopyEnergyModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t,
    ) where {FT}


A method which updates the ClimaCore.Fields.Field `fa_energy` in place
with  the energy flux associated with the root-soil
water flux for the `CanopyModel` run in standalone mode,
with a `PrescribedSoil` model.This value is ignored and set to zero
in this case.

Background information: This energy
flux is not typically included in land surface
models. We account for it when the soil model is prognostic because
the soil model includes the energy in the soil water in its energy
balance; therefore, in order to conserve energy, the canopy model
must account for it as well.
"""
function root_energy_flux_per_ground_area!(
    fa_energy::ClimaCore.Fields.Field,
    s::PrescribedSoil{FT},
    model::AbstractCanopyEnergyModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
) where {FT}
    fa_energy .= FT(0)
end
