export PrescribedCanopyTempModel,
    BigLeafEnergyModel,
    AbstractCanopyEnergyModel,
    canopy_temperature,
    root_energy_flux_per_ground_area!,
    BigLeafEnergyParameters

abstract type AbstractCanopyEnergyModel{FT} <: AbstractCanopyComponent{FT} end

ClimaLSM.name(model::AbstractCanopyEnergyModel) = :energy

"""
    PrescribedCanopyTempModel{FT} <: AbstractCanopyEnergyModel{FT}

A model for the energy of the canopy which assumes the canopy temperature
is the same as the atmosphere temperature prescribed in the
`PrescribedAtmos` struct. 

No equation for the energy of the canopy is solved.
"""
struct PrescribedCanopyTempModel{FT} <: AbstractCanopyEnergyModel{FT} end

ClimaLSM.auxiliary_vars(model::PrescribedCanopyTempModel) = (:shf, :lhf)
ClimaLSM.auxiliary_types(model::PrescribedCanopyTempModel{FT}) where {FT} =
    (FT, FT)
ClimaLSM.auxiliary_domain_names(model::PrescribedCanopyTempModel) =
    (:surface, :surface)
"""
    canopy_temperature(model::PrescribedCanopyTempModel, canopy, Y, p, t)

Returns the canopy temperature under the `PrescribedCanopyTemp` model, 
where the canopy temperature is assumed to be the same as the atmosphere
temperature.
"""
canopy_temperature(model::PrescribedCanopyTempModel, canopy, Y, p, t) =
    canopy.atmos.T(t)

## Prognostic Canopy Temperature
"""
    BigLeafEnergyParameters{FT <: AbstractFloat}
"""
struct BigLeafEnergyParameters{FT <: AbstractFloat}
    ρc_canopy::FT
end

function BigLeafEnergyParameters{FT}(; ρc_canopy = FT(2e6)) where {FT}
    return BigLeafEnergyParameters{FT}(ρc_canopy)
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

ClimaLSM.auxiliary_vars(model::AbstractCanopyEnergyModel) =
    (:shf, :lhf, :fa_energy_roots)
ClimaLSM.auxiliary_types(model::AbstractCanopyEnergyModel{FT}) where {FT} =
    (FT, FT, FT)
ClimaLSM.auxiliary_domain_names(model::AbstractCanopyEnergyModel) =
    (:surface, :surface, :surface)
"""
    canopy_temperature(model::BigLeafEnergyModel, canopy, Y, p, t)

Returns the canopy temperature under the `BigLeafEnergyModel` model, 
where the canopy temperature is modeled prognostically.
"""
canopy_temperature(model::BigLeafEnergyModel, canopy, Y, p, t) =
    Y.canopy.energy.T

function make_compute_exp_tendency(model::BigLeafEnergyModel, canopy)
    function compute_exp_tendency!(dY, Y, p, t)
        area_index = p.canopy.hydraulics.area_index
        FT = eltype(t)
        h_c = canopy.hydraulics.compartment_surfaces[end]
        ρc_canopy = model.parameters.ρc_canopy

        net_energy_flux = @. -p.canopy.radiative_transfer.LW_n -
           p.canopy.radiative_transfer.SW_n +
           p.canopy.energy.shf +
           p.canopy.energy.lhf - p.canopy.energy.fa_energy_roots

        # To prevent dividing by zero, change 1/(AI x h_c)" to
        # "1/max(AI x h_c, eps(FT))"
        @. dY.canopy.energy.T =
            -net_energy_flux / (
                max(
                    (
                        getproperty(area_index, :leaf) +
                        getproperty(area_index, :stem)
                    ) * h_c,
                    eps(FT),
                ) * ρc_canopy
            )
    end
    return compute_exp_tendency!
end

"""
    root_energy_flux_per_ground_area!(
        fa_energy::ClimaCore.Fields.Field,
        s::PrescribedSoil{FT},
        model::BigLeafEnergyModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t::FT,
    ) where {FT}


A method which updates the ClimaCore.Fields.Field `fa_energy` in place
with  the energy flux associated with the root-soil
water flux for the `BigLeafEnergyModel` run in standalone mode,
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
    model::BigLeafEnergyModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t::FT,
) where {FT}
    fa_energy .= FT(0)
end
