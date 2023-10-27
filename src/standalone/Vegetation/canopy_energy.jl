export PrescribedCanopyTempModel, AbstractCanopyEnergyModel, canopy_temperature
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

ClimaLSM.auxiliary_vars(model::AbstractCanopyEnergyModel) = (:shf, :lhf)
ClimaLSM.auxiliary_types(model::AbstractCanopyEnergyModel{FT}) where {FT} =
    (FT, FT)
ClimaLSM.auxiliary_domain_names(model::AbstractCanopyEnergyModel) =
    (:surface, :surface)

"""
    canopy_temperature(model::PrescribedCanopyTempModel, canopy, Y, p, t)

Returns the canopy temperature under the `PrescribedCanopyTemp` model, 
where the canopy temperature is assumed to be the same as the atmosphere
temperature.
"""
canopy_temperature(model::PrescribedCanopyTempModel, canopy, Y, p, t) =
    canopy.atmos.T(t)
