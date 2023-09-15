export PrescribedCanopyTempModel
abstract type AbstractCanopyEnergyModel{FT} <: AbstractCanopyComponent{FT} end

ClimaLSM.name(model::AbstractCanopyEnergyModel) = :energy

## PrescribedCanopyTemperature equal to atmos temperature
struct PrescribedCanopyTempModel{FT} <: AbstractCanopyEnergyModel{FT} end

ClimaLSM.auxiliary_vars(model::AbstractCanopyEnergyModel) = (:shf, :lhf)
ClimaLSM.auxiliary_types(model::AbstractCanopyEnergyModel{FT}) where {FT} =
    (FT, FT)
canopy_temperature(model::PrescribedCanopyTempModel, canopy, Y, p, t) =
    canopy.atmos.T(t)
