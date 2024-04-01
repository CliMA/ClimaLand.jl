export MedlynConductanceParameters, MedlynConductanceModel

abstract type AbstractStomatalConductanceModel{FT} <:
              AbstractCanopyComponent{FT} end

"""
    MedlynConductanceParameters{FT <: AbstractFloat}

The required parameters for the Medlyn stomatal conductance model.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct MedlynConductanceParameters{FT <: AbstractFloat}
    "Relative diffusivity of water vapor (unitless)"
    Drel::FT
    "Minimum stomatal conductance mol/m^2/s"
    g0::FT
    "Slope parameter, inversely proportional to the square root of marginal water use efficiency (Pa^{1/2})"
    g1::FT
end

Base.eltype(::MedlynConductanceParameters{FT}) where {FT} = FT

struct MedlynConductanceModel{FT, MCP <: MedlynConductanceParameters{FT}} <:
       AbstractStomatalConductanceModel{FT}
    parameters::MCP
end

function MedlynConductanceModel{FT}(
    parameters::MedlynConductanceParameters{FT},
) where {FT <: AbstractFloat}
    return MedlynConductanceModel{eltype(parameters), typeof(parameters)}(
        parameters,
    )
end

ClimaLand.name(model::AbstractStomatalConductanceModel) = :conductance
ClimaLand.auxiliary_vars(model::MedlynConductanceModel) =
    (:medlyn_term, :gs, :transpiration)
ClimaLand.auxiliary_types(model::MedlynConductanceModel{FT}) where {FT} =
    (FT, FT, FT)
ClimaLand.auxiliary_domain_names(::MedlynConductanceModel) =
    (:surface, :surface, :surface)
