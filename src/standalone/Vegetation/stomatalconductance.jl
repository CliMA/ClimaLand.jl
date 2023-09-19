export MedlynConductanceParameters, MedlynConductanceModel

abstract type AbstractStomatalConductanceModel{FT} <:
              AbstractCanopyComponent{FT} end

"""
    MedlynConductanceParameters{FT <: AbstractFloat}

The required parameters for the Medlyn stomatal conductance model.
$(DocStringExtensions.FIELDS)
"""
struct MedlynConductanceParameters{FT <: AbstractFloat}
    "Relative diffusivity of water vapor (unitless)"
    # TODO: move to CLIMAParameters"
    Drel::FT
    "Minimum stomatal conductance mol/m^2/s"
    g0::FT
    "Slope parameter, inversely proportional to the square root of marginal water use efficiency (Pa^{1/2})"
    g1::FT
end

"""
    function MedlynConductanceParameters{FT}(;
        Drel = FT(1.6), # unitless
        g0 =  FT(1e-4), # mol/m^2/s 
        g1 = FT(790) # converted from 5 √kPa to units of √Pa
) where{FT}

A constructor supplying default values for the MedlynConductanceParameters struct.
"""
function MedlynConductanceParameters{FT}(;
    Drel = FT(1.6),
    g0 = FT(1e-4),
    g1 = FT(790),
) where {FT}
    return MedlynConductanceParameters{FT}(Drel, g0, g1)
end


struct MedlynConductanceModel{FT} <: AbstractStomatalConductanceModel{FT}
    parameters::MedlynConductanceParameters{FT}
end

ClimaLSM.name(model::AbstractStomatalConductanceModel) = :conductance
ClimaLSM.auxiliary_vars(model::MedlynConductanceModel) =
    (:medlyn_term, :gs, :transpiration)
ClimaLSM.auxiliary_types(model::MedlynConductanceModel{FT}) where {FT} =
    (FT, FT, FT)
ClimaLSM.auxiliary_domain_names(::MedlynConductanceModel) =
    (:surface, :surface, :surface)
