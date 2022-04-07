module Carbon
using DocStringExtensions
using LinearAlgebra

using ClimaLSM
using ClimaLSM.Domains
import ClimaLSM: name, make_rhs, prognostic_vars
import ClimaLSM.Domains: coordinates

export BulkNPools, Clima2PoolDefault, PrescribedNPPModel

"""
    AbstractCarbonModel{FT} <: AbstractModel{FT}

An abstract type for carbon models.
"""
abstract type AbstractCarbonModel{FT} <: AbstractModel{FT} end

"""
    AbstractNPPModel{FT <: AbstractFloat}

An abstract type of NPP models.
"""
abstract type AbstractNPPModel{FT <: AbstractFloat} end

"""
    PrescribedNPPModel{FT} <: AbstractNPPModel{FT}

The concrete type of NPPModel for use when running the carbon model
in standalone mode, with a prescribed NPP function.
"""
struct PrescribedNPPModel{FT} <: AbstractNPPModel{FT}
    NPP::Function
end

"""
    compute_NPP(npp_model::PrescribedNPPModel, _, _, t)

A function which computes and returns the NPP.
"""
function compute_NPP(npp_model::PrescribedNPPModel, _, _, t)
    return npp_model.NPP(t)
end

"""
    BulkNPools{FT} <: AbstractCarbonModel{FT}

A bulk carbon model with N pools.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct BulkNPools{FT} <: AbstractCarbonModel{FT}
    "The number of carbon poolsN"
    number_of_pools::Int
    "The names of the pools"
    labels::Vector{String}
    "The turnover time for the pools (years)"
    turnover_time::Vector{FT}
    "The allocation fraction of NPP into each pool (unitless)"
    allocation_fraction::Vector{FT}
    "Environmental limitation vector (unitless)"
    environmental_limitations::Vector{FT}
    "The transfer matrix between pools, of size NxN (unitless)"
    transfer_matrix::Matrix{FT}
    "The NPP model"
    NPP::AbstractNPPModel{FT}
end

"""
    ClimaLSM.name(model::BulkNPools)

Returns the name of the bulk carbon model.
"""
ClimaLSM.name(model::BulkNPools) = :carbon

"""
    ClimaLSM.prognostic_vars(model::BulkNPools)

Returns the prognostic variable names of the bulk carbon model.
"""
ClimaLSM.prognostic_vars(::BulkNPools) = (:pool,)

"""
    ClimaLSM.Domains.coordinates(model::BulkNPools{FT})

A stand-in for initializing the carbon prognostic state, which is
an array at each coordinate point.

Currently, we only support one array/one coordinate point.
"""
ClimaLSM.Domains.coordinates(model::BulkNPools{FT}) where {FT} =
    FT.(zeros(model.number_of_pools))

"""
    ClimaLSM.make_rhs(model::BulkNPools{FT}) where {FT}

Returns the rhs! function for the bulk carbon model.
"""
function ClimaLSM.make_rhs(model::BulkNPools{FT}) where {FT}
    function rhs!(dY, Y, p, t)
        NPP = compute_NPP(model.NPP, Y, p, t)
        external_input = model.allocation_fraction * NPP
        dY.carbon.pool .=
            external_input .+
            model.transfer_matrix * (
                model.environmental_limitations ./ model.turnover_time .*
                Y.carbon.pool
            )
    end
    return rhs!
end

"""
    Clima2PoolDefault(::Type{FT}, NPP::AbstractNPPModel{FT}) where {FT}

Returns the default Clima 2 pool bulk carbon model.
"""
function Clima2PoolDefault(::Type{FT}, NPP::AbstractNPPModel{FT}) where {FT}
    labels = ["fast", "slow"]
    turnover_time = [10.0, 100.0] # years
    environmental_limitations = [1.0, 1.0]# making dependent on T, θ
    allocation_fraction = [1.0, 0.0] # constant or not
    respired_fraction = 0.75 # unitless, model for this?
    transfer_matrix = [-1.0 0.0; (1.0-respired_fraction) -1.0] # unitless
    return BulkNPools{FT}(
        2,
        labels,
        turnover_time,
        allocation_fraction,
        environmental_limitations,
        transfer_matrix,
        NPP,
    )
end

end # module
