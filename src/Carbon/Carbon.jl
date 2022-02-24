module Carbon

using LinearAlgebra
using ClimaLSM
using ClimaLSM.Domains

import ClimaLSM: name, make_rhs, prognostic_vars
import ClimaLSM.Domains: coordinates

export BulkNPools, Clima2PoolDefault, PrescribedNPPModel

abstract type AbstractCarbonModel{FT} <: AbstractModel{FT} end

abstract type AbstractNPPModel{FT <: AbstractFloat} end
struct PrescribedNPPModel{FT} <: AbstractNPPModel{FT}
    NPP::Function
end

function compute_NPP(npp_model::PrescribedNPPModel, _, _, t)
    return npp_model.NPP(t)
end


struct BulkNPools{FT} <: AbstractCarbonModel{FT}
    number_of_pools::Int
    labels::Vector{String}
    turnover_time::Vector{FT}
    allocation_fraction::Vector{FT}
    environmental_limitations::Vector{FT}
    transfer_matrix::Matrix{FT}
    NPP::AbstractNPPModel{FT}
end

ClimaLSM.name(model::BulkNPools) = :carbon

ClimaLSM.prognostic_vars(::BulkNPools) = (:pool,)

ClimaLSM.Domains.coordinates(model::BulkNPools{FT}) where {FT} =
    FT.(zeros(model.number_of_pools))

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
