using ClimaCore
using LinearAlgebra
using StaticArrays

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end

using ClimaLSM
using ClimaLSM.Domains
# Import the functions we are extending for our model:
import ClimaLSM: name, make_rhs, prognostic_vars, prognostic_types

abstract type AbstractCarbonModel{FT} <: AbstractModel{FT} end

struct BulkNPools{FT} <: AbstractCarbonModel{FT}
    domain::AbstractDomain{FT}
    number_of_pools::Int
    labels::Vector{String}
    turnover_time::Vector{FT}
    allocation_fraction::Vector{FT}
    environmental_limitations::Vector{FT}
    transfer_matrix::Matrix{FT}
    NPP::FT
end

ClimaLSM.name(model::BulkNPools) = :carbon

ClimaLSM.prognostic_vars(::BulkNPools) = (:pool,)
ClimaLSM.prognostic_types(model::BulkNPools{FT}) where {FT} =
    (SVector{model.number_of_pools, FT},)

function ClimaLSM.make_rhs(model::BulkNPools{FT}) where {FT}
    function rhs!(dY, Y, p, t) # gets the entire Y
        external_input = model.allocation_fraction * model.NPP
        matrix =
            model.transfer_matrix *
            Diagonal(model.environmental_limitations) *
            Diagonal(1.0 ./ model.turnover_time)
        parent(dY.carbon.pool) .= transpose(
            matrix * transpose(parent(Y.carbon.pool)) .+ external_input,
        )
    end
    return rhs!
end

function Clima2PoolDefault(::Type{FT}, domain) where {FT}
    # Default values here, some of these could also be passed
    # in, e.g. parameters
    labels = ["fast", "slow"]
    turnover_time = [10.0, 100.0] # years
    environmental_limitations = [1.0, 1.0]
    allocation_fraction = [1.0, 0.0]
    respired_fraction = 0.75 # unitless
    transfer_matrix = [-1.0 0.0; (1.0-respired_fraction) -1.0] # unitless
    NPP = 60 # PgC/year
    return BulkNPools{FT}(
        domain,
        2,
        labels,
        turnover_time,
        allocation_fraction,
        environmental_limitations,
        transfer_matrix,
        NPP,
    )
end

carbon =
    Clima2PoolDefault(Float64, Column(Float64; zlim = (-1, 0), nelements = 3))
Y, p, coords = initialize(carbon)
function init_carbon(c)
    return SVector{2, Float64}([600.0, 1500.0])
end

Y.carbon.pool = init_carbon.(coords)
ode_function! = make_ode_function(carbon);
dY = similar(Y)
ode_function!(dY, Y, p, 0.0)
@test sum(sum(abs.(parent(dY.carbon.pool)))) < 1e-15
