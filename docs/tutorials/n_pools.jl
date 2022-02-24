using OrdinaryDiffEq: ODEProblem, solve, RK4
using Plots
using Test
using ClimaCore
using LinearAlgebra
using Plots
using DifferentialEquations
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains

# Import the functions we are extending for our model:
import ClimaLSM: name, make_rhs, prognostic_vars
import ClimaLSM.Domains: coordinates

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
    FT.(zeros(model.number_of_pools));

function ClimaLSM.make_rhs(model::BulkNPools{FT}) where {FT}
    function rhs!(dY,Y,p,t) # gets the entire Y
        NPP = compute_NPP(model.NPP, Y,p,t)
        external_input = model.allocation_fraction*NPP
        dY.carbon.pool .= external_input .+ 
        model.transfer_matrix*(model.environmental_limitations ./ model.turnover_time .* Y.carbon.pool)
    end
    return rhs!
end

function Clima2PoolDefault(::Type{FT}, NPP::AbstractNPPModel{FT}) where {FT}
    # Default values here, some of these could also be passed
    # in, e.g. parameters
    labels = ["fast","slow"]
    turnover_time = [10.0,100.0] # years
    environmental_limitations = [1.0,1.0]# making dependent on T, θ
    allocation_fraction = [1.0,0.0] # constant or not?
    respired_fraction = 0.75 # unitless # model for this?
    transfer_matrix = [-1.0 0.0; (1.0-respired_fraction) -1.0] # unitless
    return BulkNPools{FT}(2,
    labels,
    turnover_time,
    allocation_fraction,
    environmental_limitations,
    transfer_matrix,
    NPP)
end
NPP_model = PrescribedNPPModel{Float64}((t) -> 60+ 5*sin(2*pi*t/10)) # PgC/year
carbon = Clima2PoolDefault(Float64, NPP_model)
Y,p,coords = initialize(carbon)

Y.carbon.pool[1] = 600.0
Y.carbon.pool[2] = 1500.0
ode_function! = make_ode_function(carbon);
#dY/dt = ode_function(dY,Y,p,t) # updates dY in place
dY = similar(Y)
ode_function!(dY,Y,p,0.0)
@test sum(abs.(dY.carbon.pool .- [0.0,0.0])) < 1e-14




### to run a simulation:
#=
t0 = 0.0;
tf = 3000.0; # years
dt = 1.0; # years
prob = ODEProblem(ode_function!, Y, (t0, tf), p);
sol = solve(prob, RK4(); dt = dt, reltol = 1e-6, abstol = 1e-6);

expected = [600.0,1500.0]

pool1 = [sol.u[k].carbon.pool[1] for k in 1:1:length(sol.t)] # time series
pool2 = [sol.u[k].carbon.pool[2] for k in 1:1:length(sol.t)]
plot(sol.t, pool1)
=#
# Follow up
# [X] Test steady state solution, look at time evolution plot
# time-varying NPP
# turnover times not constant, but parameterized
# time-varying soil moisture and temp 
# [X] rhs in matrix form?? Test with greater number of pools
# Hooking up with soil?

### Notes from 3/24
#=
           AbstractCarbonModel
- AbstractBulkModel                                       AbstractResolvedModel

- BulkNPools (standalone or wrapped in LSM if they write hook up)
- ClimaDefault3Pool (hookups already written)
- ClimaDefault10Pool


# Bulk Carbon Models
# - Default with set N, set labels, and pre defined matrices, allocation fraction, expected
# - Option for user to make their own (BulkNPools)

# Vertically Resolved Carbon Models (N pools x m layers)

struct Default3Pool <: AbstractCarbonModel
    # defaults but can change with keyword arguments
    number_of_pools = 3
    turnover_time = [1,2,3]


omdel = Defualt3pool() # 
model = Default3Pool(; turnover_time = [4,5,6])
    
Each parameter guess
    τ = [3.0,4,5]
    model = Default3Pool(;turnover_time = τ, )
    #run model
    # get output
    # compute error metric (likelihood)
    # obtain next parameter guess
end

Y.carbon.pool["slow_carbon"] #NamedTuple 
Y.carbon.pool[1] # indices 1,2,3,...N
("slow_carbon" => 1, "fast_carbon" => 2)
slow_carbon_index = where(labels == "slow_carbon")
=#