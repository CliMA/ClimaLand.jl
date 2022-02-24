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


# N pools, and the carbon is distrbuted over layers of soil. But to start is
# 2 bulk pools/single layer to start. eventually vertical flow of carbon
# within a pool.

#CLM: four soil temps observed, 25 layers are simulated, 20 hydrologically connected layers

abstract type AbstractCarbonModel{FT} <: AbstractModel{FT} end

struct BulkTwoPools{FT} <: AbstractCarbonModel{FT}
    turnover_time1::FT
    turnover_time2::FT
    NPP::FT
    respired_fraction::FT
end


ClimaLSM.name(model::BulkTwoPools) = :carbon

ClimaLSM.prognostic_vars(::BulkTwoPools) = (:pool,)

ClimaLSM.Domains.coordinates(model::BulkTwoPools{FT}) where {FT} =
    FT.([0.0, 0.0]);

function ClimaLSM.make_rhs(model::BulkTwoPools{FT}) where {FT}
    function rhs!(dY, Y, p, t) # gets the entire Y
        #=
        M = compute_matrix(model, Y, p, t)

        M[1,1] = -1.0/model.turnover_time1(p.soil.T)
        M[1,2] = 0.0
        M[2,1] = FT(1.0) - model.respired_fraction)/model.turnover_time1
        M[2,2] = -1.0/carbon_turnover_time2

        Input = [model.NPP, 0.0]

        dY.carbon.pool  = Input .+ M * Y.carbon.pool
        =#
        dY.carbon.pool[1] = model.NPP - Y.carbon.pool[1] / model.turnover_time1
        dY.carbon.pool[2] =
            (FT(1.0) - model.respired_fraction) * Y.carbon.pool[1] /
            model.turnover_time1 - Y.carbon.pool[2] / model.turnover_time2
    end
    return rhs!
end
turnover_time1 = 10.0 # years
turnover_time2 = 100.0 # years
respired_fraction = 0.75 # unitless
NPP = 60 # PgC/year
carbon = BulkTwoPools{Float64}(
    turnover_time1,
    turnover_time2,
    NPP,
    respired_fraction,
) # sets all the attributes in carbon model
Y, p, coords = initialize(carbon)

Y.carbon.pool[1] = 600.0
Y.carbon.pool[2] = 1500.0
ode_function! = make_ode_function(carbon);
#dY/dt = ode_function(dY,Y,p,t) # updates dY in place
dY = similar(Y)
ode_function!(dY, Y, p, 0.0)
@test sum(abs.(dY.carbon.pool .- [0.0, 0.0])) < 1e-14


t0 = 0.0;
tf = 1000.0; # years
dt = 1.0; # years
prob = ODEProblem(ode_function!, Y, (t0, tf), p);
sol = solve(prob, RK4(); dt = dt, reltol = 1e-6, abstol = 1e-6);

expected = [600.0, 1500.0]

pool1 = [sol.u[k].carbon.pool[1] for k in 1:1:length(sol.t)] # time series
pool2 = [sol.u[k].carbon.pool[2] for k in 1:1:length(sol.t)]
plot(sol.t, pool1)
# Follow up
# [X] Test steady state solution, look at time evolution plot
# time-varying NPP
# turnover times not constant, but parameterized
# time-varying soil moisture and temp 
# rhs in matrix form
# Hooking up with soil?

#=
function NPP(model::BulkTwoPools{FT},npp_model::Standalone)
    return model.NPP
end

function NPP(model::BulkTwoPools{FT}, npp_model::Coupled)
    p.NPP
end
=#
