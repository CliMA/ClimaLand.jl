import ClimaTimeSteppers
import LinearAlgebra

mutable struct ClimaLandSimulation{
    STEPPER,
    U,
    P,
    T,
    CALL,
    }
    stepper::STEPPER
    u::U
    p::P
    t::T
    dt::T
    tf::T
    callbacks::CALL
end

function Base.show(io::IO, sim::ClimaLandSimulation)
    return print(
        io,
        "ClimaLand Simulation\n",
        "└── Current time: $(sim.t)\n",
    )
end

struct TimeStepper{
    ALGO,
    FUNC,
    CACHE,
    }
    algo::ALGO
    func::FUNC
    cache::CACHE
end

# TODO: Remove type restriction from init_cache in CTS so that we don't need to
# copy over this function
function init_cache(prob, alg::ClimaTimeSteppers.IMEXAlgorithm{ClimaTimeSteppers.Unconstrained}; kwargs...)
    (; u0, f) = prob
    (; T_imp!) = f
    (; tableau, newtons_method) = alg
    (; a_exp, b_exp, a_imp, b_imp) = tableau
    s = length(b_exp)
    inds = ntuple(i -> i, s)
    inds_T_exp = filter(i -> !all(iszero, a_exp[:, i]) || !iszero(b_exp[i]), inds)
    inds_T_imp = filter(i -> !all(iszero, a_imp[:, i]) || !iszero(b_imp[i]), inds)
    U = zero(u0)
    T_lim = ClimaTimeSteppers.SparseContainer(map(i -> zero(u0), collect(1:length(inds_T_exp))), inds_T_exp)
    T_exp = ClimaTimeSteppers.SparseContainer(map(i -> zero(u0), collect(1:length(inds_T_exp))), inds_T_exp)
    T_imp = ClimaTimeSteppers.SparseContainer(map(i -> zero(u0), collect(1:length(inds_T_imp))), inds_T_imp)
    temp = zero(u0)
    γs = unique(filter(!iszero, LinearAlgebra.diag(a_imp)))
    γ = length(γs) == 1 ? γs[1] : nothing # TODO: This could just be a constant.
    jac_prototype = ClimaTimeSteppers.has_jac(T_imp!) ? T_imp!.jac_prototype : nothing
    newtons_method_cache =
        isnothing(T_imp!) || isnothing(newtons_method) ? nothing : ClimaTimeSteppers.allocate_cache(newtons_method, u0, jac_prototype)
    return ClimaTimeSteppers.IMEXARKCache(U, T_lim, T_exp, T_imp, temp, γ, newtons_method_cache)
end

# u0 is used as a prototype
function TimeStepper(u0, t0, alg, func)
    prob = (; u0, f = func)
    cache = init_cache(prob, alg)

    isnothing(func.cache!) || func.cache!(u0, cache, t0)
    return TimeStepper(alg, func, cache)
end

function step!(simulation::ClimaLandSimulation)
    simulation.t += simulation.dt
    integrator = (; simulation.u, simulation.p, simulation.t, simulation.dt, alg = simulation.stepper.algo, sol = (; prob = (; f = simulation.stepper.func)))
    ClimaTimeSteppers.step_u!(integrator, simulation.stepper.cache)
    for callback in simulation.callbacks
        callback.condition(simulation.u, simulation.t, integrator) && callback.affect!(integrator)
    end
end

function solve!(simulation::ClimaLandSimulation)
    while simulation.t < simulation.tf
        step!(simulation)
    end
end
