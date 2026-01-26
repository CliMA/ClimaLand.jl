import DataStructures

"""
    DistributedODEIntegrator <: AbstractODEIntegrator

A simplified variant of [`ODEIntegrator`](https://github.com/SciML/OrdinaryDiffEq.jl/blob/6ec5a55bda26efae596bf99bea1a1d729636f412/src/integrators/type.jl#L77-L123).

"""
mutable struct DistributedODEIntegrator{
    algType,
    uType,
    tType,
    pType,
    heapType,
    tstopsType,
    saveatType,
    callbackType,
    cacheType,
    solType,
} <: DiffEqBase.AbstractODEIntegrator{algType, true, uType, tType}
    alg::algType
    u::uType
    p::pType
    t::tType
    dt::tType
    _dt::tType # argument to __init used to set dt in step!
    dtchangeable::Bool
    tstops::heapType
    _tstops::tstopsType # argument to __init used as default argument to reinit!
    saveat::heapType
    _saveat::saveatType # argument to __init used as default argument to reinit!
    step::Int
    stepstop::Int
    callback::callbackType
    advance_to_tstop::Bool
    u_modified::Bool # not used; field is required for compatibility with
    # DiffEqBase.initialize! and DiffEqBase.finalize!
    cache::cacheType
    sol::solType
    tdir::tType # see https://docs.sciml.ai/DiffEqCallbacks/stable/output_saving/#DiffEqCallbacks.SavingCallback
end

"""
    SavedValues{tType<:Real, savevalType}

A struct used to save values of the time in `t::Vector{tType}` and
additional values in `saveval::Vector{savevalType}`.

From DiffEqCallbacks.
"""
struct SavedValues{tType, savevalType}
    t::Vector{tType}
    saveval::Vector{savevalType}
end

"""
    SavedValues(tType::DataType, savevalType::DataType)

Return `SavedValues{tType, savevalType}` with empty storage vectors.

From DiffEqCallbacks.
"""
function SavedValues(::Type{tType}, ::Type{savevalType}) where {tType, savevalType}
    SavedValues{tType, savevalType}(Vector{tType}(), Vector{savevalType}())
end


# helper function for setting up min/max heaps for tstops and saveat
function tstops_and_saveat_heaps(t0, tf, tstops, saveat = [])
    # We promote to a common type to ensure that t0 and tf have the same type
    FT = typeof(first(promote(t0, tf)))
    ordering = tf > t0 ? DataStructures.FasterForward : DataStructures.FasterReverse

    # ensure that tstops includes tf and only has values ahead of t0
    tstops = [filter(t -> t0 < t < tf || tf < t < t0, tstops)..., tf]
    tstops = DataStructures.BinaryHeap{FT, ordering}(tstops)

    isnothing(saveat) && (saveat = [t0, tf])

    saveat = DataStructures.BinaryHeap{FT, ordering}(collect(saveat))

    return tstops, saveat
end

compute_tdir(ts) = ts[1] > ts[end] ? sign(ts[end] - ts[1]) : oneunit(ts[1])

# called by DiffEqBase.init and DiffEqBase.solve
function DiffEqBase.__init(
    prob::DiffEqBase.AbstractODEProblem,
    alg::DistributedODEAlgorithm,
    args...;
    dt,
    tstops = (),
    saveat = nothing,
    save_everystep = false,
    callback = nothing,
    advance_to_tstop = false,
    save_func = (u, t) -> copy(u),   # custom kwarg
    dtchangeable = true,             # custom kwarg
    stepstop = -1,                   # custom kwarg
    tdir = compute_tdir(prob.tspan), #
    kwargs...,
)
    (; u0, p) = prob
    t0, tf = prob.tspan
    t0, tf, dt = promote(t0, tf, dt)

    # We need zero(oneunit()) because there's no zerounit
    dt > zero(oneunit(dt)) || error("dt must be positive")
    _dt = dt
    dt = tf > t0 ? dt : -dt

    _tstops = tstops
    _saveat = saveat
    tstops, saveat = tstops_and_saveat_heaps(t0, tf, tstops, saveat)

    sol = DiffEqBase.build_solution(prob, alg, typeof(t0)[], typeof(save_func(u0, t0))[])
    saving_callback = NonInterpolatingSavingCallback(save_func, SavedValues(sol.t, sol.u), save_everystep)
    callback = DiffEqBase.CallbackSet(callback, saving_callback)
    isempty(callback.continuous_callbacks) || error("Continuous callbacks are not supported")

    integrator = DistributedODEIntegrator(
        alg,
        u0,
        p,
        t0,
        dt,
        _dt,
        dtchangeable,
        tstops,
        _tstops,
        saveat,
        _saveat,
        0,
        stepstop,
        callback,
        advance_to_tstop,
        false,
        init_cache(prob, alg; dt, kwargs...),
        sol,
        tdir,
    )
    if prob.f isa ClimaODEFunction
        (; cache!) = prob.f
        isnothing(cache!) || cache!(u0, p, t0)
    end
    DiffEqBase.initialize!(callback, u0, t0, integrator)
    return integrator
end

DiffEqBase.has_reinit(integrator::DistributedODEIntegrator) = true
function DiffEqBase.reinit!(
    integrator::DistributedODEIntegrator,
    u0 = integrator.sol.prob.u0;
    t0 = integrator.sol.prob.tspan[1],
    tf = integrator.sol.prob.tspan[2],
    erase_sol = true,
    tstops = integrator._tstops,
    saveat = integrator._saveat,
    reinit_callbacks = true,
)
    integrator.u .= u0
    integrator.t = t0
    integrator.tstops, integrator.saveat = tstops_and_saveat_heaps(t0, tf, tstops, saveat)
    integrator.step = 0
    if erase_sol
        resize!(integrator.sol.t, 0)
        resize!(integrator.sol.u, 0)
    end
    if reinit_callbacks
        DiffEqBase.initialize!(integrator.callback, u0, t0, integrator)
    else # always reinit the saving callback so that t0 can be saved if needed
        saving_callback = integrator.callback.discrete_callbacks[end]
        DiffEqBase.initialize!(saving_callback, u0, t0, integrator)
    end
end

# called by DiffEqBase.solve
function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem, alg::DistributedODEAlgorithm, args...; kwargs...)
    integrator = DiffEqBase.__init(prob, alg, args...; kwargs...)
    DiffEqBase.solve!(integrator)
end

# either called directly (after init), or by DiffEqBase.solve (via __solve)
NVTX.@annotate function DiffEqBase.solve!(integrator::DistributedODEIntegrator)
    while !isempty(integrator.tstops) && integrator.step != integrator.stepstop
        __step!(integrator)
    end
    DiffEqBase.finalize!(integrator.callback, integrator.u, integrator.t, integrator)
    return integrator.sol
end

function DiffEqBase.step!(integrator::DistributedODEIntegrator)
    if integrator.advance_to_tstop
        tstop = first(integrator.tstops)
        while !reached_tstop(integrator, tstop)
            __step!(integrator)
        end
    else
        __step!(integrator)
    end
end

function DiffEqBase.step!(integrator::DistributedODEIntegrator, dt, stop_at_tdt = false)
    # OridinaryDiffEq lets dt be negative if tdir is -1, but that's inconsistent
    dt <= zero(dt) && error("dt must be positive")
    stop_at_tdt && !integrator.dtchangeable && error("Cannot stop at t + dt if dtchangeable is false")
    t_plus_dt = integrator.t + tdir(integrator) * dt
    stop_at_tdt && DiffEqBase.add_tstop!(integrator, t_plus_dt)
    while !reached_tstop(integrator, t_plus_dt, stop_at_tdt)
        __step!(integrator)
    end
end

# helper functions for dealing with time-reversed integrators in the same way
# that OrdinaryDiffEq.jl does
tdir(integrator) = integrator.tstops.ordering isa DataStructures.FasterForward ? 1 : -1
is_past_t(integrator, t) = tdir(integrator) * (t - integrator.t) < zero(integrator.t)
reached_tstop(integrator, tstop, stop_at_tstop = integrator.dtchangeable) =
    integrator.t == tstop || (!stop_at_tstop && is_past_t(integrator, tstop))


@inline unrolled_foreach(::Tuple{}, integrator) = nothing
@inline unrolled_foreach(callback, integrator) =
    callback.condition(integrator.u, integrator.t, integrator) ? callback.affect!(integrator) : nothing
@inline unrolled_foreach(discrete_callbacks::Tuple{Any}, integrator) =
    unrolled_foreach(first(discrete_callbacks), integrator)
@inline function unrolled_foreach(discrete_callbacks::Tuple, integrator)
    unrolled_foreach(first(discrete_callbacks), integrator)
    unrolled_foreach(Base.tail(discrete_callbacks), integrator)
end

function __step!(integrator)
    (; _dt, dtchangeable, tstops) = integrator

    # update step and dt before incrementing u; if dt is changeable and there is
    # a tstop within dt, reduce dt to tstop - t
    integrator.step += 1
    integrator.dt =
        !isempty(tstops) && dtchangeable ? tdir(integrator) * min(_dt, abs(first(tstops) - integrator.t)) :
        tdir(integrator) * _dt
    step_u!(integrator)

    # increment t by dt, rounding to the first tstop if that is roughly
    # equivalent up to machine precision; the specific bound of 100 * eps...
    # is taken from OrdinaryDiffEq.jl
    t_plus_dt = integrator.t + integrator.dt
    t_unit = oneunit(integrator.t)
    max_t_error = 100 * eps(float(integrator.t / t_unit)) * float(t_unit)
    integrator.t =
        !isempty(tstops) && abs(float(first(tstops)) - float(t_plus_dt)) < max_t_error ? first(tstops) : t_plus_dt

    # apply callbacks
    discrete_callbacks = integrator.callback.discrete_callbacks
    unrolled_foreach(discrete_callbacks, integrator)

    # remove tstops that were just reached
    while !isempty(tstops) && reached_tstop(integrator, first(tstops))
        pop!(tstops)
    end
end

# solvers need to define this interface
NVTX.@annotate function step_u!(integrator)
    step_u!(integrator, integrator.cache)
end

DiffEqBase.get_dt(integrator::DistributedODEIntegrator) = integrator._dt
function set_dt!(integrator::DistributedODEIntegrator, dt)
    # TODO: figure out interface for recomputing other objects (linear operators, etc)
    dt <= zero(dt) && error("dt must be positive")
    integrator._dt = dt
end

function DiffEqBase.add_tstop!(integrator::DistributedODEIntegrator, t)
    is_past_t(integrator, t) && error("Cannot add a tstop at $t because that is behind the current \
                                       integrator time $(integrator.t)")
    push!(integrator.tstops, t)
end

function DiffEqBase.add_saveat!(integrator::DistributedODEIntegrator, t)
    is_past_t(integrator, t) && error("Cannot add a saveat point at $t because that is behind the \
                                       current integrator time $(integrator.t)")
    push!(integrator.saveat, t)
end

# not sure what this should do?
# defined as default initialize: https://github.com/SciML/DiffEqBase.jl/blob/master/src/callbacks.jl#L3
DiffEqBase.u_modified!(i::DistributedODEIntegrator, bool) = nothing

# this is roughly based on SavingCallback from DiffEqCallbacks, except that it
# doesn't interpolate; instead it will save the first step after
function NonInterpolatingSavingCallback(save_func, saved_values, save_everystep)
    if save_everystep
        condition = (u, t, integrator) -> true
    else
        condition =
            (u, t, integrator) -> begin
                cond = false
                while !isempty(integrator.saveat) &&
                    (first(integrator.saveat) == integrator.t || is_past_t(integrator, first(integrator.saveat)))
                    cond = true
                    pop!(integrator.saveat)
                end
                return cond
            end
    end
    function affect!(integrator)
        push!(saved_values.t, integrator.t)
        push!(saved_values.saveval, save_func(integrator.u, integrator.t))
    end
    initialize(cb, u, t, integrator) = condition(u, t, integrator) && affect!(integrator)
    finalize(cb, u, t, integrator) = !save_everystep && !isempty(integrator.saveat) && affect!(integrator)
    SciMLBase.DiscreteCallback(condition, affect!; initialize, finalize)
end
