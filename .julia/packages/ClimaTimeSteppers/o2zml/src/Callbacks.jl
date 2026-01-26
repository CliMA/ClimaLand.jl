"""
    ClimaTimeSteppers.Callbacks

A suite of callback functions to be used with the ClimaTimeSteppers.jl ODE solvers.
"""
module Callbacks

import ClimaComms, DiffEqBase
import SciMLBase

"""
    ClimaTimeSteppers.Callbacks.initialize!(f!::F, integrator)

Initialize a callback event for callbacks of type `F`. By default this does nothing, but
can be extended for new callback events.
"""
function initialize!(f!, integrator) end

"""
    ClimaTimeSteppers.Callbacks.finalize!(f!::F, integrator)

Finalize a callback event for callbacks of type `F`. By default this does nothing, but
can be extended for new callback events.
"""
function finalize!(f!, integrator) end


export EveryXWallTimeSeconds, EveryXSimulationTime, EveryXSimulationSteps

"""
    EveryXWallTimeSeconds(
        f!,
        Δwt,
        comm_ctx::ClimaComms.AbstractCommsContext;
        atinit=false
    )

Trigger `f!(integrator)` every `Δwt` wallclock seconds.

An [ClimaComms context](https://clima.github.io/ClimaComms.jl/) must be provided to synchronize timing across all ranks.

[`Callbacks.initialize!`](@ref) and [`Callbacks.finalize!`](@ref) can be defined for `f!`.

If `atinit=true`, then `f!(integrator)` will additionally be triggered at initialization,
otherwise the first trigger will be after `Δwt` seconds.
"""
function EveryXWallTimeSeconds(f!, Δwt, comm_ctx::ClimaComms.AbstractCommsContext; atinit = false)
    wt_next = 0.0

    function _initialize(c, u, t, integrator)
        wt = ClimaComms.allreduce(comm_ctx, time(), max)
        wt_next = wt + Δwt
        initialize!(c.affect!, integrator)
        if atinit
            c.affect!(integrator)
        end
    end

    function _finalize(c, u, t, integrator)
        finalize!(c.affect!, integrator)
    end

    function condition(u, t, integrator)
        wt = ClimaComms.allreduce(comm_ctx, time(), max)
        if wt >= wt_next
            while wt >= wt_next
                wt_next += Δwt
            end
            return true
        else
            return false
        end
    end

    if isdefined(DiffEqBase, :finalize!)
        SciMLBase.DiscreteCallback(condition, f!; initialize = _initialize, finalize = _finalize)
    else
        SciMLBase.DiscreteCallback(condition, f!; initialize = _initialize)
    end
end

"""
    EveryXSimulationTime(f!, Δt; atinit=false)

Trigger `f!(integrator)` every `Δt` simulation time.

[`Callbacks.initialize!`](@ref) and [`Callbacks.finalize!`](@ref) can be defined for `f!`.

If `atinit=true`, then `f!` will additionally be triggered at initialization. Otherwise
the first trigger will be after `Δt` simulation time.
"""
function EveryXSimulationTime(f!, Δt; atinit = false)
    t_next = zero(Δt)

    function _initialize(c, u, t, integrator)
        t_next = Δt
        initialize!(c.affect!, integrator)
        if atinit
            c.affect!(integrator)
        end
    end

    function _finalize(c, u, t, integrator)
        finalize!(c.affect!, integrator)
    end


    function condition(u, t, integrator)
        if t >= t_next
            while t >= t_next
                t_next += Δt
            end
            return true
        else
            return false
        end
    end
    if isdefined(DiffEqBase, :finalize!)
        SciMLBase.DiscreteCallback(condition, f!; initialize = _initialize, finalize = _finalize)
    else
        SciMLBase.DiscreteCallback(condition, f!; initialize = _initialize)
    end
end

"""
    EveryXSimulationSteps(f!, Δsteps; atinit=false)

Trigger `f!(integrator)` every `Δsteps` simulation steps.

[`Callbacks.initialize!`](@ref) and [`Callbacks.finalize!`](@ref) can be defined for `f!`.

If `atinit==true`, then `f!` will additionally be triggered at initialization. Otherwise
the first trigger will be after `Δsteps`.
"""
function EveryXSimulationSteps(f!, Δsteps; atinit = false)
    steps = 0
    steps_next = 0

    function _initialize(c, u, t, integrator)
        steps = 0
        steps_next = Δsteps
        initialize!(c.affect!, integrator)
        if atinit
            c.affect!(integrator)
        end
    end

    function _finalize(c, u, t, integrator)
        finalize!(c.affect!, integrator)
    end

    function condition(u, t, integrator)
        steps += 1
        if steps >= steps_next
            steps_next += Δsteps
            return true
        else
            return false
        end
    end

    if isdefined(DiffEqBase, :finalize!)
        SciMLBase.DiscreteCallback(condition, f!; initialize = _initialize, finalize = _finalize)
    else
        SciMLBase.DiscreteCallback(condition, f!; initialize = _initialize)
    end
end

end # module
