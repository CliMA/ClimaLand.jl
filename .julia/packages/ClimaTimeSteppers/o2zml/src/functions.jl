import DiffEqBase
export AbstractClimaODEFunction
export ClimaODEFunction, ForwardEulerODEFunction

abstract type AbstractClimaODEFunction <: DiffEqBase.AbstractODEFunction{true} end

"""
    ClimaODEFunction(; T_imp!, [dss!], [cache!], [cache_imp!])
    ClimaODEFunction(; T_exp!, T_lim!, [T_imp!], [lim!], [dss!], [cache!], [cache_imp!])
    ClimaODEFunction(; T_exp_lim!, [T_imp!], [lim!], [dss!], [cache!], [cache_imp!])

Container for all functions used to advance through a timestep:
    - `T_imp!(T_imp, u, p, t)`: sets the implicit tendency
    - `T_exp!(T_exp, u, p, t)`: sets the component of the explicit tendency that
      is not passed through the limiter
    - `T_lim!(T_lim, u, p, t)`: sets the component of the explicit tendency that
      is passed through the limiter
    - `T_exp_lim!(T_exp, T_lim, u, p, t)`: fused alternative to the separate
      functions `T_exp!` and `T_lim!`
    - `lim!(u, p, t, u_ref)`: applies the limiter to every state `u` that has
      been incremented from `u_ref` by the explicit tendency component `T_lim!`
    - `dss!(u, p, t)`: applies direct stiffness summation to every state `u`,
      except for intermediate states generated within the implicit solver
    - `cache!(u, p, t)`: updates the cache `p` to reflect the state `u` before
      the first timestep and on every subsequent timestepping stage
    - `cache_imp!(u, p, t)`: updates the components of the cache `p` that are
      required to evaluate `T_imp!` and its Jacobian within the implicit solver
By default, `lim!`, `dss!`, and `cache!` all do nothing, and `cache_imp!` is
identical to `cache!`. Any of the tendency functions can be set to `nothing` in
order to avoid corresponding allocations in the integrator.
"""
struct ClimaODEFunction{TEL, TL, TE, TI, L, D, C, CI} <: AbstractClimaODEFunction
    T_exp_T_lim!::TEL
    T_lim!::TL
    T_exp!::TE
    T_imp!::TI
    lim!::L
    dss!::D
    cache!::C
    cache_imp!::CI
    function ClimaODEFunction(;
        T_exp_T_lim! = nothing,
        T_lim! = nothing,
        T_exp! = nothing,
        T_imp! = nothing,
        lim! = (u, p, t, u_ref) -> nothing,
        dss! = (u, p, t) -> nothing,
        cache! = (u, p, t) -> nothing,
        cache_imp! = cache!,
    )
        args = (T_exp_T_lim!, T_lim!, T_exp!, T_imp!, lim!, dss!, cache!, cache_imp!)

        if !isnothing(T_exp_T_lim!)
            @assert isnothing(T_exp!) "`T_exp_T_lim!` was passed, `T_exp!` must be `nothing`"
            @assert isnothing(T_lim!) "`T_exp_T_lim!` was passed, `T_lim!` must be `nothing`"
        end
        if !isnothing(T_exp!) && !isnothing(T_lim!)
            @warn "Both `T_exp!` and `T_lim!` are not `nothing`, please use `T_exp_T_lim!` instead."
        end
        return new{typeof.(args)...}(args...)
    end
end

has_T_exp(f::ClimaODEFunction) = !isnothing(f.T_exp!) || !isnothing(f.T_exp_T_lim!)
has_T_lim(f::ClimaODEFunction) = !isnothing(f.lim!) && (!isnothing(f.T_lim!) || !isnothing(f.T_exp_T_lim!))

# Don't wrap a AbstractClimaODEFunction in an ODEFunction (makes ODEProblem work).
DiffEqBase.ODEFunction{iip}(f::AbstractClimaODEFunction) where {iip} = f
DiffEqBase.ODEFunction(f::AbstractClimaODEFunction) = f

"""
    ForwardEulerODEFunction(f; jac_prototype, Wfact, tgrad)

An ODE function wrapper where `f(un, u, p, t, dt)` provides a forward Euler update
```
un .= u .+ dt * f(u, p, t)
```

"""
struct ForwardEulerODEFunction{F, J, W, T} <: DiffEqBase.AbstractODEFunction{true}
    f::F
    jac_prototype::J
    Wfact::W
    tgrad::T
end
ForwardEulerODEFunction(f; jac_prototype = nothing, Wfact = nothing, tgrad = nothing) =
    ForwardEulerODEFunction(f, jac_prototype, Wfact, tgrad)
(f::ForwardEulerODEFunction{F})(un, u, p, t, dt) where {F} = f.f(un, u, p, t, dt)

# Don't wrap a ForwardEulerODEFunction in an ODEFunction.
DiffEqBase.ODEFunction{iip}(f::ForwardEulerODEFunction) where {iip} = f
DiffEqBase.ODEFunction(f::ForwardEulerODEFunction) = f

"""
    OffsetODEFunction(f,α,β,γ,x)

An ODE function wrapper which evaluates `f` with an offset.

Evaluates as
```math
f(u,p,α+β*t) .+ γ .* x
```

It supports the 3, 4, 5, and 6 argument forms.
"""
mutable struct OffsetODEFunction{iip, F, S, A} <: DiffEqBase.AbstractODEFunction{iip}
    f::F
    α::S
    β::S
    γ::S
    x::A
end
function OffsetODEFunction{iip}(f, α, β, γ, x) where {iip}
    α, β, γ = promote(α, β, γ)
    OffsetODEFunction{iip, typeof(f), typeof(γ), typeof(x)}(f, α, β, γ, x)
end
OffsetODEFunction(f::DiffEqBase.AbstractODEFunction{iip}, α, β, γ, x) where {iip} =
    OffsetODEFunction{iip}(f, α, β, γ, x)

function (o::OffsetODEFunction)(u, p, t)
    o.f(u, p, o.α + o.β * t) .+ o.γ .* o.x
end
function (o::OffsetODEFunction)(du, u, p, t)
    o.f(du, u, p, o.α + o.β * t)
    du .+= o.γ .* o.x
end
function (o::OffsetODEFunction)(du, u, p, t, α)
    o.f(du, u, p, o.α + o.β * t, α)
    du .+= α .* o.γ .* o.x
end
function (o::OffsetODEFunction)(du, u, p, t, α, β)
    o.f(du, u, p, o.α + o.β * t, α, β)
    du .+= α .* o.γ .* o.x
end
