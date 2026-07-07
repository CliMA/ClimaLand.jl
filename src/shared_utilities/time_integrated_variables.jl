export TimeIntegratedVariable,
    RunningMean, RunningIntegral, RunningSum, time_integrated_tendency!

import ClimaCore.RecursiveApply: ⊟, rdiv

"""
    AbstractTimeReduction

Reduction kernel applied by a [`TimeIntegratedVariable`](@ref). It sets the ODE
obeyed by the stored state `X` given an instantaneous quantity `f` and a memory
timescale `τ`.
"""
abstract type AbstractTimeReduction end

"""
    RunningMean <: AbstractTimeReduction

Exponentially-weighted running mean: `dX/dt = (f - X) / τ`. `X` carries the
units of `f` and relaxes toward `f` with e-folding memory `τ` (the smooth analog
of "the mean of `f` over the past `τ`"). Stepped explicitly with step `Δt` this
is identical to the blend `X ← (1 - Δt/τ) X + (Δt/τ) f` used by the P-model
acclimation EMA.
"""
struct RunningMean <: AbstractTimeReduction end

"""
    RunningIntegral <: AbstractTimeReduction

Exponentially-weighted running integral over ~`τ`: `dX/dt = f - X / τ`. At
steady state `X = τ f̄`, so for a flux `f` and `τ = 1 yr`, `X` is an annual total
(the smooth analog of "the accumulation of `f` over the past `τ`"). Equal to
`τ ×` [`RunningMean`](@ref).
"""
struct RunningIntegral <: AbstractTimeReduction end

"""
    RunningSum <: AbstractTimeReduction

Pure accumulator with no memory decay: `dX/dt = f` (`τ` is ignored). This is the
reduction of the soil water/energy conservation trackers (`∫F_vol_liq_water_dt`,
`∫F_e_dt`).
"""
struct RunningSum <: AbstractTimeReduction end

"""
    TimeIntegratedVariable{FT <: AbstractFloat, R <: AbstractTimeReduction, F}

Specification of a prognostic variable that stores a trailing-window reduction
(`reduction`) of an instantaneous quantity over a memory timescale `timescale`.

Because the variable lives in the prognostic state `Y`, it is advanced by the
model's time-stepper on every stage — smoothly, with no periodic callback and no
year-boundary discontinuity — and it is serialized with checkpoints/restarts. A
model that owns such variables declares them among its prognostic variables via
[`time_integrated_prognostic_vars`](@ref) (and the `_types` / `_domain_names`
companions), and adds their tendencies with [`time_integrated_tendency!`](@ref).

# Fields
$(DocStringExtensions.FIELDS)

# Example
```julia
# cumulative precipitation over ~the past year, seeded from a climatology
TimeIntegratedVariable(;
    name = :precip_annual,
    reduction = RunningIntegral(),
    timescale = FT(365 * 86400),
    compute_instantaneous! = (dst, Y, p, t) -> (@. dst = p.drivers.P_liq),
)
```
"""
Base.@kwdef struct TimeIntegratedVariable{
    FT <: AbstractFloat,
    R <: AbstractTimeReduction,
    F,
}
    "Name of the prognostic variable (the symbol used in `Y.<component>`)."
    name::Symbol
    "Reduction kernel: `RunningMean`, `RunningIntegral`, or `RunningSum`."
    reduction::R
    "Memory timescale `τ` in seconds (ignored by `RunningSum`)."
    timescale::FT
    "In-place instantaneous quantity: `compute_instantaneous!(dst, Y, p, t)` writes `f` into the field `dst`."
    compute_instantaneous!::F
    "Domain the variable lives on (`:surface` or `:subsurface`)."
    domain_name::Symbol = :surface
end

"""
    time_integrated_prognostic_vars(specs)
    time_integrated_prognostic_types(specs)
    time_integrated_prognostic_domain_names(specs)

Expand a tuple of [`TimeIntegratedVariable`](@ref)s into the parallel tuples
expected by [`prognostic_vars`](@ref), [`prognostic_types`](@ref), and
[`prognostic_domain_names`](@ref). A model concatenates these with its own
prognostic-variable tuples.
"""
time_integrated_prognostic_vars(specs) = map(s -> s.name, Tuple(specs))
time_integrated_prognostic_types(specs) =
    map(s -> typeof(s.timescale), Tuple(specs))
time_integrated_prognostic_domain_names(specs) =
    map(s -> s.domain_name, Tuple(specs))

# Convert the instantaneous value already stored in `dst` into the tendency for
# the chosen reduction, in place. `dst` enters holding `f`; on exit it holds dX/dt.
# Written with ClimaCore.RecursiveApply operators (`⊟`, `rdiv`) so a variable may
# be scalar-valued or have a structured element type (e.g. a NamedTuple-valued
# field): for scalar `FT` these are bit-identical to `(dst - X)/τ` and `dst - X/τ`.
@inline apply_time_reduction!(dst, X, τ, ::RunningMean) =
    (@. dst = rdiv(dst ⊟ X, τ))
@inline apply_time_reduction!(dst, X, τ, ::RunningIntegral) =
    (@. dst = dst ⊟ rdiv(X, τ))
@inline apply_time_reduction!(dst, X, τ, ::RunningSum) = dst

"""
    time_integrated_tendency!(dY_component, Y_component, Y, p, t, specs)

Add the tendencies of a collection of [`TimeIntegratedVariable`](@ref)s to
`dY_component`, the sub-FieldVector of `dY` owned by the calling model (e.g.
`dY.canopy.biomass`). `Y_component` is the matching sub-FieldVector of `Y`; the
full `Y` and cache `p` are forwarded to each spec's `compute_instantaneous!` so
that `f` may depend on drivers or other state.

For each variable the instantaneous value `f` is written into the destination
tendency field (reused as scratch), then converted in place to the
reduction-specific rate, so no additional field is allocated.
"""
function time_integrated_tendency!(dY_c, Y_c, Y, p, t, specs)
    foreach(specs) do s
        dst = getproperty(dY_c, s.name)
        X = getproperty(Y_c, s.name)
        s.compute_instantaneous!(dst, Y, p, t)
        apply_time_reduction!(dst, X, s.timescale, s.reduction)
    end
    return nothing
end
