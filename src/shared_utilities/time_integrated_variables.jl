export TimeIntegratedVariable,
    RunningMean, RunningSum, TimeIntegral, time_integrated_tendency!

import ClimaCore.RecursiveApply: ⊟, ⊠, rdiv

"""
    AbstractTimeReduction

Reduction kernel applied by a [`TimeIntegratedVariable`](@ref). It sets the ODE
obeyed by the stored running statistic `X` given an instantaneous quantity `f` and
a memory timescale `τ`.
"""
abstract type AbstractTimeReduction end

"""
    RunningMean <: AbstractTimeReduction

Exponentially-weighted running mean: `dX/dt = (f - X) / τ`. `X` carries the units
of `f` and relaxes toward `f` with e-folding memory `τ` — the smooth analog of
"the mean of `f` over the past `τ`".

Explicit-Euler stepping of this ODE is the blend `X ← (1 - Δt/τ) X + (Δt/τ) f`, the
exponential moving average used by the P-model acclimation. The two are identical
for constant `f`; for time-varying `f` they agree to `O((Δt/τ)²)` when `Δt/τ` is
small (the EMA uses `f` at the end of the step, the ODE its value during it).
"""
struct RunningMean <: AbstractTimeReduction end

"""
    RunningSum <: AbstractTimeReduction

Exponentially-weighted running sum over ~`τ`: `dX/dt = f - X / τ`. At steady state
`X = τ f̄`, so for a flux `f` and `τ = 1 yr`, `X` is a trailing annual total (the
smooth analog of "the accumulation of `f` over the past `τ`"). Equal to `τ ×`
[`RunningMean`](@ref). The trailing window (the decay term `X/τ`) is what makes it
"running"; contrast [`TimeIntegral`](@ref), which has none.
"""
struct RunningSum <: AbstractTimeReduction end

"""
    TimeIntegral <: AbstractTimeReduction

Pure time integral with no memory decay: `dX/dt = f` (`τ` is ignored). Unlike the
running reductions there is no trailing window — `X` accumulates `f` over all of
time.
"""
struct TimeIntegral <: AbstractTimeReduction end

"""
    TimeIntegratedVariable

Specification of a prognostic variable that stores a running reduction
(`reduction`) of an instantaneous quantity `f` over a memory timescale
`timescale`.

Because the variable lives in the prognostic state `Y`, it is advanced by the
model's time-stepper on every stage — smoothly, with no periodic callback and no
year-boundary discontinuity — and it is serialized with checkpoints/restarts. A
model that owns such variables derives its `prognostic_vars` (and the matching
`_types` / `_domain_names`) from a tuple of specs via
[`time_integrated_prognostic_vars`](@ref) and companions, and adds their
tendencies with [`time_integrated_tendency!`](@ref).

`f` (and hence `X`) may be scalar-valued or have a structured element type, e.g. a
`NamedTuple`-valued field like the P-model acclimated capacities. `element_type`
records that element type for the prognostic-state allocation; it defaults to the
type of `timescale` (i.e. `FT`) for the common scalar case.

A spec is usually built in two steps: the model declares the metadata (`name`,
`reduction`, `timescale`, ...) where the cache and grid are not yet in scope — this
drives the `prognostic_vars` helpers — and later attaches the runtime
`compute_instantaneous` (and optional `weight`) closures via the copy constructor
`TimeIntegratedVariable(base; compute_instantaneous, weight)`.

# Fields
$(DocStringExtensions.FIELDS)

# Example
```julia
# trailing ~1-year precipitation total, seeded from a climatology
TimeIntegratedVariable(;
    name = :precip_annual,
    reduction = RunningSum(),
    timescale = FT(365 * 86400),
    compute_instantaneous = (Y, p, t) -> p.drivers.P_liq,
)
```
"""
struct TimeIntegratedVariable{
    FT <: AbstractFloat,
    R <: AbstractTimeReduction,
    F,
    W,
}
    "Name of the prognostic variable (the symbol used in `Y.<component>.<name>`, e.g. `Y.canopy.biomass.A0_annual`)."
    name::Symbol
    "Reduction kernel: `RunningMean`, `RunningSum`, or `TimeIntegral`."
    reduction::R
    "Memory timescale `τ` in seconds (ignored by `TimeIntegral`)."
    timescale::FT
    "Instantaneous quantity `f`: `compute_instantaneous(Y, p, t)` returns the field (or `lazy` broadcast) whose running reduction is stored. `nothing` for a metadata-only spec (declaration without a tendency)."
    compute_instantaneous::F
    "Optional weight `weight(Y, p, t)` returning a scalar field that multiplies the reduction rate (e.g. a solar-noon window whose daily mean is 1, so it reshapes the sub-daily forcing without changing `τ`). `nothing` applies no weighting."
    weight::W
    "Domain the variable lives on (`:surface` or `:subsurface`)."
    domain_name::Symbol
    "Element type of the stored field (defaults to `typeof(timescale)`; set explicitly for structured, e.g. `NamedTuple`-valued, variables)."
    element_type::DataType
end

function TimeIntegratedVariable(;
    name::Symbol,
    reduction::AbstractTimeReduction,
    timescale::FT,
    compute_instantaneous = nothing,
    weight = nothing,
    domain_name::Symbol = :surface,
    element_type::DataType = typeof(timescale),
) where {FT <: AbstractFloat}
    return TimeIntegratedVariable{
        FT,
        typeof(reduction),
        typeof(compute_instantaneous),
        typeof(weight),
    }(
        name,
        reduction,
        timescale,
        compute_instantaneous,
        weight,
        domain_name,
        element_type,
    )
end

"""
    TimeIntegratedVariable(base::TimeIntegratedVariable; compute_instantaneous, weight)

Copy `base`'s declaration (name, reduction, timescale, domain, element type) and
attach the runtime `compute_instantaneous` (and optional `weight`) closures. This
lets a model declare its specs' metadata once — for the `prognostic_vars` helpers —
and supply the tendency kernels later, where the cache and grid are in scope.
"""
function TimeIntegratedVariable(
    base::TimeIntegratedVariable;
    compute_instantaneous = base.compute_instantaneous,
    weight = base.weight,
)
    return TimeIntegratedVariable(;
        name = base.name,
        reduction = base.reduction,
        timescale = base.timescale,
        compute_instantaneous,
        weight,
        domain_name = base.domain_name,
        element_type = base.element_type,
    )
end

"""
    time_integrated_prognostic_vars(specs)
    time_integrated_prognostic_types(specs)
    time_integrated_prognostic_domain_names(specs)

Expand a tuple of [`TimeIntegratedVariable`](@ref)s into the parallel tuples
expected by [`prognostic_vars`](@ref), [`prognostic_types`](@ref), and
[`prognostic_domain_names`](@ref). A model concatenates these with its own
prognostic-variable tuples. Only the specs' metadata is used, so the specs may be
built without a `compute_instantaneous` closure at declaration time.
"""
time_integrated_prognostic_vars(specs) = map(s -> s.name, Tuple(specs))
time_integrated_prognostic_types(specs) = map(s -> s.element_type, Tuple(specs))
time_integrated_prognostic_domain_names(specs) =
    map(s -> s.domain_name, Tuple(specs))

# Write the tendency dX/dt for the chosen reduction into the destination `dY` field
# `dst`, given the instantaneous quantity `f` and current state `X`. Written with
# ClimaCore.RecursiveApply operators (`⊟`, `rdiv`) so `f` and `X` may be scalar-valued
# or have a structured element type (e.g. a NamedTuple-valued field, as for the P-model
# capacities); for scalar `FT` these reduce to the ordinary `(f - X)/τ` and `f - X/τ`.
# `f` may be a `lazy` broadcast, which fuses into the assignment with no allocation.
@inline apply_time_reduction!(dst, f, X, τ, ::RunningMean) =
    (@. dst = rdiv(f ⊟ X, τ))
@inline apply_time_reduction!(dst, f, X, τ, ::RunningSum) =
    (@. dst = f ⊟ rdiv(X, τ))
@inline apply_time_reduction!(dst, f, X, τ, ::TimeIntegral) = (@. dst = f)

# Multiply the just-computed tendency `dst` by an optional weight `w` (a scalar field
# returned by the spec's `weight`). Dispatch on `nothing` keeps the common unweighted
# path branch-free. Written with `⊠` so `dst` may be scalar- or structured-valued.
@inline _apply_tendency_weight!(dst, ::Nothing, Y, p, t) = nothing
@inline function _apply_tendency_weight!(dst, weight, Y, p, t)
    w = weight(Y, p, t)
    @. dst = dst ⊠ w
    return nothing
end

"""
    time_integrated_tendency!(dY_component, Y_component, Y, p, t, specs)

Add the tendencies of a collection of [`TimeIntegratedVariable`](@ref)s to
`dY_component`, the sub-FieldVector of `dY` owned by the calling model (e.g.
`dY.canopy.biomass`). `Y_component` is the matching sub-FieldVector of `Y`; the
full `Y` and cache `p` are forwarded to each spec's `compute_instantaneous` so that
`f` may depend on drivers or other state.

For each variable the instantaneous `f = compute_instantaneous(Y, p, t)` (a field or
`lazy` broadcast, so no field is allocated) is reduced into the destination tendency
field, which a spec's (daily-mean-1) `weight` then optionally rescales.
"""
function time_integrated_tendency!(dY_c, Y_c, Y, p, t, specs)
    foreach(specs) do s
        dst = getproperty(dY_c, s.name)
        X = getproperty(Y_c, s.name)
        f = s.compute_instantaneous(Y, p, t)
        apply_time_reduction!(dst, f, X, s.timescale, s.reduction)
        _apply_tendency_weight!(dst, s.weight, Y, p, t)
    end
    return nothing
end
