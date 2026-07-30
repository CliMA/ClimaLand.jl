export TimeIntegratedVariable, RunningMean, RunningSum, TimeIntegral

import ClimaCore.RecursiveApply: ⊟, ⊠, rdiv

"""
    AbstractTimeReduction

Reduction kernel applied by a [`TimeIntegratedVariable`](@ref). It sets the ODE
obeyed by the stored running statistic `X` given an instantaneous quantity `f` and
a memory timescale `τ`.
"""
abstract type AbstractTimeReduction end

# Reductions are passed to `apply_time_reduction` inside `@.` expressions over fields.
Base.broadcastable(r::AbstractTimeReduction) = tuple(r)

"""
    RunningMean{FT <: AbstractFloat} <: AbstractTimeReduction

Exponentially-weighted running mean: `dX/dt = (f - X) / τ`. `X` carries the units
of `f` and relaxes toward `f` with e-folding memory `τ` — the smooth analog of
"the mean of `f` over the past `τ`".

Explicit-Euler stepping of this ODE is the blend `X ← (1 - Δt/τ) X + (Δt/τ) f`, the
exponential moving average used by the P-model acclimation. The two are identical
for constant `f`; for time-varying `f` they agree to `O((Δt/τ)²)` when `Δt/τ` is
small (the EMA uses `f` at the end of the step, the ODE its value during it).
"""
struct RunningMean{FT <: AbstractFloat} <: AbstractTimeReduction
    τ::FT
end

"""
    RunningSum{FT <: AbstractFloat} <: AbstractTimeReduction

Exponentially-weighted running sum over ~`τ`: `dX/dt = (fτ - X) / τ_long`. At steady state
`X = τ f̄`, so for a flux `f` and `τ = 1 yr`, `X` is a trailing annual total (the
smooth analog of "the accumulation of `f` over the past `τ`"). Equal to `τ ×`
[`RunningMean`](@ref). The use of `τ_long > τ` is meant to reduce the aliasing of seasonality.
"""
struct RunningSum{FT <: AbstractFloat} <: AbstractTimeReduction
    τ::FT
    τ_long::FT
end

"""
    TimeIntegral <: AbstractTimeReduction

Pure time integral: `dX/dt = f`. Unlike the
running reductions there is no trailing window — `X` accumulates `f` over all of
time.
"""
struct TimeIntegral <: AbstractTimeReduction end

"""
    TimeIntegratedVariable

Declaration of a prognostic variable that stores a running reduction
(`reduction`) of an instantaneous quantity `f`.

Because the variable lives in the prognostic state `Y`, it is advanced by the
model's time-stepper. A model that owns such variables stores them (see
[`time_integrated_variables`](@ref)) and derives its `prognostic_vars` (and the
matching `_types` / `_domain_names`) from them via
[`time_integrated_prognostic_vars`](@ref) and companions. Its
`make_compute_exp_tendency` writes each tendency with `apply_time_reduction`.

`f` (and hence `X`) may be scalar-valued or have a structured element type, e.g. a
`NamedTuple`-valued field like the P-model acclimated capacities. `element_type`
records that element type for the prognostic-state allocation; it defaults to `FT`
for the common scalar case.

A spec holds only the declaration, so it is built once by the model that owns it
(from its parameters) and stored in it. How `f` is computed is not part of the
declaration: it is written where the tendency is, since it generally needs the
cache, grid and drivers, which are not in scope where a component declares its
variables.

# Fields
$(DocStringExtensions.FIELDS)

# Example
```julia
# trailing ~1-year precipitation mean
TimeIntegratedVariable{FT}(;
    name = :precip_mean,
    reduction = RunningMean(FT(365 * 24 * 3600)),
)
```
"""
struct TimeIntegratedVariable{FT <: AbstractFloat, R <: AbstractTimeReduction}
    "Name of the prognostic variable (the symbol used in `Y.<component>.<name>`, e.g. `Y.canopy.biomass.A0_annual`)."
    name::Symbol
    "Reduction kernel: `RunningMean`, `RunningSum`, or `TimeIntegral`."
    reduction::R
    "Domain the variable lives on (`:surface` or `:subsurface`)."
    domain_name::Symbol
    "Element type of the stored field (defaults to `FT`; set explicitly for structured, e.g. `NamedTuple`-valued, variables)."
    element_type::DataType
end

function TimeIntegratedVariable{FT}(;
    name::Symbol,
    reduction::AbstractTimeReduction,
    domain_name::Symbol = :surface,
    element_type::DataType = FT,
) where {FT <: AbstractFloat}
    return TimeIntegratedVariable{FT, typeof(reduction)}(
        name,
        reduction,
        domain_name,
        element_type,
    )
end

"""
    time_integrated_variables(specs::TimeIntegratedVariable...)

Collect the [`TimeIntegratedVariable`](@ref)s a model owns into a NamedTuple keyed
by their `name`s, the form a model stores them in: the declaration helpers below
expand it into the prognostic-variable tuples, and the model's tendency reaches a
single spec as `model.time_integrated_vars.<name>`.
"""
time_integrated_variables(specs::TimeIntegratedVariable...) =
    NamedTuple{map(s -> s.name, specs)}(specs)

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
time_integrated_prognostic_types(specs) = map(s -> s.element_type, Tuple(specs))
time_integrated_prognostic_domain_names(specs) =
    map(s -> s.domain_name, Tuple(specs))

# Return the tendency dX/dt of the chosen reduction, given the instantaneous quantity
# `f` and the current state `X`. Meant to be broadcast into the destination `dY` field:
#     @. dY.component.name = apply_time_reduction(f, Y.component.name, r)
# Written with ClimaCore.RecursiveApply operators (`⊟`, `⊠`, `rdiv`) so `f` and `X` may be
# scalar-valued or have a structured element type (e.g. a NamedTuple-valued field, as for
# the P-model capacities); for scalar `FT` these are the ordinary `(f - X)/τ` and
# `(f τ - X)/τ_long`.
@inline apply_time_reduction(f, X, r::RunningMean) = rdiv(f ⊟ X, r.τ)
@inline apply_time_reduction(f, X, r::RunningSum) = rdiv(f ⊠ r.τ ⊟ X, r.τ_long)
@inline apply_time_reduction(f, X, ::TimeIntegral) = f
