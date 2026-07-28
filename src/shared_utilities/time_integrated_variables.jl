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
    TimeIntegral{FT <: AbstractFloat} <: AbstractTimeReduction

Pure time integral: `dX/dt = f`. Unlike the
running reductions there is no trailing window — `X` accumulates `f` over all of
time.
"""
struct TimeIntegral <: AbstractTimeReduction

"""
    TimeIntegratedVariable

Declaration of a prognostic variable that stores a running reduction
(`reduction`) of an instantaneous quantity `f`.

Because the variable lives in the prognostic state `Y`, it is advanced by the
model's time-stepper. A
model that owns such variables derives its `prognostic_vars` (and the matching
`_types` / `_domain_names`) from a tuple of specs via
[`time_integrated_prognostic_vars`](@ref) and companions, and adds their
tendencies with [`time_integrated_tendency!`](@ref).

`f` (and hence `X`) may be scalar-valued or have a structured element type, e.g. a
`NamedTuple`-valued field like the P-model acclimated capacities. `element_type`
records that element type for the prognostic-state allocation; it defaults to the
type of `timescale` (i.e. `FT`) for the common scalar case.

A spec holds only the declaration, so it is built once by the model that owns it
(from its parameters) and used both to declare the prognostic variables and to add
their tendencies. How `f` is computed is supplied separately, to
[`time_integrated_tendency!`](@ref), since it generally closes over the parent
model — the cache, grid and drivers are not in scope where a component declares
its variables.

# Fields
$(DocStringExtensions.FIELDS)

# Example
```julia
# trailing ~1-year precipitation mean
TimeIntegratedVariable(;
    name = :precip_mean,
    reduction = RunningMean(365*24*3600),
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
    "Element type of the stored field (defaults to `typeof(timescale)`; set explicitly for structured, e.g. `NamedTuple`-valued, variables)."
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

# Write the tendency dX/dt for the chosen reduction into the destination `dY` field
# `dst`, given the instantaneous quantity `f` and current state `X`. Written with
# ClimaCore.RecursiveApply operators (`⊟`, `rdiv`) so `f` and `X` may be scalar-valued
# or have a structured element type (e.g. a NamedTuple-valued field, as for the P-model
# capacities); for scalar `FT` these reduce to the ordinary `(f - X)/τ` and `f - X/τ`.
# `f` may be a `lazy` broadcast, which fuses into the assignment with no allocation.
@inline apply_time_reduction!(dst, f, X, r::RunningMean) =
    (dst = rdiv(f ⊟ X, r.τ))
@inline apply_time_reduction!(dst, f, X, r::RunningSum) =
    (dst = rmul(f,r.τ) ⊟ rdiv(X, r.τ_long))
@inline apply_time_reduction!(dst, f, X, ::TimeIntegral) = (dst = f)
