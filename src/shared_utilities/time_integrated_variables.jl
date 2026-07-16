export TimeIntegratedVariable,
    RunningMean, RunningSum, TimeIntegral, time_integrated_tendency!

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

Exponentially-weighted running mean: `dX/dt = (f - X) / τ`. `X` carries the units
of `f` and relaxes toward `f` with e-folding memory `τ` — the smooth analog of
"the mean of `f` over the past `τ`".
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
time, as in the soil water/energy conservation trackers (`∫F_vol_liq_water_dt`,
`∫F_e_dt`).
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

# Fields
$(DocStringExtensions.FIELDS)

# Example
```julia
# trailing ~1-year precipitation total, seeded from a climatology
TimeIntegratedVariable(;
    name = :precip_annual,
    reduction = RunningSum(),
    timescale = FT(365 * 86400),
    compute_instantaneous! = (dst, Y, p, t) -> (@. dst = p.drivers.P_liq),
)
```
"""
struct TimeIntegratedVariable{
    FT <: AbstractFloat,
    R <: AbstractTimeReduction,
    F,
}
    "Name of the prognostic variable (the symbol used in `Y.<component>.<name>`)."
    name::Symbol
    "Reduction kernel: `RunningMean`, `RunningSum`, or `TimeIntegral`."
    reduction::R
    "Memory timescale `τ` in seconds (ignored by `TimeIntegral`)."
    timescale::FT
    "In-place instantaneous quantity: `compute_instantaneous!(dst, Y, p, t)` writes `f` into `dst`. `nothing` for a metadata-only spec (declaration without a tendency)."
    compute_instantaneous!::F
    "Domain the variable lives on (`:surface` or `:subsurface`)."
    domain_name::Symbol
    "Element type of the stored field (defaults to `typeof(timescale)`; set explicitly for structured, e.g. `NamedTuple`-valued, variables)."
    element_type::DataType
end

function TimeIntegratedVariable(;
    name::Symbol,
    reduction::AbstractTimeReduction,
    timescale::FT,
    compute_instantaneous! = nothing,
    domain_name::Symbol = :surface,
    element_type::DataType = typeof(timescale),
) where {FT <: AbstractFloat}
    return TimeIntegratedVariable{
        FT,
        typeof(reduction),
        typeof(compute_instantaneous!),
    }(
        name,
        reduction,
        timescale,
        compute_instantaneous!,
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
prognostic-variable tuples. Only the specs' metadata is used, so the specs may be
built without a `compute_instantaneous!` closure at declaration time.
"""
time_integrated_prognostic_vars(specs) = map(s -> s.name, Tuple(specs))
time_integrated_prognostic_types(specs) = map(s -> s.element_type, Tuple(specs))
time_integrated_prognostic_domain_names(specs) =
    map(s -> s.domain_name, Tuple(specs))

# Convert the instantaneous value already stored in `dst` into the tendency dX/dt
# for the chosen reduction, in place: `dst` enters holding `f` and exits holding
# dX/dt. Written with ClimaCore.RecursiveApply operators (`⊟`, `rdiv`) so that `X`
# may be scalar-valued or have a structured element type (e.g. a NamedTuple-valued
# field, as for the P-model capacities); for scalar `FT` these reduce to the
# ordinary `(f - X)/τ` and `f - X/τ`.
@inline apply_time_reduction!(dst, X, τ, ::RunningMean) =
    (@. dst = rdiv(dst ⊟ X, τ))
@inline apply_time_reduction!(dst, X, τ, ::RunningSum) =
    (@. dst = dst ⊟ rdiv(X, τ))
@inline apply_time_reduction!(dst, X, τ, ::TimeIntegral) = dst

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
