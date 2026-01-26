module ClimaTimeSteppersBenchmarkToolsExt


import StatsBase
import SciMLBase
import PrettyTables
import OrderedCollections
import LinearAlgebra as LA
import ClimaTimeSteppers as CTS
using CUDA
import ClimaComms

import BenchmarkTools
Base.:*(t::BenchmarkTools.Trial, n::Int) = n * t
Base.:*(n::Int, t::BenchmarkTools.Trial) =
    BenchmarkTools.Trial(t.params, t.times .* n, t.gctimes .* n, t.memory * n, t.allocs * n)

include("benchmark_utils.jl")

function get_n_calls_per_step(integrator::CTS.DistributedODEIntegrator)
    (; alg) = integrator
    if alg isa CTS.RosenbrockAlgorithm
        return n_calls_per_step(alg)
    end
    (; newtons_method, name) = alg
    (; max_iters) = newtons_method
    n_calls_per_step(name, max_iters)
end

# TODO: generalize
n_calls_per_step(::CTS.ARS343, max_newton_iters) = Dict(
    "Wfact" => 3 * max_newton_iters,
    "ldiv!" => 3 * max_newton_iters,
    "T_imp!" => 3 * max_newton_iters,
    "T_exp_T_lim!" => 4,
    "lim!" => 4,
    "dss!" => 4,
    "cache!" => 3,
    "cache_imp!" => 4,
    "step!" => 1,
)
function n_calls_per_step(alg::CTS.RosenbrockAlgorithm)
    return Dict(
        "Wfact" => 1,
        "ldiv!" => CTS.n_stages(alg.tableau),
        "T_imp!" => CTS.n_stages(alg.tableau),
        "T_exp_T_lim!" => CTS.n_stages(alg.tableau),
        "lim!" => 0,
        "dss!" => CTS.n_stages(alg.tableau),
        "cache!" => 0,
        "cache_imp!" => CTS.n_stages(alg.tableau),
        "step!" => 1,
    )
end

function maybe_push!(trials₀, name, f!, args, kwargs, only)
    if isnothing(only) || name in only
        trials₀[name] = get_trial(f!, args, name; kwargs...)
    end
end

const allowed_names = ["Wfact", "ldiv!", "T_imp!", "T_exp_T_lim!", "lim!", "dss!", "cache!", "cache_imp!", "step!"]

"""
    benchmark_step(
        integrator::DistributedODEIntegrator,
        device::ClimaComms.AbstractDevice;
        with_cu_prof::Symbol = :bfrofile, # [:bprofile, :profile]
        trace::Bool = false,
        crop::Bool = false,
        hcrop::Union{Nothing, Int} = nothing,
        only::Union{Nothing, Vector{String}} = nothing,
    )

Benchmark a DistributedODEIntegrator given:
 - `integrator` the `DistributedODEIntegrator`.
 - `device` the `ClimaComms` device.
 - `with_cu_prof`, `:profile` or `:bprofile`, to call `CUDA.@profile` or `CUDA.@bprofile` respectively.
 - `trace`, Bool passed to `CUDA.@profile` (see CUDA docs)
 - `crop`, Bool indicating whether or not to crop the `CUDA.@profile` printed table.
 - `hcrop`, Number of horizontal characters to include in the table before cropping.
 - `only, list of functions to benchmarks (benchmark all by default)

`only` may contain:
 - "Wfact"
 - "ldiv!"
 - "T_imp!"
 - "T_exp_T_lim!"
 - "lim!"
 - "dss!"
 - "cache!"
 - "cache_imp!"
 - "step!"
"""
function CTS.benchmark_step(
    integrator::CTS.DistributedODEIntegrator,
    device::ClimaComms.AbstractDevice;
    with_cu_prof::Symbol = :bprofile,
    trace::Bool = false,
    crop::Bool = false,
    hcrop::Union{Nothing, Int} = nothing,
    only::Union{Nothing, Vector{String}} = nothing,
)
    (; u, p, t, dt, sol, alg) = integrator
    (; f) = sol.prob
    if f isa CTS.ClimaODEFunction
        if !isnothing(only)
            @assert all(x -> x in allowed_names, only) "Allowed names in `only` are: $allowed_names"
        end

        X = similar(u)
        Xlim = similar(u)
        @. X = u
        @. Xlim = u
        trials₀ = OrderedCollections.OrderedDict()
        kwargs = (; device, with_cu_prof, trace, crop, hcrop)
#! format: off
        if !isnothing(implicit_fun(integrator))
            W = get_W(integrator)
            maybe_push!(trials₀, "Wfact", wfact_fun(integrator), (W, u, p, dt, t), kwargs, only)
            maybe_push!(trials₀, "ldiv!", LA.ldiv!, (X, W, u), kwargs, only)
            maybe_push!(trials₀, "T_imp!", implicit_fun(integrator), implicit_args(integrator), kwargs, only)
        end
        if !isnothing(remaining_fun(integrator))
            maybe_push!(trials₀, "T_exp_T_lim!", remaining_fun(integrator), remaining_args(integrator), kwargs, only)
        end
        maybe_push!(trials₀, "lim!", f.lim!, (Xlim, p, t, u), kwargs, only)
        maybe_push!(trials₀, "dss!", f.dss!, (u, p, t), kwargs, only)
        maybe_push!(trials₀, "cache!", f.cache!, (u, p, t), kwargs, only)
        maybe_push!(trials₀, "cache_imp!", f.cache_imp!, (u, p, t), kwargs, only)
        maybe_push!(trials₀, "step!", SciMLBase.step!, (integrator, ), kwargs, only)
#! format: on

        trials = OrderedCollections.OrderedDict()

        n_calls_per_step = get_n_calls_per_step(integrator)
        for k in keys(trials₀)
            isnothing(trials₀[k]) && continue
            trials[k] = trials₀[k] * n_calls_per_step[k]
        end

        table_summary = OrderedCollections.OrderedDict()
        for k in keys(trials)
            isnothing(trials[k]) && continue
            trial_step = haskey(trials, "step!") ? trials["step!"] : nothing
            table_summary[k] = get_summary(trials[k], trial_step)
        end

        if !isnothing(only)
            @warn "Percentages are only based on $only, pass `only = nothing` for accurately reported percentages"
        end
        tabulate_summary(table_summary; n_calls_per_step)

        return (; table_summary, trials)
    else
        @warn "`ClimaTimeSteppers.benchmark` is not yet supported for $f."
        return (; table_summary = nothing, trials = nothing)
    end
end


end
