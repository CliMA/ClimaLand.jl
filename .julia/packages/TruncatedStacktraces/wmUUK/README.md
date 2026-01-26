# TruncatedStacktraces.jl: Truncated and Simpler Stacktraces for the Julia Programming Language

Don't you wish Julia stacktraces were simpler? Introducing TruncatedStacktraces.jl! The purpose of this
package is to give package authors a single uniform system for implementing truncation of type printing
in stack traces.

> **Note**
> Starting v1.10 a similar feature is inbuilt into julia. Starting julia v1.10, this package does nothing!

## Enabling TruncatedStacktraces.jl

TruncatedStacktraces.jl is currently disabled by default, as it causes invalidations which will slow down package loading.

It can be enabled using Preferences.jl. To enable it, create a `LocalPreferences.toml` with the following entry:

```toml
[TruncatedStacktraces]
disable = false
```

Alternatively, you can generate the `LocalPreferences.toml` using:

```julia
using Preferences, UUIDs

using TruncatedStacktraces
Preferences.set_preferences!(TruncatedStacktraces, "disable" => false)

# OR if you don't want to load TruncatedStacktraces.jl

Preferences.set_preferences!(UUID("781d530d-4396-4725-bb49-402e4bee1e77"), "disable" => false)
```

In either case, you need to reload your packages (depending on TruncatedStacktraces) for the
change to take effect.

**TruncatedStacktraces is known to create invalidations, to remove these simply set the preference to disable it!**

## Users: How to Interact with TruncatedStacktraces.jl

If a package you are using is making use of TruncatedStacktraces.jl, you will see shorter stack traces. Everything
is easier to read by default! This looks like:

```julia
 [14] initialize!(integrator::ODEIntegrator{true, Tsit5{Static.False, …}, Vector{Float64}, Float64, …}, cache::Tsit5Cache{Vector{Float64}, …})
    @ OrdinaryDiffEq C:\Users\accou\.julia\packages\OrdinaryDiffEq\0Pm1I\src\perform_step\low_order_rk_perform_step.jl:766
```

But if you want to see the type in full glory, say to share with developers on Discourse, then you can opt to show
the entire stacktrace via simply running:

```julia
TruncatedStacktraces.VERBOSE[] = true
```

then if you run the code to error again, it will print out exactly what everyone wants to read:

```julia
 [14] initialize!(integrator::OrdinaryDiffEq.ODEIntegrator{Tsit5{typeof(OrdinaryDiffEq.trivial_limiter!), typeof(OrdinaryDiffEq.trivial_limiter!), Static.False}, true, Vector{Float64}, Nothing, Float64, SciMLBase.NullParameters, Float64, Float64, Float64, Float64, Vector{Vector{Float64}}, ODESolution{Float64, 2, Vector{Vector{Float64}}, Nothing, Nothing, Vector{Float64}, Vector{Vector{Vector{Float64}}}, ODEProblem{Vector{Float64}, Tuple{Float64, Float64}, true, SciMLBase.NullParameters, ODEFunction{true, SciMLBase.AutoSpecialize, FunctionWrappersWrappers.FunctionWrappersWrapper{Tuple{FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{Float64}, Vector{Float64}, SciMLBase.NullParameters, Float64}}, FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, SciMLBase.NullParameters, Float64}}, FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, Vector{Float64}, SciMLBase.NullParameters, ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}}, FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, SciMLBase.NullParameters, ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}}}, false}, LinearAlgebra.UniformScaling{Bool}, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, typeof(SciMLBase.DEFAULT_OBSERVED), Nothing, Nothing}, Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}}, SciMLBase.StandardODEProblem}, Tsit5{typeof(OrdinaryDiffEq.trivial_limiter!), typeof(OrdinaryDiffEq.trivial_limiter!), Static.False}, OrdinaryDiffEq.InterpolationData{ODEFunction{true, SciMLBase.AutoSpecialize, FunctionWrappersWrappers.FunctionWrappersWrapper{Tuple{FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{Float64}, Vector{Float64}, SciMLBase.NullParameters, Float64}}, FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, SciMLBase.NullParameters, Float64}}, FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, Vector{Float64}, SciMLBase.NullParameters, ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}}, FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, SciMLBase.NullParameters, ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}}}, false}, LinearAlgebra.UniformScaling{Bool}, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, typeof(SciMLBase.DEFAULT_OBSERVED), Nothing, Nothing}, Vector{Vector{Float64}}, Vector{Float64}, Vector{Vector{Vector{Float64}}}, OrdinaryDiffEq.Tsit5Cache{Vector{Float64}, Vector{Float64}, Vector{Float64}, typeof(OrdinaryDiffEq.trivial_limiter!), typeof(OrdinaryDiffEq.trivial_limiter!), Static.False}}, DiffEqBase.DEStats, Nothing}, ODEFunction{true, SciMLBase.AutoSpecialize, FunctionWrappersWrappers.FunctionWrappersWrapper{Tuple{FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{Float64}, Vector{Float64}, SciMLBase.NullParameters, Float64}}, FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, SciMLBase.NullParameters, Float64}}, FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, Vector{Float64}, SciMLBase.NullParameters, ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}}, FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}, SciMLBase.NullParameters, ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}}}}, false}, LinearAlgebra.UniformScaling{Bool}, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, typeof(SciMLBase.DEFAULT_OBSERVED), Nothing, Nothing}, OrdinaryDiffEq.Tsit5Cache{Vector{Float64}, Vector{Float64}, Vector{Float64}, typeof(OrdinaryDiffEq.trivial_limiter!), typeof(OrdinaryDiffEq.trivial_limiter!), Static.False}, OrdinaryDiffEq.DEOptions{Float64, Float64, Float64, Float64, PIController{Rational{Int64}}, typeof(DiffEqBase.ODE_DEFAULT_NORM), typeof(LinearAlgebra.opnorm), Nothing, CallbackSet{Tuple{}, Tuple{}}, typeof(DiffEqBase.ODE_DEFAULT_ISOUTOFDOMAIN), typeof(DiffEqBase.ODE_DEFAULT_PROG_MESSAGE), typeof(DiffEqBase.ODE_DEFAULT_UNSTABLE_CHECK), DataStructures.BinaryHeap{Float64, DataStructures.FasterForward}, DataStructures.BinaryHeap{Float64, DataStructures.FasterForward}, Nothing, Nothing, Int64, Tuple{}, Tuple{}, Tuple{}}, Vector{Float64}, Float64, Nothing, OrdinaryDiffEq.DefaultInit}, cache::OrdinaryDiffEq.Tsit5Cache{Vector{Float64}, Vector{Float64}, Vector{Float64}, typeof(OrdinaryDiffEq.trivial_limiter!), typeof(OrdinaryDiffEq.trivial_limiter!), Static.False})
    @ OrdinaryDiffEq C:\Users\accou\.julia\packages\OrdinaryDiffEq\0Pm1I\src\perform_step\low_order_rk_perform_step.jl:766
 ```
 
 Beautiful. You can turn it back into the not beautiful short stacktrace with the command:

 ```julia
TruncatedStacktraces.VERBOSE[] = false
```

## How to Opt A Package Into TruncatedStacktraces.jl

Opting into TruncatedStacktraces.jl is easy: for every type that you want to omit the printing of something,
use the macro `TruncatedStacktraces.@truncate_stacktrace` like:

```julia
TruncatedStacktraces.@truncate_stacktrace ODEProblem 3 1 2
```

where `3 1 2` gives the order of the types to print, with indices corresponding to the original type. For example,
on a type `MyType{T1,T2,T3,T4}`, this will change the stacktrace printing to default to `MyType{T3,T1,T2,…}`.

For any new error exception you add to your package, make sure to include the note from TruncatedStacktraces.jl on
how to effect the type printing. This is done by adding `println(io, VERBOSE_MSG)` to the bottom of any error message.

## Default values

* `TruncatedStacktraces.VERBOSE[]` defaults to `false` for non-CI workflows and to `true` for CI jobs.
* `TruncatedStacktraces.DISABLE` defaults to `true`.

## How It's Implemented

This is done by writing an overload on `Base.show` on the DataType which is conditional on `TruncatedStacktraces.VERBOSE[]`.
For example, the following does this for the `SciMLBase.ODEProblem`:

```julia
@static if !TruncatedStacktraces.DISABLE
function Base.show(io::IO,
                   t::Type{<:ODEProblem{uType, tType, isinplace}}) where {uType, tType, isinplace}
    if TruncatedStacktraces.VERBOSE[]
        invoke(show, Tuple{IO, Type}, io, t)
    else
        print(io, "ODEProblem{$isinplace,$uType,$tType,…}")
    end
end
end
```

## FAQ: Why is this not in Base Julia?

There are attempts like https://github.com/JuliaLang/julia/pull/48444, but no one agrees on what exactly to do and how to
make it perfect. So until people agree, we can use this solution as a nice hack that gets the job done 90%.

## Related Projects

Check out https://github.com/BioTurboNick/AbbreviatedStackTraces.jl which doesn't change type printing but instead the
number of calls which are shown.
