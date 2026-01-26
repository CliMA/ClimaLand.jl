# PtrArrays

[![Build Status](https://github.com/LilithHafner/PtrArrays.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/LilithHafner/PtrArrays.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/LilithHafner/PtrArrays.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/LilithHafner/PtrArrays.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/P/PtrArrays.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/P/PtrArrays.html)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![deps](https://juliahub.com/docs/General/PtrArrays/stable/deps.svg)](https://juliahub.com/ui/Packages/General/PtrArrays?t=2)

Do you miss playing hide and seek with memory leaks? Do you find GC overhead problematic?
PtrArrays.jl can take you back to the good old days of manual memory management.
See also [Bumper.jl](https://github.com/MasonProtter/Bumper.jl) if you want to avoid GC
overhead and don't like hide and seek.

This package provides `malloc(T, dims...)` which allocates an `AbstractArray{T}` with the
provided `dims`. If you want, you can call `free` on the array once you're done using it
but it can be more fun to see what happens if you don't.

Example usage

```julia
julia> malloc(Int, 4)
4-element PtrArray{Int64, 1}:
 1053122630
          0
  936098496
  936099008

julia> free(ans)

julia> malloc(Int, 4, 4)
4×4 PtrArray{Int64, 2}:
       923300075       1046634192       1046634192       1046634408
               0              120              124              152
               0                0                0                0
 281474587621896  281474587621899  281474587621900  281474587621896

julia> free(ans)
```

Benchmarks:

```julia
using PtrArrays
function f(n)
    x = malloc(Int, n)
    try
        sum(x) # Let's see what we get!
    finally
        free(x) # Putting the `free` call in a finally block makes memory leaks less common
    end
end

f(1000)
# 6474266410623015

using Chairmarks
@b f(1000)
# 101.317 ns

function g(n)
    x = Vector{Int}(undef, n)
    sum(x) # Let's see what we get!
end

@b g(1000)
# 130.125 ns (3 allocs: 7.875 KiB)
```

The whole package's source code is only about 44 lines (excluding comments and whitespace),
half of which is re-implementing Julia's buggy `Core.checked_dims` function.
[Read it here](https://github.com/LilithHafner/PtrArrays.jl/blob/main/src/PtrArrays.jl)

## Alternatives

[Bumper.jl](https://github.com/MasonProtter/Bumper.jl) provides bump allocators which allow
you to manage your own allocation stack, bypassing the Julia GC.

[StaticTools.jl](https://github.com/brenhinkeller/) is a much larger package which provides,
among other things, a
[MallocArray](https://brenhinkeller.github.io/StaticTools.jl/dev/#StaticTools.MallocArray)
type that behaves similarly to PtrArray. As StaticTools.jl is meant to run without the Julia
runtime, it does not make use of features such as exceptions. Consequently, all usages of
`StaticTools.MallocArray` are implicitly annotated with `@inbounds` and out of bounds
accesses are UB instead of errors, StaticTools.MallocArray with invalid indices will return
a null pointer with size zero while PtrArrays will throw, etc. In general, StaticTools.jl's
`MallocArray` can be thought of as "unsafe" while `PtrArray` is "safe" (though memory leaks
will still occur if you call `malloc` and fail to call `free`)

[MallocArrays.jl](https://github.com/LilithHafner/PtrArrays.jl/tree/0b6dbdc012e1058b2b64d0f94863eff4120def85)
is the original name of this package. It is obsolete.
