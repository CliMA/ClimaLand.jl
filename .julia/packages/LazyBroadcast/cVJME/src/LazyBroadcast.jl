module LazyBroadcast

using Base.Broadcast: materialize, instantiate

export lazy_broadcast, materialize

"""
    function lazy_broadcast end # exported
    const lazy = lazy_broadcast # not exported

This function has no methods and is only meant to be used to consume a broadcast
expression, causing it to not materialize, allowing it to be used in a lazy
manner and be consumed later.

For example, consider the situation where one wants to break up a complicated
broadcast expression into multiple steps, and then sum up all of the
components:

```julia
julia> function foo(x)
           y = x .+ x
           z = 2 .* y
           sum(z)
       end;

julia> @benchmark foo(v) setup=(v=rand(10))
BenchmarkTools.Trial: 10000 samples with 995 evaluations.
 Range (min … max):  31.405 ns …  4.801 μs  ┊ GC (min … max):  0.00% … 98.56%
 Time  (median):     34.809 ns              ┊ GC (median):     0.00%
 Time  (mean ± σ):   45.504 ns ± 93.354 ns  ┊ GC (mean ± σ):  20.30% ± 11.70%

  █▅▃▂                                                        ▁
  █████▆▅▅▅▅▅▅▁▄▅▅▇▅▃▁▁▁▄▅▅▁▁▁▁▁▃▁▁▁▁▁▄▁▅▁▁▁▁▁▃▅▆▅▄▅▄▁▄▄▆▆▅▄▅ █
  31.4 ns      Histogram: log(frequency) by time       298 ns <

 Memory estimate: 288 bytes, allocs estimate: 4.
 ```

 This is significantly slower than it needs to be because new arrays need to be
 allocated for `y` and `z`, and the data needs to be passed over multiple times
 because the broadcast kernels are not 'fused'.

`DontMaterialize` gives a simple way to avoid these allocations and retain
broadcast fusion:
```julia
julia> function bar(x)
           y = lazy_broadcast.(x .+ x)
           z = lazy_broadcast.(2 .* y)
           sum(z)
       end;

julia> @benchmark bar(v) setup=(v=rand(10))
BenchmarkTools.Trial: 10000 samples with 1000 evaluations.
 Range (min … max):  5.931 ns … 59.562 ns  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     6.252 ns              ┊ GC (median):    0.00%
 Time  (mean ± σ):   6.435 ns ±  2.767 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

        ▁█
  ▃▂▃▃▄▇██▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▁▂▂▂▂▂▂▂▂▂▂ ▂
  5.93 ns        Histogram: frequency by time        8.75 ns <

 Memory estimate: 0 bytes, allocs estimate: 0.
```

the result of a `lazy_broadcast` call can be collected into an array with the
`materialize` function (re-exported here from `Base.Broadcast`):

```julia
julia> lazy_broadcast.(2 .* [1,2,3])
Broadcasted{Base.Broadcast.DefaultArrayStyle{1}}(*, (2, [1, 2, 3]))

julia> materialize(ans)
3-element Vector{Int64}:
 2
 4
 6
```
"""
function lazy_broadcast end
const lazy = lazy_broadcast

struct LazyBroadcasted{T}
    value::T
end
Base.Broadcast.broadcasted(::typeof(lazy_broadcast), x) = LazyBroadcasted(x)
Base.materialize(x::LazyBroadcasted) = instantiate(x.value)

import Base.Broadcast: instantiate, broadcasted
macro lazy_broadcast(expr)
    return quote
        LazyBroadcast.lazy_broadcast.($(esc(expr)))
    end
end

macro lazy(expr)
    return quote
        LazyBroadcast.lazy_broadcast.($(esc(expr)))
    end
end

include("code_lowered_single_expression.jl")

end # module LazyBroadcast
