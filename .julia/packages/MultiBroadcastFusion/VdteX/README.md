# MultiBroadcastFusion.jl

A Julia package for fusing multiple broadcast expressions together.

A motivating example of this package is the following:

```julia
x1 = rand(3,3)
x2 = rand(3,3)
x3 = rand(3,3)
x4 = rand(3,3)
y1 = rand(3,3)
y2 = rand(3,3)

# 2 writes, 4 unique reads, but 8 reads including redundant ones
@. y1 = x1 * x2 + x3 * x4
@. y2 = x1 * x3 + x2 * x4
```

In this example, there are 4 unique reads, and 2 writes. However, because the reads are in two separate broadcast expressions, there are 8 reads total, including redundant ones. Another important note is that `y1` and `y2` are stored separately in memory. Fusing these operations can be achieved by changing the memory layout, and adjusting the implementation. For example:

```julia
X = map(x->Tuple(rand(4)),zeros(3,3));
Y = map(x->Tuple(rand(2)),zeros(3,3));

foo(x) = (x[1] * x[2] + x[3] * x[4], x[1] * x[3] + x[2] * x[4])
# 4 reads, 2 writes
@. Y = foo(X)
```

However, this is not an objectively better solution:

 - The memory layout, and code implementation, had to be changed in order to make this work, and this can be very difficult for a complex codebase.
 - Memory acces of `X` and `Y` is now _strided_, which could result in less performant code than a single fused loop with more contiguous memory.

Ideally, we would like for the loops to be fused with the more contiguous data layouts:

```julia
x1 = rand(3,3)
x2 = rand(3,3)
x3 = rand(3,3)
x4 = rand(3,3)
y1 = rand(3,3)
y2 = rand(3,3)

# 2 writes, 4 unique reads. The compiler can hoist the redundant memory reads here.
for i in eachindex(x1,x2,x3,x4,y1,y2)
  y1[i] = x1[i] * x2[i] + x3[i] * x4[i]
  y2[i] = x1[i] * x3[i] + x2[i] * x4[i]
end
```

With this package, we can apply `@fused_direct` to reduce the number of reads and preserve the memory layout:

```julia
import MultiBroadcastFusion as MBF
x1 = rand(3,3)
x2 = rand(3,3)
x3 = rand(3,3)
x4 = rand(3,3)
y1 = rand(3,3)
y2 = rand(3,3)

# 4 reads, 2 writes
MBF.@fused_direct begin
  @. y1 = x1 * x2 + x3 * x4
  @. y2 = x1 * x3 + x2 * x4
end
```

This is achieved by fusing the loops and inlining with the given data, resulting in the compiler being able to perform Common-SubExpression Elimination (CSE) on the memory loads.

## Custom implementations

Users can write custom implementations, using the `@make_type` and `@make_fused` macros, and then defining `Base.copyto!` on the type you've defined

```julia
import MultiBroadcastFusion as MBF
import MultiBroadcastFusion: fused_direct

MBF.@make_type MyFusedMultiBroadcast
MBF.@make_fused fused_direct MyFusedMultiBroadcast my_fused
# Now, `@fused_direct` will call `Base.copyto!(::MyFusedMultiBroadcast)`. Let's define it:
function Base.copyto!(fmb::MyFusedMultiBroadcast)
    pairs = fmb.pairs
    destinations = map(x->x.first, pairs)
    @inbounds for i in eachindex(destinations)
        # does `@inline pair.first[i] = pair.second[i]` for all pairs
        MBF.rcopyto_at!(pairs, i)
    end
    return nothing
end

x1 = rand(3,3)
x2 = rand(3,3)
x3 = rand(3,3)
x4 = rand(3,3)
y1 = rand(3,3)
y2 = rand(3,3)

# 4 reads, 2 writes
@my_fused begin
  @. y1 = x1 * x2 + x3 * x4
  @. y2 = x1 * x3 + x2 * x4
end
```

## Writing custom macros

Users can also write custom macros with, for example,

```julia
import MultiBroadcastFusion as MBF

struct FusedMultiBroadcast{T}
    pairs::T
end
macro get_fused_multi_broadcast(expr)
    _pairs = gensym()
    quote
        $_pairs = $(esc(MBF.fused_direct(expr)))
        FusedMultiBroadcast($_pairs)
    end
end
```

This can be helpful for inspecting multibroadcast objects.
