# LazyBroadcast.jl

LazyBroadcast.jl provides a function, `lazy_broadcast`, and macro
`@lazy_broadcast`, to transform a given Julia broadcast expression into a
`Base.Broadcast.Broadcasted` object, without materializing it.

For more information about Julia broadcasting, please see
https://docs.julialang.org/en/v1/manual/arrays/#Broadcasting.

This utility is useful in a few situations:

 - Debugging broadcast machinery
 - Fusing operations in multiple broadcast expressions (alternatively, see
   [MultiBroadcastFusion.jl](https://github.com/CliMA/MultiBroadcastFusion.jl))
 - Delaying execution of a broadcast expression

For not-in-place expressions, `lazy_broadcast`/`@lazy_broadcast` simply returns
the instantiated broadcasted object, via `Base.Broadcast.instantiate
(Base.Broadcast.broadcasted(x))`, of the right-hand-side:

```julia
using Test
import LazyBroadcast: @lazy_broadcast, lazy_broadcast
import Base.Broadcast: instantiate, broadcasted, materialize

a = rand(3,3)
b = rand(3,3)

@testset "lazy_broadcast" begin
    bc = lazy_broadcast.(a .+ b) # get the broadcasted object
    @test instantiate(broadcasted(+, a, b)) == bc
    @test materialize(bc) == @. a + b # materialize the broadcasted object
end

@testset "@lazy_broadcast" begin
    bc = @lazy_broadcast @. a + b # get the broadcasted object
    @test instantiate(broadcasted(+, a, b)) == bc
    @test materialize(bc) == @. a + b # materialize the broadcasted object
end
```

`lazy_broadcast` does not support in-place expressions (as is supported in
`MultiBroadcastFusion.jl`).

## Acknowledgement

The original implementation of `LazyBroadcast` involved a similar recipe to
`MultiBroadcastFusion`, which has different yet justified needs for its
implementation. Since then, [DontMaterialize.jl](https://github.com/MasonProtter/DontMaterialize.jl) was developed, which
satisfied most needs by LazyBroadcast using a significantly more elegant approach, discussed [here](https://github.com/CliMA/LazyBroadcast.jl/issues/14). So, we've borrowed that
implementation in order to be more consistent with `Base.Broadcast` semantics,
and provide users with less surprising behavior. We've also kept `code_lowered_single_expression`, as this is one feature that `DontMaterialize.jl` does not offer (transforming broadcasted `Expr` to `Expr`).

