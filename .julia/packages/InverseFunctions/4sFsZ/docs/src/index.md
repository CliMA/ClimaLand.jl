# InverseFunctions.jl

```@docs
InverseFunctions
```

This package defines the function [`inverse`](@ref). `inverse(f)` returns the inverse function of a function `f`, so that `inverse(f)(f(x)) ≈ x`.

`inverse` supports mapped/broadcasted functions (via `Base.Broadcast.BroadcastFunction` or `Base.Fix1`) and function composition.

Implementations of `inverse(f)` for `identity`, `inv`, `adjoint` and `transpose` as well as for `exp`, `log`, `exp2`, `log2`, `exp10`, `log10`, `expm1`, `log1p` and `sqrt` are included.
