# DiffResults

```@meta
CurrentModule = DiffResults
```

Many differentiation techniques can calculate primal values and multiple orders of
derivatives simultaneously. In other words, there are techniques for computing `f(x)`,
`∇f(x)` and `H(f(x))` in one fell swoop!

For this purpose, DiffResults provides the `DiffResult` type, which can be passed
to in-place differentiation methods instead of an output buffer. The method
then loads all computed results into the given `DiffResult`, which the user
can then query afterwards.

Here's an example of `DiffResult` in action using ForwardDiff:

```julia
julia> using ForwardDiff, DiffResults

julia> f(x) = sum(sin, x) + prod(tan, x) * sum(sqrt, x);

julia> x = rand(4);

# construct a `DiffResult` with storage for a Hessian, gradient,
# and primal value based on the type and shape of `x`.
julia> result = DiffResults.HessianResult(x)

# Instead of passing an output buffer to `hessian!`, we pass `result`.
# Note that we re-alias to `result` - this is important! See `hessian!`
# docs for why we do this.
julia> result = ForwardDiff.hessian!(result, f, x);

# ...and now we can get all the computed data from `result`
julia> DiffResults.value(result) == f(x)
true

julia> DiffResults.gradient(result) == ForwardDiff.gradient(f, x)
true

julia> DiffResults.hessian(result) == ForwardDiff.hessian(f, x)
true
```

The rest of this document describes the API for constructing, accessing, and mutating
`DiffResult` instances. For details on how to use a `DiffResult` with a specific
package's methods, please consult that package's documentation.

## Constructing a `DiffResult`

```@docs
DiffResults.DiffResult
DiffResults.JacobianResult
DiffResults.GradientResult
DiffResults.HessianResult
```

## Accessing data from a `DiffResult`

```@docs
DiffResults.value
DiffResults.derivative
DiffResults.gradient
DiffResults.jacobian
DiffResults.hessian
```

## Mutating a `DiffResult`

```@docs
DiffResults.value!
DiffResults.derivative!
DiffResults.gradient!
DiffResults.jacobian!
DiffResults.hessian!
```
