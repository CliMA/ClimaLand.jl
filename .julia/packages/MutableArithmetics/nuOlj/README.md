# MutableArithmetics.jl

[![Stable][docs-stable-img]][docs-stable-url]
[![Dev][docs-dev-img]][docs-dev-url]
[![Build Status][build-img]][build-url]
[![Codecov][codecov-img]][codecov-url]
[![deps][deps-img]][deps-url]
[![version][version-img]][version-url]
[![pkgeval][pkgeval-img]][pkgeval-url]
[![DOI](https://proceedings.juliacon.org/papers/10.21105/jcon.00093/status.svg)](https://doi.org/10.21105/jcon.00093)


**MutableArithmetics** (MA for short) is a [Julia](http://julialang.org) package
which allows:
 * for mutable types to implement mutable arithmetics
 * for algorithms that could exploit mutable arithmetics to exploit them while
   still being completely generic.

While in some cases, similar features have been included in packages
idiosyncratically, the goal of MutableArithmetics is to provide a generic
interface to allow anyone to make use of mutability when desired.

The package allows a type to declare itself mutable through the `MA.mutability`
trait. Then the user can use the `MA.operate!!` function to write generic code
that works for arbitrary type while exploiting mutability of the type
if possible. More precisely:
 * The `MA.operate!!(op::Function, x, args...)` redirects to `op(x, args...)`
   if `x` is not mutable or if the result of the operation cannot be stored in
   `x`. Otherwise, it redirects to `MA.operate!(op, x, args...)`.
 * `MA.operate!(op::Function, x, args...)` stores the result of the operation in
   `x`. It is a `MethodError` if `x` is not mutable or if the result of the
   operation cannot be stored in `x`.

So from a generic code perspective, `MA.operate!!` can be used when the value of
`x` is not used anywhere else. This allows the code to both work for mutable and
for non-mutable type.

When the type is known to be mutable, `MA.operate!` can be used to make sure the
operation is done in-place. If it is not possible, the `MethodError` allows one
to easily fix the issue while `MA.operate!!` would have silently fallen back to
the non-mutating function.

In conclusion, the distinction between `MA.operate!!` and `MA.operate!` covers
all use case while having an universal convention accross all operations.

## Implementations

The following types and packages implement the MutableArithmetics API:

 * `Base.BigInt` in `src/interfaces/BigInt.jl`.
 * `Base.BigFloat` in `src/interfaces/BigFloat.jl`.
 * `Base.Array` in `src/interfaces/LinearAlgebra.jl`.
 * [Polynomials.jl](https://github.com/JuliaMath/Polynomials.jl) uses MA for its
   `Polynomial` type
 * [MultivariatePolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl)
   uses MA for its multivariate polynomials, as well as its two implementations
   in [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl)
   and [TypedPolynomials](https://github.com/JuliaAlgebra/TypedPolynomials.jl)
 * [JuMP](https://github.com/jump-dev/JuMP.jl) and
   [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl) use
   MA for the scalar and quadratic functions used to define an optimization
   program

In addition, the implementation of the following `Base` functionalities are
reimplemented using the MA API:
 * Matrix-matrix, matrix-vector and array-scalar multiplication including
   `SparseArrays.AbstractSparseArray`, `LinearAlgebra.Adjoint`,
   `LinearAlgebra.Transpose`, `LinearAlgebra.Symmetric`.
 * `Base.sum`, `LinearAlgebra.dot` and `LinearAlgebra.diagm`.

These methods are reimplemented in this package for several reasons:
* The implementation in `Base` does not exploit the mutability of the type
  (except for `sum(::Vector{BigInt})` which has a specialized method) and
  are hence much slower.
* Some implementations in `Base` assume the following for the types `S`, `T` used satisfy:
  - `typeof(zero(T)) == T`, `typeof(one(T)) == T`, `typeof(S + T) == promote_type(S, T)`
    or `typeof(S * T) == promote_type(S, T)` which is not true for
    instance if `T` is a polynomial variable or the decision variable of an
    optimization model.
  - The multiplication between elements of type `S` and `T` is commutative which
    is not true for matrices or non-commutative polynomial variables.

The trait defined in this package cannot make the methods for the functions
defined in Base to be dispatched to the implementations of this package.
For these to be used for a given type, it needs to inherit from `MA.AbstractMutable`.
Not that subtypes of `MA.AbstractMutable` are not necessarily mutable,
for instance, polynomial variables and the decision variable of an optimization
model are subtypes of `MA.AbstractMutable` but are not mutable.
The only purpose of this abstract type is to have `Base` methods to be dispatched
to the implementations of this package. See `src/dispatch.jl` for more details.

## Quick Example & Benchmark

```julia
using BenchmarkTools
using MutableArithmetics
const MA = MutableArithmetics

n = 200
A = rand(-10:10, n, n)
b = rand(-10:10, n)
c = rand(-10:10, n)

# MA.mul works for arbitrary types
MA.mul(A, b)

A2 = big.(A)
b2 = big.(b)
c2 = big.(c)
```

The default implementation `LinearAlgebra.generic_matvecmul!` does not exploit
the mutability of `BigInt` is quite slow and allocates a lot:
```julia
using LinearAlgebra
trial = @benchmark LinearAlgebra.mul!($c2, $A2, $b2)
display(trial)

# output

BenchmarkTools.Trial: 407 samples with 1 evaluation.
 Range (min … max):   5.268 ms … 161.929 ms  ┊ GC (min … max):  0.00% … 73.90%
 Time  (median):      5.900 ms               ┊ GC (median):     0.00%
 Time  (mean ± σ):   12.286 ms ±  21.539 ms  ┊ GC (mean ± σ):  29.47% ± 14.50%

  █▃
  ██▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅█▆▇▅▅ ▆
  5.27 ms       Histogram: log(frequency) by time      80.6 ms <

 Memory estimate: 3.66 MiB, allocs estimate: 197732.
```

In `MA.operate!(::typeof(MA.add_mul), ::Vector, ::Matrix, ::Vector)`, we
exploit the mutability of `BigInt` through the MutableArithmetics API.
This provides a significant speedup and a drastic reduction of memory usage:
```julia
trial2 = @benchmark MA.add_mul!!($c2, $A2, $b2)
display(trial2)

# output

BenchmarkTools.Trial: 4878 samples with 1 evaluation.
 Range (min … max):  908.860 μs …   1.758 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):       1.001 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):     1.021 ms ± 102.381 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

  █▅
  ██▂▂▂▇▅▇▇▅▅▅▇▅▆▄▄▅▄▄▃▄▄▃▃▂▃▃▂▃▂▂▂▂▂▂▂▂▂▂▂▁▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
  909 μs           Histogram: frequency by time         1.36 ms <

 Memory estimate: 48 bytes, allocs estimate: 3.
```

There is still 48 bytes that are allocated, where does this come from ?
`MA.operate!(::typeof(MA.add_mul), ::BigInt, ::BigInt, ::BigInt)`
allocates a temporary `BigInt` to hold the result of the multiplication.
This buffer is allocated only once for the whole matrix-vector multiplication
through the system of buffers of MutableArithmetics.
If may Matrix-Vector products need to be computed, the buffer can even be allocated
outside of the matrix-vector product as follows:
```julia
buffer = MA.buffer_for(MA.add_mul, typeof(c2), typeof(A2), typeof(b2))
trial3 = @benchmark MA.buffered_operate!!($buffer, MA.add_mul, $c2, $A2, $b2)
display(trial3)

# output

BenchmarkTools.Trial: 4910 samples with 1 evaluation.
 Range (min … max):  908.414 μs …   1.774 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     990.964 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):     1.014 ms ± 103.364 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

  █▂
  ██▃▂▂▄▄▅▆▃▄▄▅▄▄▃▃▄▃▃▃▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
  908 μs           Histogram: frequency by time         1.35 ms <

 Memory estimate: 0 bytes, allocs estimate: 0.
```
Note that there are now 0 allocations.

## Related projects

For related projects, see:

 * https://github.com/JuliaFolds/BangBang.jl
 * https://github.com/tkluck/InPlace.jl

## Citing MutableArithmetics

If you find MutableArithmetics useful in your work, we kindly request that you cite the
following [paper](https://doi.org/10.21105/jcon.00093):
```bibtex
@article{Legat2024,
  doi = {10.21105/jcon.00093},
  url = {https://doi.org/10.21105/jcon.00093},
  year = {2024},
  publisher = {The Open Journal},
  volume = {6},
  number = {58},
  pages = {93},
  author = {Benoît Legat},
  title = {MutableArithmetics: An API for mutable operations},
  journal = {Proceedings of the JuliaCon Conferences}
}
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-stable-url]: https://jump.dev/MutableArithmetics.jl/stable
[docs-dev-url]: https://jump.dev/MutableArithmetics.jl/dev

[build-img]: https://github.com/jump-dev/MutableArithmetics.jl/actions/workflows/ci.yml/badge.svg?branch=master
[build-url]: https://github.com/jump-dev/MutableArithmetics.jl/actions?query=workflow%3ACI
[codecov-img]: https://codecov.io/gh/jump-dev/MutableArithmetics.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/jump-dev/MutableArithmetics.jl

[deps-img]: https://juliahub.com/docs/MutableArithmetics/deps.svg
[deps-url]: https://juliahub.com/ui/Packages/MutableArithmetics/EoEec?t=2
[version-img]: https://juliahub.com/docs/MutableArithmetics/version.svg
[version-url]: https://juliahub.com/ui/Packages/MutableArithmetics/EoEec
[pkgeval-img]: https://juliahub.com/docs/MutableArithmetics/pkgeval.svg
[pkgeval-url]: https://juliahub.com/ui/Packages/MutableArithmetics/EoEec
