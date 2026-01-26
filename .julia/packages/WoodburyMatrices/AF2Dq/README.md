# WoodburyMatrices

[![Build Status](https://github.com/JuliaLinearAlgebra/WoodburyMatrices.jl/workflows/CI/badge.svg)](https://github.com/JuliaLinearAlgebra/WoodburyMatrices.jl/actions?query=workflow%3ACI+branch%3Amaster)
[![Coverage Status](http://codecov.io/github/JuliaLinearAlgebra/WoodburyMatrices.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaLinearAlgebra/WoodburyMatrices.jl?branch=master)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

This package provides support for the [Woodbury matrix identity](http://en.wikipedia.org/wiki/Woodbury_matrix_identity) for the Julia programming language.  This is a generalization of the [Sherman-Morrison formula](http://en.wikipedia.org/wiki/Sherman%E2%80%93Morrison_formula). Note that the Woodbury matrix identity is notorious for floating-point roundoff errors, so be prepared for a certain amount of inaccuracy in the result.

## Usage

### Woodbury Matrices
```julia
using WoodburyMatrices
W = Woodbury(A, U, C, V)
```
creates a `Woodbury` matrix from the `A`, `U`, `C`, and `V` matrices representing `A + U*C*V`. These matrices can be dense or sparse (or generally any type of `AbstractMatrix`), with the caveat that
`inv(inv(C) + V*(A\U))` will be calculated explicitly and hence needs to be representable with the available resources.
(In many applications, this is a fairly small matrix.)

Here are some of the things you can do with a Woodbury matrix:

- `Matrix(W)` converts to its dense representation.
- `W\b` solves the equation `W*x = b` for `x`.
- `W*x` computes the product.
- `det(W)` computes the determinant of `W`.
- `α*W` and `W1 + W2`.

It's worth emphasizing that `A` can be supplied as a factorization, which makes `W\b` and `det(W)` more efficient.

### SymWoodbury Matrices

```julia
using WoodburyMatrices
W = SymWoodbury(A, B, D)
```
creates a `SymWoodbury` matrix, a symmetric version of a Woodbury matrix representing `A + B*D*B'`. In addition to the features above, `SymWoodbury` also supports "squaring" `W*W`.

## Efficiency and thread-safety

If passed the keyword argument `allocatetmp=true`, (Sym)Woodbury allocates internal temporary storage for intermediate results in `W*v` and `W\v` where `v` is a vector.
This eliminates memory allocation for these common operations.

However, using the same `W` across multiple threads can lead to race conditions. Hence, this optimization is opt-in and should only be used if you know it is safe.
