# AxisAlgorithms

[![CI](https://github.com/timholy/AxisAlgorithms.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/timholy/AxisAlgorithms.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/timholy/AxisAlgorithms.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/timholy/AxisAlgorithms.jl)

AxisAlgorithms is a collection of filtering and linear algebra algorithms for multidimensional arrays.
For algorithms that would typically apply along the columns of a matrix, you can instead pick an arbitrary axis (dimension).

Note that all functions come in two variants, a `!` version that uses pre-allocated output (where the output is
the first argument) and a version that allocates the output. Below, the `!` versions will be described.

### Tridiagonal and Woodbury inversion

If `F` is an LU-factorization of a tridiagonal matrix, or a [Woodbury matrix](WoodburyMatrices.jl) created from such a factorization,
then `A_ldiv_B_md!(dest, F, src, axis)` will solve the equation `F\b` for 1-dimensional slices
along dimension `axis`.
Unlike many linear algebra algorithms, this one is safe to use as a mutating algorithm with `dest=src`.
The tridiagonal case does not create temporaries, and it has excellent cache behavior.

### Matrix multiplication

Multiply a matrix `M` to all 1-dimensional slices along a particular dimension.
Here you have two algorithms to choose from:

- `A_mul_B_perm!(dest, M, src, axis)` uses `permutedims` and standard BLAS-accelerated routines; it allocates temporary storage.
- `A_mul_B_md!(dest, M, src, axis)` is a non-allocating naive routine. This also has optimized implementations for sparse `M` and 2x2 matrices.

In general it is very difficult to get efficient cache behavior for multidimensional multiplication, and often using `A_mul_B_perm!` is the best strategy.
However, there are cases where `A_mul_B_md!` is faster.
It's a good idea to time both and see which works better for your case.
