# AMD

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3783808.svg)](https://doi.org/10.5281/zenodo.3381898)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSmoothOptimizers.github.io/AMD.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSmoothOptimizers.github.io/AMD.jl/dev)
[![CI](https://github.com/JuliaSmoothOptimizers/AMD.jl/workflows/CI/badge.svg?branch=main)](https://github.com/JuliaSmoothOptimizers/AMD.jl/actions)
[![Build Status](https://api.cirrus-ci.com/github/JuliaSmoothOptimizers/AMD.jl.svg)](https://cirrus-ci.com/github/JuliaSmoothOptimizers/AMD.jl)
[![codecov.io](https://codecov.io/github/JuliaSmoothOptimizers/AMD.jl/coverage.svg?branch=main)](https://codecov.io/github/JuliaSmoothOptimizers/AMD.jl?branch=main)

Given a square sparse matrix, compute an approximate minimum degree ordering.
This package is an interface to the AMD library (Amestoy, Davis and Duff) and COLAMD library (Liromore, Davis, Gilberg an Ng).

### Installing

```julia
julia> ]
pkg> add AMD
pkg> test AMD
```

### Algorithms

**amd**: an approximate minimum degree ordering algorithm for the sparse factorization of square matrices.

**symamd**: an approximate minimum degree ordering algorithm for the sparse factorization of symmetric matrices.

**colamd**: an approximate minimum degree column ordering algorithm for the sparse factorization of arbitrary, square or rectangular, matrices.

### Examples

In the simplest case:

```julia
using AMD
A = sprand(10, 10, .5)
p_amd = amd(A)
p_symamd = symamd(A)
p_colamd = colamd(A)
```

If statistics on the permutation are of interest and/or for changing the
default control parameters:

```julia
julia> meta = Amd{Clong}();  # because A's index type is Int64 on my platform
julia> # optionally change meta.control: ?Amd
julia> p = amd(A, meta)
julia> print(meta)
Control:
  dense row parameter: 10.0
  aggressive absorption: 1.0
Info:
  status: ok
  matrix size: 10.0
  number of nonzeros: 54.0
  pattern symmetry: 0.5
  number of nonzeros on diagonal: 6.0
  number of nonzeros in A + A': 72.0
  number of dense columns: 0.0
  memory used: 1408.0
  number of garbage collections: 0.0
  approx number of nonzers in factor: 38.0
  number of float divides: 38.0
  number of float * or - for LDL: 114.0
  number of float * or - for LU: 190.0
  max nonzeros in any column of factor: 8.0
```

The **amd** algorithm computes a fill-reducing permutation based on the sparsity pattern of A + Aᵀ. The input pattern can be anything: diagonal entries will be ignored and the rest will be used to implicitly work on the pattern of A + Aᵀ. Thus if A is symmetric, it is sufficient to supply the strict lower or upper triangle only.

### References

1. P. R. Amestoy, T. A. Davis and I. S. Duff. An Approximate Minimum Degree
   Ordering Algorithm. *SIAM Journal on Matrix Analysis and Applications*, 17(4),
   pp. 886&ndash;905, 1996.
   DOI [10.1137/S0895479894278952](https://doi.org/10.1137/S0895479894278952)
2. P. R. Amestoy, T. A. Davis, and I. S. Duff. Algorithm 837: An approximate
   minimum degree ordering algorithm. *ACM Transactions on Mathematical
   Software*, 30(3), pp. 381&ndash;388, 2004.
   DOI [10.1145/1024074.1024081](https://doi.org/10.1145/1024074.1024081)
3. T. A. Davis, J. R. Gilbert, S. Larimore, E. Ng. Algorithm 836: COLAMD,
  an approximate column minimum degree ordering algorithm, *ACM
  Transactions on Mathematical Software*, 30(3), pp. 377&ndash;380, 2004.
  DOI [10.1145/1024074.1024080](https://doi.org/10.1145/1024074.1024080)
