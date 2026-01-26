# LDLFactorizations: Factorization of Symmetric Matrices

A translation of Tim Davis's Concise LDLᵀ Factorization, part of [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) with several improvements.

Please cite this repository if you use LDLFactorizations.jl in your work: see [`CITATION.bib`](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl/blob/main/CITATION.bib).

[![DOI](https://zenodo.org/badge/98073166.svg)](https://zenodo.org/badge/latestdoi/98073166)
![CI](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl/workflows/CI/badge.svg?branch=main)
[![Build Status](https://api.cirrus-ci.com/github/JuliaSmoothOptimizers/LDLFactorizations.jl.svg)](https://cirrus-ci.com/github/JuliaSmoothOptimizers/LDLFactorizations.jl)
[![codecov](https://codecov.io/gh/JuliaSmoothOptimizers/LDLFactorizations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSmoothOptimizers/LDLFactorizations.jl)
[![Documentation/stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSmoothOptimizers.github.io/LDLFactorizations.jl/stable)
[![Documentation/dev](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaSmoothOptimizers.github.io/LDLFactorizations.jl/dev)

This package is appropriate for matrices A that possess a factorization of the
form LDLᵀ without pivoting, where L is unit lower triangular and D is *diagonal* (indefinite in general), including definite and quasi-definite matrices.

LDLFactorizations.jl should not be expected to be as fast, as robust or as accurate as factorization
packages such as [HSL.jl](https://github.com/JuliaSmoothOptimizers/HSL.jl), [MUMPS.jl](https://github.com/JuliaSmoothOptimizers/MUMPS.jl) or [Pardiso.jl](https://github.com/JuliaSparse/Pardiso.jl).
Those are multifrontal and/or implement various forms of parallelism, and
employ sophisticated pivot strategies.

The main advantages of LDLFactorizations.jl are that

1. it is very short and has a small footprint;
2. it is in pure Julia, and so

   2.a. it does not require external compiled dependencies;

   2.b. it will work with multiple input data types.

Whereas MUMPS.jl, HSL.jl and Pardiso.jl only work with single and double precision
reals and complex data types, LDLFactorizations.jl accepts any numerical data type.

# Installing

```julia
julia> ]
pkg> add LDLFactorizations
```

# Usage

The only exported functions are `ldl()`, `\` and `ldiv!`.
Calling `ldl()` with a dense array converts it to a sparse matrix.
A permutation ordering can be supplied: `ldl(A, p)` where `p` is an `Int`
array representing a permutation of the integers between 1 and the order
of `A`.
If no permutation is supplied, one is automatically computed using [AMD.jl](https://github.com/JuliaSmoothOptimizers/AMD.jl).
Only the upper triangle of `A` is accessed.

`ldl` returns a factorization in the form of a `LDLFactorization` object.
The `\` and `ldiv!` methods are implemented for objects of type `LDLFactorization` so that
solving a linear system is as easy as
```julia
LDLT = ldl(A)  # LDLᵀ factorization of A

x = LDLT \ b   # solves Ax = b

ldiv!(LDLT, b)     # computes LDLT \ b in-place and overwriting b to store the result
y = similar(b)
ldiv!(y, LDLT, b)  # computes LDLT \ b in-place and store the result in y
```
The factorization can of course be reused to solve for multiple right-hand
sides.

Factors can be accessed as `LDLT.L` and `LDLT.D`, and the permutation vector as `LDLT.P`.
Because the L factor is unit lower triangular, its diagonal is not stored.
Thus the factors satisfy: PAPᵀ = (L + I) D (L + I)ᵀ.

# References

Timothy A. Davis. 2005. Algorithm 849: A concise sparse Cholesky factorization package. ACM Trans. Math. Softw. 31, 4 (December 2005), 587-591. DOI [10.1145/1114268.1114277](http://dx.doi.org/10.1145/1114268.1114277).

Like the original LDL, this package is distributed under the LGPL.

[![LGPLv3](http://www.gnu.org/graphics/lgplv3-88x31.png)](http://www.gnu.org/licenses/lgpl.html "LGPLv3")

## Bug reports and discussions

If you think you found a bug, feel free to open an [issue](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl/issues).
Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.

If you want to ask a question not suited for a bug report, feel free to start a discussion [here](https://github.com/JuliaSmoothOptimizers/Organization/discussions). This forum is for general discussion about this repository and the [JuliaSmoothOptimizers](https://github.com/JuliaSmoothOptimizers) organization, so questions about any of our packages are welcome.
