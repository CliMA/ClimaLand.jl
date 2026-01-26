# PositiveFactorizations

[![Build Status](https://travis-ci.org/timholy/PositiveFactorizations.jl.svg?branch=master)](https://travis-ci.org/timholy/PositiveFactorizations.jl)

## Overview

PositiveFactorizations is a package for computing a positive definite
matrix decomposition (factorization) from an arbitrary symmetric
input.  The motivating application is optimization (Newton or
quasi-Newton methods), in which the canonical search direction `-H\g`
(`H` being the Hessian and `g` the gradient) may not be a descent
direction if `H` is not positive definite.  This package provides an
efficient approach to computing `-Htilde\g`, where `Htilde` is equal
to `H` if `H` is positive definite, and otherwise is a
positive definite matrix that is "spiritually like `H`."

The approach favored here is different from the well-known
Gill-Murray-Wright approach of computing the Cholesky factorization of
`H+E`, where `E` is a minimal correction needed to make `H+E`
positive-definite (sometimes known as GMW81).  See the discussion
starting
[here](https://github.com/JuliaOpt/Optim.jl/issues/153#issuecomment-161268535)
for justification; briefly, the idea of a small correction conflates
large negative eigenvalues with numerical roundoff error, which (when
stated that way) seems like a puzzling choice.  In contrast, this
package provides methods that are largely equivalent to taking the
absolute value of the diagonals D in an LDLT factorization, and setting
any "tiny" diagonals (those consistent with roundoff error) to 1.  For
a diagonal matrix with some entries negative, this results in
approximately twice the correction used in GMW81.

## Usage

Given a symmetric matrix `H`, compute a positive factorization `F` like this:

```jl
F = cholesky(Positive, H, [pivot=Val{false}])
```

Pivoting (turned on with `Val{true}`) can make the correction smaller
and increase accuracy, but is not necessary for existence or stability.

For a little more information, call `ldlt` instead:

```jl
F, d = ldlt(Positive, H, [pivot=Val{false}])
```

`F` will be the same as for `cholesky`, but this also returns `d`, a
vector of `Int8` with values +1, 0, or -1 indicating the sign of the
diagonal as encountered during processing (so in order of rows/columns
if not using pivoting, in order of pivot if using pivoting).  This
output can be useful for determining whether the original matrix was
already positive (semi)definite.

Note that `cholesky`/`ldlt` can be used with any matrix, even
those which lack a conventional LDLT factorization.  For example, the
matrix `[0 1; 1 0]` is factored as `L = [1 0; 0 1]` (the identity matrix),
with all entries of `d` being 0.  Symmetry is assumed but not checked;
only the lower-triangle of the input is used.

`cholesky` is recommended because it is very efficient.  A slower alternative is to use `eigen`:

```jl
F = eigen(Positive, H)
```

which may be easier to reason about from the standpoint of fundamental linear algebra.
