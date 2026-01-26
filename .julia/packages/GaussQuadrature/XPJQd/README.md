[![Build Status](https://travis-ci.org/billmclean/GaussQuadrature.jl.svg?branch=master)](https://travis-ci.org/billmclean/GaussQuadrature.jl)

GaussQuadrature.jl
==================

Julia package to compute points and weights for Gauss quadrature rules
using the Golub-Welsch algorithm.

Handles the classical Legendre, Chebyshev, Jacobi, Laguerre and Hermite 
weights, as well as a logarithmic weight function.  In fact, the 
Gauss rule is available for any custom weight function such that the 
coefficients are known for the three-term recurrence relation satisfied 
by the associated orthogonal polynomials.  The modified Chebyshev 
algorithm is provided to determine these coefficients from the modified 
moments of the weight function.

The Lobatto and Radau variants of all these rules are also provided by
appropriate choice of the `endpt` argument: `neither` (the default), 
`both`, `left` or `right`.

For example, to obtain a plain 5-point Gauss-Legendre rule with weight
function `w(x)=1` on the interval `-1 < x < -1` do

    julia> using GaussQuadrature
    julia> x, w = legendre(5)

whereas for the Lobatto version do

    julia> x, w = legendre(5, both)

The package supports `BigFloat`s; for example,

    julia> x, w = legendre(BigFloat, 5)

gives a plain 5-point Gauss-Legendre rule accurate to about 75 significant
figures.

Read the initial comments in the src/GaussQuadrature.jl module
for full details, or read the help documentation for the individual
functions called `legendre`, `chebyshev`, `jacobi`, `laguerre`, `hermite`, 
`logweight` and `custom_gauss_rule`.
