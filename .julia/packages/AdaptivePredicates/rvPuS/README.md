# AdaptivePredicates

A Julia port of "Routines for Arbitrary Precision Floating-point Arithmetic and Fast Robust Geometric Predicates"
by Jonathan Richard Shewchuk. https://www.cs.cmu.edu/~quake/robust.html

The package provides four predicates. In the functions below, all the inputs should be `NTuple`s with either `Float64` or `Float32` coordinates; complex inputs can be used for `orient2` and `incircle`.

- `orient2(pa, pb, pc)`: Given three points `pa`, `pb`, and `pc` in two dimensions, returns a positive value if the points are in counter-clockwise order; a negative value if they occur in clockwise order; and zero if they are collinear. Equivalently, returns a positive value if `pc` is left of the oriented line from `pa` to `pb`; a negative value if `pc` is right of this line; and zero if they are collinear.
- `orient3(pa, pb, pc, pd)`: Given four points `pa`, `pb`, `pc`, and `pd` in three dimensions, define the oriented plane on which the triangle `(pa, pb, pc)` is positively oriented. Returns a positive value if `pd` is below this plane; a negative value if `pd` is above this plane; and zero if the points are coplanar.
- `incircle(pa, pb, pc, pd)`: Given four points `pa`, `pb`, `pc`, and `pd` in two dimensions, returns a positive value if `pd` is inside the circle through `pa`, `pb`, and `pc`; a negative value if `pd` is outside this circle; and zero if `pd` is on the circle.
- `insphere(pa, pb, pc, pd, pe)`: Given five points `pa`, `pb`, `pc`, `pd`, and `pe` in three dimensions, returns a positive value if `pe` inside of the sphere through `pa`, `pb`, `pc`, and `pd`; a negative value if `pe` is outside this sphere; and zero if `pe` is on the sphere.

We also define the functions `orient2p`, `orient3p`, `incirclep`, and `inspherep` which simply return the sign of the corresponding predicate. For example,
```julia-repl
julia> using AdaptivePredicates

julia> pa, pb, pc = (0.2, 0.3), (0.1, -0.5), (0.7, 0.3);

julia> orient2(pa, pb, pc)
0.39999999999999997

julia> orient2p(pa, pb, pc)
1

julia> pa, pb, pc, pd, pe = (0.3f0, 0.3f0, 0.17f0), (-0.3f0, 1.71f0, 0.0f0), (0.0f0, 0.0f0, 5.0f0), (1.1f0, -0.53f0, 1.2f0), (0.5f0, 0.50f0, 0.5f0);

julia> insphere(pa, pb, pc, pd, pe)
-5.021922f0

julia> inspherep(pa, pb, pc, pd, pe)
-1
```

## Installation

If you want to use the package, you can do
```julia
using Pkg
Pkg.add("AdaptivePredicates")
using AdaptivePredicates
```

## Other Functions 

All the functions from the `predicates.c` file from Shewchuk's original code have been included in this package. This includes,

- All macros have been implemented as functions, e.g. `Fast_Two_Sum` and `Four_Four_Sum`.
- All arithmetic functions have been implemented, e.g. `grow_expansion` and `scale_expansion_zeroelim`.
- All the predicates have been implemented. In particular, not only have `orient2`, `orient3`, `incircle`, and `insphere` been implemented, but also the forms with the suffixes `fast`, `exact`, and `slow` (and `adapt`, but this is what `orient2`, `orient3`, `incircle`, and `insphere` use anyway).

Only the functions `orient2`, `orient3`, `incircle`, and `insphere` have been marked as `public`, as well as their `p` and `fast` counterparts.

## Caveats

Shewchuk's original paper gives no analysis in the presence of underflow or overflow. The only mention of it is:

> This article does not address issues of overflow and underflow, so I allow the exponent to be an integer in the range  $[-\infty, \infty]$. (Fortunately, many applications have inputs whose exponents fall within a circumscribed range. The four predicates implemented for this article will not overflow nor underflow if their inputs have exponents in the range $[-142, 201]$ and IEEE 754 double precision arithmetic is used.)

- Richard Shewchuk, J. Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates. Discrete Comput Geom 18(3), 305–363 (1997)

Note that this range comes from the `insphere` predicate. The number range is much wider for `orient2`, for example, since it requires far fewer additions, subtractions, and multiplications.

Thus, for some numbers, the values returned from these predicates may be invalid due to underflow or overflow. In ranges where this is a concern, you should use [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl) instead. If you need the values of the predicates and not just their signs, but are outside of the range valid for AdaptivePredicates.jl, you are unfortunately out of luck.

If you need more information about how these predicates work, you should refer to Shewchuk's paper.

## License

The original code is in the public domain and this Julia port is under the MIT License.
