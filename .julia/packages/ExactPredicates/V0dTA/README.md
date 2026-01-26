# ExactPredicates.jl

→ [Documentation](https://lairez.github.io/ExactPredicates.jl/)

This package provides fast and robust predicates for Euclidean geometry. It is
implemented as a code generator that transform any homogeneous polynomial into
an algorithm that rigorously determines its sign on any given input points.

### Installation

```julia
]add ExactPredicates
```

## Features

* Planar predicates: `orient` (left or right of a line), `incircle` (inside or
  outside a circle), `closest` (closer to this or that point), `meet` (intersection of line segments), etc.
* Spatial predicates: `orient` (above or below a plane), `insphere` (inside or
  outside a ball), `closest` (closer to this or that point), etc.
* Simplistic API, points are just `Tuple{Float64,Float64}` or `Tuple{Float64,Float64,Float64}`
* Accepts anything convertible to `Tuple`
* Extensible
* MIT license


### Robust

Robust means that the code:

- raises an exception on `NaN` and `Inf` arguments;
- gives a correct answer on all other inputs with `Float64` coordinates, no matter what (overflow, underflow, etc.);
- in particular, no restriction on the coordinate range.


#### Why robustness matters?

When the input data is approximate it looks absurd to insist on exact
computation. To determine if a point is inside or outside a circle, one may need
more precision than the data really convey. So what is the point?

Robust computation is important because it guarantees *soundness* with respect
to some combinatorial properties of the predicates, on which many algorithms
rely. For example
> `orient(a, b, c) == orient(b, c, a) == orient(c, a, b)`,

this
is a very basic geometric observation, but a floating computation may fail to
see this.


> “Inexact versions of these tests *[orient and incircle]* are vulnerable to roundoff error, and the wrong
> answers they produce can cause geometric algorithms to hang, crash, or produce
> incorrect output.”
>
> —Jonathan Shewchuk, *Robust Adaptive Floating-point Geometric Predicates*

### Fast

The reference point is the cost of the pure floating point evaluation, without any certification.
To evaluate the performance, one need to distinguish several scenarios.

#### The quick path

For normal input (for example random input), the floating point computation is
correct. The only overhead is the computation of the error bound. Expect a 2×
slowdown with respect to the reference. On complicated predicates, like
`insphere`, the floating point computation dominates and the slowdown is only
1.3×.

The emphasis in the package is to make the quick path as fast as possible.

#### The slow path

If the error estimation fails to certify the floating point computation, the
code falls back to interval arithmetic, using [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl). It will work especially well if the input points have small integer coordinates.
Expect a 50× slowdown.

In this scenario, adaptive methods, like the one famously implemented by [Shewchuk](https://www.cs.cmu.edu/~quake/robust.html), may be more performant.

#### The worst path

When the input points are maximally degenerate, the code will resort to exact arithmetic using `Rational{BigInt}`.
The computation will always succeed, but expect a 50-100× slowdown.


## Usage

### Basic usage 

The type for representing points is `NTuple{N, Float64}`, where `N` is 2 or 3. Very concretely, that is `Tuple{Float64,Float64}` or `Tuple{Float64,Float64,Float64}`.

```julia
p, q, r, a = (1.0, 3.0), (1.5, 10.0), (-87.0, 1e64), (1e-100, 3.0)

incircle(p, q, r, a)
# -> 1
```

### Working with other types of points

All the predicates will work with a type `T`, if one of the function `Tuple(::T)` or
`coord(::T)` is defined and outputs a `NTuple{N, Float64}` that contains the coordinates the
coordinates. Naturally, the computation is only robust if the conversion is robust too.
Here is a typical use.

```julia
using ExactPredicates
struct Point
    x :: Float64
    y :: Float64
end

Tuple(p :: Point) = (p.x, p.y)
incircle(Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0), Point(.5, .5))
```

A nice type to represent points in the plane is `Complex{Float64}`.
It is not desirable to redefine `Tuple(::Complex)`, so we overload `coord` instead.

```julia
coord(p :: Float64) = (p, 0.0)
coord(p :: Complex) = reim(p)
incircle(0.0, 1.0, complex(0.0, 1.0), complex(.5, .5))
```

Another interesting type for points is `SVector` from [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl).
`Tuple` is already defined for this type, so we can use them readily.

```julia
using StaticArrays
incircle(SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(0.0, 1.0), SVector(.5, .5))
```

In all the examples above, the conversion has no overhead.

## Extensibility

The code generator turns any piece of code that evaluates a homogeneous polynomial
into a robust piece of code that evaluates the sign of the same polynomial.
It only needs hint to group variables into homogeneous groups.

For example, the discriminant of a degree 2 polynomial ``a x^2 + bx + c`` is ``b^2 - 4ac``, which is a homogeneous polynomial in ``a``, ``b`` and ``c``.
With ExactPredicates, you can write

```julia
using ExactPredicates
using ExactPredicates.Codegen

@genpredicate function discriminant(a, b, c)
    Codegen.group!(a, b, c)
    b*b - 4*a*c
end
```

It will define a function `discriminant(a :: Float64, b :: Float64, c :: Float64)` that returns 1, –1 or 0 depending on the sign of the determinant.


## How it works

This package implements the method used by CGAL and described in:

* Olivier Devillers, Sylvain Pion. “Efficient Exact Geometric Predicates for Delaunay Triangulations”. RR-4351, INRIA. 2002. [⟨inria-00072237⟩](https://hal.inria.fr/inria-00072237)
* Guillaume Melquiond, Sylvain Pion. “Formally Certified Floating-Point Filters For Homogeneous Geometric Predicates”. RAIRO, EDP Sciences, 2007, 41, pp. 57-69. [⟨10.1051/ita:2007005⟩](https://dx.doi.org/10.1051/ita:2007005) [⟨inria-00071232v2⟩](https://hal.inria.fr/inria-00071232)
* Andreas Meyer, Sylvain Pion. “FPG: A code generator for fast and certified geometric predicates”. Real Numbers and Computers, Jun 2008, Santiago de Compostela, Spain. pp.47-60. [⟨inria-00344297⟩](https://hal.inria.fr/inria-00344297)

The implementation relies on Julia's facilities for code generation and
evaluates *a priori* the precision of the floating evaluation, relatively to the
magnitude of the input variables.
