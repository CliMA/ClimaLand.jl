# ColorVectorSpace

[![Build Status](https://github.com/JuliaGraphics/ColorVectorSpace.jl/workflows/Unit%20test/badge.svg)](https://github.com/JuliaGraphics/ColorVectorSpace.jl/actions)
[![codecov](https://codecov.io/github/JuliaGraphics/ColorVectorSpace.jl/graph/badge.svg?token=NFLslfrU3p)](https://codecov.io/github/JuliaGraphics/ColorVectorSpace.jl)

This package is an add-on to [ColorTypes](https://github.com/JuliaGraphics/ColorTypes.jl), and provides fast
mathematical operations for objects with types such as `RGB` and
`Gray`.
Specifically, with this package both grayscale and `RGB` colors are treated as if they are points
in a normed vector space.

## Introduction

Colorspaces such as RGB, unlike XYZ, are technically non-linear;
perhaps the most "colorimetrically correct" approach when averaging two RGBs is to
first convert each to XYZ, average them, and then convert back to RGB.
Nor is there a clear definition of computing the sum of two colors.
As a consequence, Julia's base color package,
[ColorTypes](https://github.com/JuliaGraphics/ColorTypes.jl),
does not support mathematical operations on colors.

However, particularly in image processing it is common to ignore this
concern, and for the sake of performance treat an RGB as if it were a
3-vector.  The role of this package is to extend ColorTypes to support such mathematical operations.
Specifically, it defines `+` and multiplication by a scalar (and by extension, `-` and division by a scalar) for grayscale and `AbstractRGB` colors.
These are the requirements of a [vector space](https://en.wikipedia.org/wiki/Vector_space).

If you're curious about how much the "colorimetrically correct" and
"vector space" views differ, the following
diagram might help. The first 10 `distinguishable_colors` were
generated, and all pairs were averaged. Each box represents the
average of the pair of diagonal elements intersected by tracing
vertically and horizontally; within each box, the upper diagonal is
the "colorimetrically-correct" version, while the lower diagonal
represents the "RGB vector space" version.

![ColorVectorSpace](images/comparison.png "Comparison")

This package also defines `norm(c)` for RGB and grayscale colors.
This makes these color spaces [normed vector spaces](https://en.wikipedia.org/wiki/Normed_vector_space).
Note that `norm` has been designed to satisfy **equivalence** of grayscale and RGB representations: if
`x` is a scalar, then `norm(x) == norm(Gray(x)) == norm(RGB(x, x, x))`.
Effectively, there's a division-by-3 in the `norm(::RGB)` case compared to the Euclidean interpretation of
the RGB vector space.
Equivalence is an important principle for the Colors ecosystem, and violations should be reported as likely bugs.
One violation is `abs2`; see the section below for more detail.

## Usage

```julia
using ColorTypes, ColorVectorSpace
```

For the most part, that's it; just by loading `ColorVectorSpace`, most basic mathematical
operations will "just work" on `AbstractRGB`, `AbstractGray`
(`Color{T,1}`), `TransparentRGB`, and `TransparentGray` objects.
(See definitions for the latter inside of `ColorTypes`).

However, there are some additional operations that you may need to distinguish carefully.

### Multiplication

Grayscale values are conceptually similar to scalars, and consequently it seems straightforward to define multiplication of two grayscale values.
RGB values present more options.
This package supports three different notions of multiplication: the inner product, the hadamard (elementwise) product, and the tensor product.

```julia
julia> c1, c2 = RGB(0.2, 0.3, 0.4), RGB(0.5, 0.3, 0.2)
(RGB{Float64}(0.2,0.3,0.4), RGB{Float64}(0.5,0.3,0.2))

julia> c1⋅c2     # \cdot<TAB> # or dot(c1, c2)
0.09000000000000001

# This is equivelant to `mapc(*, c1, c2)`
julia> c1⊙c2     # \odot<TAB> # or hadamard(c1, c2)
RGB{Float64}(0.1,0.09,0.08000000000000002)

julia> c1⊗c2    # \otimes<TAB> # or tensor(c1, c2)
RGBRGB{Float64}:
 0.1   0.06  0.04
 0.15  0.09  0.06
 0.2   0.12  0.08
```

Note that `c1⋅c2 = (c1.r*c2.r + c1.g*c2.g + c1.b*c2.b)/3`, where the division by 3 ensures the equivalence `norm(x) == norm(Gray(x)) == norm(RGB(x, x, x))`.

Ordinary multiplication `*` is not supported because it is not obvious which one of these should be the default option.

However, `*` is defined for grayscale since all these three multiplication operations (i.e., `⋅`, `⊙` and `⊗`) are equivalent in the 1D vector space.

### Variance

The variance `v = E((c - μ)^2)` (or its bias-corrected version) involves a multiplication,
and to be consistent with the above you must specify which sense of multiplication you wish to use:

```julia
julia> cs = [c1, c2]
2-element Array{RGB{Float64},1} with eltype RGB{Float64}:
 RGB{Float64}(0.2,0.3,0.4)
 RGB{Float64}(0.5,0.3,0.2)

julia> varmult(⋅, cs)
0.021666666666666667

julia> varmult(⊙, cs)
RGB{Float64}(0.045,0.0,0.020000000000000004)

julia> varmult(⊗, cs)
RGBRGB{Float64}:
  0.045  0.0  -0.03
  0.0    0.0   0.0
 -0.03   0.0   0.02
```

The corresponding `stdmult` computes standard deviation.

### `abs` and `abs2`

To begin with, there is no general and straightforward definition of the
absolute value of a vector.
There are two reasonably intuitive definitions of `abs`/`abs2`: as a channel-wise
operator or as a function which returns a real number based on the norm.
For the latter, there are also variations in the definition of norm.

In ColorVectorSpace v0.9 and later, `abs` is defined as a channel-wise operator.
`abs2` returns a real-valued scalar. In previous versions of ColorVectorSpace,
for `g = Gray(0.3)`, ColorVectorSpace returned different values for `abs2(g)` and
`abs2(RGB(g))`. This breaks the equivalence of `g` and `RGB(g)`.
This behavior is retained, with a deprecation warning, starting with
ColorVectorSpace 0.9.6.

We anticipate the following transition schedule:

- Sept 21, 2021: release ColorVectorSpace 0.9.7 with both `abs2` and `Future.abs2`.
  `abs2` will have a "quiet" deprecation warning (visible with `--depwarn=yes`
  or when running `Pkg.test`)
- May 19, 2022: make the deprecation warning "noisy" (cannot be turned off).
  This is designed to catch user-level scripts that may need to update thresholds
  or other constants.
- *July 20, 2023: transition `abs2` to the new definition
  and make `Future.abs2` give a "noisy" depwarn to revert to regular `abs2`.
- *Dec 1, 2023: remove `Future.abs2`.

The line marked with `*` denotes a breaking release.
