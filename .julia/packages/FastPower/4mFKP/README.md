# FastPower.jl

A faster approximation to floating point power, at the trade-off of some accuracy. While Julia's
built-in floating point `^` tries to achieve ~1ulp accuracy, this version of floating point power
approximation achieves much fewer digits of accuracy (approximately 10 digits of accuracy) while
being much faster. This is developed as a library in order to make the choice to opt-in as a
replacement to `^` very easy but deliberate on the side of the user.

## Installation

```julia
using Pkg
Pkg.add("FastPower")
```

## Using FastPower.jl

Using FastPower.jl is dead simple: instead of `x^y`, do the following:

```julia
using FastPower
FastPower.fastpower(x,y)
```

That's it. That's all there is. 

## FastPower vs FastMath (`@fastmath`)

The name simply derives from the Julia standard of `x_fast` for things that are approximations.
FastPower is simply the the `^_fast` or `pow_fast` function, following the standard conventions
developed from Base. However, this differs from the `pow_fast` you get from Base which is still
a lot more accurate. `FastPower.fastpower` loses about 12 digits of accuracy on Float64, so it's
about 3-4 digits of accuracy. For many applications, such as the adaptivity algorithm when 
solving differential equations, this can be a sufficient amount of accuracy for a power 
function approximation.

This approximation is tested for the range of `x^y` where `x>=0` and `y>=0`. If `x<=1`, then
the approximation is accurate for very large values of `y`. If `x<100`, then the approximation
is accurate for `y<1`. If `x>100`, or if `x` is large and `y` is large, caution should be used. 

## What about FastPow.jl?

These two packages are completely unrelated since [FastPow.jl](https://github.com/JuliaMath/FastPow.jl) is a specialization for *literal* integer powers: powers that are not only
integers but appear as literal constants in the source code.
It does things like:

```julia
x^5
```

which is transformed by the `@fastpow` macro to be computed via:

```julia
sq = x^2
fourth = sq^2
fourth * 2
```

This is faster than `^(::AbstractFloat, Integer)` but with a bit of accuracy loss compared to
what LLVM generates by default for `x^5`.

Meanwhile, FastPower.jl is all about floating-point powers (whose value may only be known at runtime). 

## Why is this not in Base?

Maybe it could be. If you wish to change `pow_fast` to this, open a PR. There can be a debate
as to which one is better. However, as separate package there is no debate: use this one if
the less accuracy and more speed is appropriate for your needs. Arguably, FastMath should be
split from Base.

Also, one major purpose of this package is to allow for the bithacking to be held in a place 
that allows for extensions which enable compatbility with automatic differentiation. Extensions
for:

* Enzyme
* ForwardDiff
* Tracker
* ReverseDiff
* MonteCarloMeasurements
* Measurements

currently exist for this function. More can be added on-demand. This allows for `pow_fast` to be
safe in most AD contexts, though in some cases improved extensions could be created which
improve performance.
