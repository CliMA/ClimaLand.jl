<h1 align="center">
CRlibm
</h1>

CRlibm.jl is a Julia package wrapping the [CRlibm library](http://lipforge.ens-lyon.fr/www/crlibm/). This library provides correctly-rounded mathematical functions, as described on the library's home page:
- implementations of the double-precision C99 standard elementary functions
- correctly rounded in the four IEEE-754 rounding modes
- with a comprehensive proof of both the algorithms used and their implementation
- sufficiently efficient in average time, worst-case time, and memory consumption to replace existing libms transparently

CRlibm is distributed under the GNU Lesser General Public License (LGPL), while the CRlibm.jl package is distributed under the MIT license.

## Installation

The IntervalArithmetic.jl package requires to [install Julia](https://julialang.org/downloads/) (v1.3 or above).

Then, start Julia and execute the following command in the REPL:

```julia
using Pkg; Pkg.add("CRlibm")
```

## Usage

The floating-point rounding mode must be set to `RoundNearest` for the library to work correctly. Generally, nothing needs to be done, since this is the default value:

```julia
julia> rounding(Float64)
RoundingMode{:Nearest}()
```

The library provides correctly-rounded versions of elementary functions, such as `sin` and `exp` (see [below](#list-of-implemented-functions) for the complete list). The available rounding modes are `RoundNearest`, `RoundUp`, `RoundDown` and `RoundToZero`. The functions are used as follows:

```julia
julia> using CRlibm

julia> CRlibm.cos(0.5, RoundUp)
0.8775825618903728

julia> CRlibm.cos(0.5, RoundDown)
0.8775825618903726

julia> CRlibm.cos(0.5, RoundNearest)
0.8775825618903728

julia> CRlibm.cos(1.6, RoundToZero)
-0.029199522301288812

julia> CRlibm.cos(1.6, RoundDown)
-0.029199522301288815

julia> CRlibm.cos(0.5) # equivalent to `CRlibm.cos(0.5, RoundNearest)`
0.8775825618903728
```

## List of implemented functions

All functions from CRlibm are wrapped, except the power function:
- `exp`, `expm1`
- `log`, `log1p`, `log2`, `log10`
- `sin`, `cos`, `tan`, `asin`, `acos`, `atan`
- `sinpi`, `cospi`, `tanpi`, `atanpi`
- `sinh`, `cosh`

## What is correct rounding?

Suppose that we ask Julia to calculate the cosine of a number

```julia
julia> cos(0.5) # Julia build-in `cos` function
0.8775825618903728
```

using the built-in mathematics library [OpenLibm](https://github.com/JuliaLang/openlibm). The result is a floating-point number that is a very good approximation to the true value. However, we do not know if the result that Julia gives is below or above the true value, nor how far away it is.

Correctly-rounded functions guarantee that when the result is not exactly representable as a floating-point number, the value returned is the **next largest floating-point number when rounding up**, or **the next smallest when rounding down**. This is equivalent to doing the calculation in infinite precision and then performing the rounding.

## Rationale for the Julia wrapper

The `CRlibm` library is state-of-the-art as regards correctly-rounded functions of `Float64` arguments. It is required for the interval arithmetic library [`IntervalArithmetic`](https://github.com/JuliaIntervals/IntervalArithmetic.jl). Having gone to the trouble of wrapping it, it made sense to release it separately; for example, it could be used to test the quality of the OpenLibm functions.

## Lacunae

CRlibm is missing a (guaranteed) correctly-rounded power function, i.e., `x^y`. The fact that there are two arguments, instead of a single argument for functions such as `sin`, means that correct rounding is much harder; see e.g., [1].

[1] P. Kornerup, C. Lauter, V. Lefèvre, N. Louvet and J.-M. Muller, [Computing Correctly Rounded Integer Powers in Floating-Point Arithmetic](http://perso.ens-lyon.fr/jean-michel.muller/p1-Kornerup.pdf), *ACM Transactions on Mathematical Software*, **37** (1), 2010

## MPFR as an alternative to CRlibm

As far as we are aware, the only alternative package to CRlibm is [MPFR](http://www.mpfr.org/), which provides correctly-rounded functions for
floating-point numbers of **arbitrary precision**, including the power function. However, it can be slow.

MPFR is wrapped in base Julia in the `BigFloat` type, and can emulate double-precision floating point by setting the precision to 53 bits. In particular, the functions provided by CRlibm are extended to `BigFloat` via MPFR, with the same syntax:

```julia
julia> setprecision(64) # change `BigFloat` precision
64

julia> CRlibm.exp(BigFloat(0.51), RoundDown)
1.66529119494588632316

julia> CRlibm.exp(BigFloat(0.51), RoundUp)
1.66529119494588632327
```

The function `CRlibm.setup(true)` can be called to redefine the functions to use MPFR instead of CRlibm; this is automatic on 32 bit systems.
