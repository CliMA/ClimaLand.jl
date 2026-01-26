# RoundingEmulator.jl

Emulate directed rounding using only the default rounding mode. 

This package is meant to produce the exact same results of `Rounding.setrounding` ([deprecated](https://github.com/JuliaLang/julia/pull/27166)) without switching rounding modes.

## Requirements 
 - Julia 1.3 or higher
 - `Base.Rounding.get_zero_subnormals() == true`. (See [Base.Rounding.get_zero_subnormals](https://docs.julialang.org/en/v1/base/numbers/#Base.Rounding.get_zero_subnormals))

## Use

This package provides
* [`add_up`](@ref), [`add_down`](@ref) - Addition
* [`sub_up`](@ref), [`sub_down`](@ref) - Subtraction
* [`mul_up`](@ref), [`mul_down`](@ref) - Multiplication
* [`div_up`](@ref), [`div_down`](@ref) - Division
* [`sqrt_up`](@ref), [`sqrt_down`](@ref) - Square root

`up`: Round up,
`down`: Round down

```julia
julia> using RoundingEmulator

julia> add_up(0.1, 0.2)
0.30000000000000004

julia> bitstring(add_up(0.1, 0.2))
"0011111111010011001100110011001100110011001100110011001100110100"

julia> add_down(0.1, 0.2)
0.3

julia> bitstring(add_down(0.1, 0.2))
"0011111111010011001100110011001100110011001100110011001100110011"

julia> sub_up(-0.1, 0.2)
-0.3

julia> bitstring(sub_up(-0.1, 0.2))
"1011111111010011001100110011001100110011001100110011001100110011"

julia> sub_down(-0.1, 0.2)
-0.30000000000000004

julia> bitstring(sub_down(-0.1, 0.2))
"1011111111010011001100110011001100110011001100110011001100110100"

julia> mul_up(0.1, 0.2)
0.020000000000000004

julia> bitstring(mul_up(0.1, 0.2))
"0011111110010100011110101110000101000111101011100001010001111100"

julia> mul_down(0.1, 0.2)
0.02

julia> bitstring(mul_down(0.1, 0.2))
"0011111110010100011110101110000101000111101011100001010001111011"

julia> div_up(1.0, 3.0)
0.33333333333333337

julia> bitstring(div_up(1.0, 3.0))
"0011111111010101010101010101010101010101010101010101010101010110"

julia> div_down(1.0, 3.0)
0.3333333333333333

julia> bitstring(div_down(1.0, 3.0))
"0011111111010101010101010101010101010101010101010101010101010101"

julia> sqrt_up(2.0)
1.4142135623730951

julia> bitstring(sqrt_up(2.0))
"0011111111110110101000001001111001100110011111110011101111001101"

julia> sqrt_down(2.0)
1.414213562373095

julia> bitstring(sqrt_down(2.0))
"0011111111110110101000001001111001100110011111110011101111001100"
```

## Corner cases
```julia
julia> u = nextfloat(zero(Float64))
5.0e-324

julia> v = floatmax(Float64)
1.7976931348623157e308

julia> v + v
Inf

julia> add_up(v, v)
Inf

julia> add_down(v, v)
1.7976931348623157e308

julia> u * u
0.0

julia> mul_up(u, u)
5.0e-324

julia> mul_down(u, u)
0.0

julia> 1.0 / u
Inf

julia> div_up(1.0, u)
Inf

julia> div_down(1.0, u)
1.7976931348623157e308
```

## Signed zero
`RoundingEmulator` follows the special rules for signed zero specified in the chapter 6.3 of IEEE 754-2019.
```julia
julia> add_up(-1.0, 1.0)
0.0

julia> add_down(-1.0, 1.0)
-0.0

julia> add_up(-0.0, 0.0)
0.0

julia> add_down(-0.0, 0.0)
-0.0

julia> add_up(0.0, 0.0)
0.0

julia> add_down(0.0, 0.0)
0.0

julia> sqrt_up(-0.0)
-0.0

julia> sqrt_down(-0.0)
-0.0
```
