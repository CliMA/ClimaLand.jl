"""
`Ratios` provides `SimpleRatio`, a faster variant of `Rational`.
Speed comes at the cost of greater vulnerability to overflow.

API summary:

- `r = SimpleRatio(num, den)` is equivalent to `num // den`
- `common_denominator` standardizes a collection of `SimpleRatio`s to have the same denominator,
  making some arithmetic operations among them less likely to overflow.
"""
module Ratios

import Base: convert, promote_rule, *, /, +, -, ^, ==, decompose, isinteger


export SimpleRatio, common_denominator

struct SimpleRatio{T<:Integer} <: Real
    num::T
    den::T
end

"""
    SimpleRatio(num::Integer, den::Integer)

Construct the equivalent of the rational number `num // den`.

Operations with `SimpleRatio` are faster, but also more vulnerable to integer overflow,
than with `Rational`. Arithmetic with `SimpleRatio` does not perform any simplification of the resulting ratios.
The best defense against overflow is to use ratios with the same denominator: in such cases,
`+` and `-` will skip forming the product of denominators. See [`common_denominator`](@ref).

If overflow is a risk, consider constructing them using SaferIntegers.jl.

# Examples

```julia
julia> x, y, z = SimpleRatio(1, 8), SimpleRatio(1, 4), SimpleRatio(2, 8)
(SimpleRatio{$Int}(1, 8), SimpleRatio{$Int}(1, 4), SimpleRatio{$Int}(2, 8))

julia> x+y
SimpleRatio{$Int}(12, 32)

julia> x+z
SimpleRatio{$Int}(3, 8)
```
"""
SimpleRatio(num::Integer, den::Integer) = SimpleRatio(promote(num, den)...)

convert(::Type{BigFloat}, r::SimpleRatio{S}) where {S} = BigFloat(r.num)/r.den
function convert(::Type{T}, r::SimpleRatio{S}) where {T<:AbstractFloat,S}
    P = promote_type(T,S)
    convert(T, convert(P, r.num)/convert(P, r.den))
end
(::Type{T})(r::SimpleRatio) where T<:AbstractFloat = convert(T, r)

SimpleRatio{T}(i::Integer) where {T<:Integer} = SimpleRatio{T}(convert(T, i), oneunit(T))
SimpleRatio{T}(r::Rational{S}) where {T<:Integer, S<:Integer} = SimpleRatio(convert(T, r.num), convert(T, r.den))
Rational{T}(r::SimpleRatio{S}) where {T<:Integer, S<:Integer} = convert(T, r.num) // convert(T, r.den)

*(x::SimpleRatio, y::SimpleRatio) = SimpleRatio(x.num*y.num, x.den*y.den)
*(x::SimpleRatio, y::Bool) = SimpleRatio(x.num*y, x.den)
*(x::SimpleRatio, y::Integer) = SimpleRatio(x.num*y, x.den)
*(x::Bool, y::SimpleRatio) = SimpleRatio(x*y.num, y.den)
*(x::Integer, y::SimpleRatio) = SimpleRatio(x*y.num, y.den)
/(x::SimpleRatio, y::SimpleRatio) = SimpleRatio(x.num*y.den, x.den*y.num)
/(x::SimpleRatio, y::Integer) = SimpleRatio(x.num, x.den*y)
/(x::Integer, y::SimpleRatio) = SimpleRatio(x*y.den, y.num)
+(x::Integer, y::SimpleRatio) = SimpleRatio(x*y.den + y.num, y.den)
-(x::Integer, y::SimpleRatio) = SimpleRatio(x*y.den - y.num, y.den)
+(x::SimpleRatio, y::SimpleRatio) = x.den == y.den ? SimpleRatio(x.num + y.num, x.den) :
                                                     SimpleRatio(x.num*y.den + x.den*y.num, x.den*y.den)
-(x::SimpleRatio, y::SimpleRatio) = x.den == y.den ? SimpleRatio(x.num - y.num, x.den) :
                                                     SimpleRatio(x.num*y.den - x.den*y.num, x.den*y.den)
^(x::SimpleRatio, y::Integer) = SimpleRatio(x.num^y, x.den^y)

-(x::SimpleRatio{T}) where {T<:Signed} = SimpleRatio(-x.num, x.den)
-(x::SimpleRatio{T}) where {T<:Unsigned} = throw(VERSION < v"0.7.0-DEV.1269" ? OverflowError() : OverflowError("cannot negate unsigned number"))

promote_rule(::Type{SimpleRatio{T}}, ::Type{S}) where {T<:Integer,S<:Integer} = SimpleRatio{promote_type(T,S)}
promote_rule(::Type{SimpleRatio{T}}, ::Type{SimpleRatio{S}}) where {T<:Integer,S<:Integer} = SimpleRatio{promote_type(T,S)}
promote_rule(::Type{SimpleRatio{T}}, ::Type{S}) where {T<:Integer,S<:AbstractFloat} = promote_type(T,S)
promote_rule(::Type{SimpleRatio{T}}, ::Type{Rational{S}}) where {T<:Integer,S<:Integer} = Rational{promote_type(T,S)}

==(x::SimpleRatio, y::SimpleRatio) = x.num*y.den == x.den*y.num

==(x::SimpleRatio, y::Integer) = x.num == x.den*y
==(x::Integer, y::SimpleRatio) = x*y.den == y.num

function ==(x::AbstractFloat, q::SimpleRatio)
    if isfinite(x)
        (count_ones(q.den) == 1) & (x*q.den == q.num)
    else
        x == q.num/q.den
    end
end

==(q::SimpleRatio, x::AbstractFloat) = x == q

decompose(x::SimpleRatio) = x.num, 0, x.den

isinteger(x::SimpleRatio) = gcd(x.num, x.den) == abs(x.den)

"""
    common_denominator(x::SimpleRatio, ys::SimpleRatio...)

Return the equivalent of `(x, ys...)` but using a common denominator.
This can be useful to avoid overflow.

!!! info
    This function is quite slow.  In performance-sensitive contexts where the
    ratios are constructed with literal constants, it is better to ensure a common
    denominator at the time of original construction.
"""
function common_denominator(x::SimpleRatio, ys::SimpleRatio...)
    all(y -> y.den == x.den, ys) && return (x, ys...)
    cd = gcd(x.den, map(y -> y.den, ys)...)        # common divisor
    # product of the denominators
    pd = Base.Checked.checked_mul(cd, mapreduce(z -> z.den ÷ cd, Base.Checked.checked_mul, (x, ys...)))
    return map((x, ys...)) do z
        SimpleRatio(z.num * (pd ÷ z.den), pd)
    end
end


if !isdefined(Base, :get_extension)
    using Requires
end

@static if !isdefined(Base, :get_extension)
function __init__()
    @require FixedPointNumbers = "53c48c17-4a7d-5ca2-90c5-79b7896eea93" begin
        include("../ext/RatiosFixedPointNumbersExt.jl")
    end
end
end

end
