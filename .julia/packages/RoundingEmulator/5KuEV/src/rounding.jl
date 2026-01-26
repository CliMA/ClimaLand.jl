using Base: add12, mul12, significand_bits
using Base.Math: ldexp

const SysFloat = Union{Float32, Float64}

# N_min^s : The smallest positive subnormal number (=nextfloat(zero(T)))
# N_max^s : The largest positive subnormal number (=prevfloat(floatmin(T)))
# N_min^n : The smallest positive normal number (=floatmin(T))
# N_max^n : The largest positive normal number (=floatmax(T))

# constants
for T in (Float32, Float64)
    # log_2(N_min^s)
    # N_min^s = 2 * 2^{-precision(T)} * N_min^n
    @eval exponent_smallest_subnormal(::Type{$T}) = $(Int(log2(nextfloat(zero(T)))))
    @eval exponent_product_errorfree_threshold(::Type{$T}) = $(exponent_smallest_subnormal(T) + 2 * significand_bits(T))
    @eval product_errorfree_threshold(::Type{$T}) = $(ldexp(one(T), exponent_product_errorfree_threshold(T)))
    @eval exponent_product_underflow_mult(::Type{$T}) = $(ceil(Int, -exponent_smallest_subnormal(T)//2))
    @eval product_underflow_mult(::Type{$T}) = $(ldexp(one(T), exponent_product_underflow_mult(T)))
    @eval exponent_quotient_errorfree_threshold(::Type{$T}) = $(-exponent_smallest_subnormal(T) - 3 * significand_bits(T))
    @eval quotient_errorfree_threshold(::Type{$T}) = $(ldexp(one(T), exponent_quotient_errorfree_threshold(T)))
    @eval exponent_quotient_underflow_mult(::Type{$T}) = $(2 * significand_bits(T) + 1)
    @eval quotient_underflow_mult(::Type{$T}) = $(ldexp(one(T), exponent_quotient_underflow_mult(T)))
    @eval inverse_smallest_normal(::Type{$T}) = $(ldexp(one(T), precision(T)))
end

"""
    add_up(a, b)

Computes `a + b` with the rounding mode 
[`Base.Rounding.RoundUp`](https://docs.julialang.org/en/v1/base/math/#Base.Rounding.RoundUp).

```jldoctest
julia> add_up(0.1, 0.2)
0.30000000000000004

julia> add_up(10.0^308, 10.0^308)
Inf

julia> add_up(-10.0^308, -10.0^308)
-1.7976931348623157e308

julia> add_up(-0.1, 0.1)
0.0

julia> add_up(0.0, 0.0)
0.0

julia> add_up(0.0, -0.0)
0.0

julia> add_up(-0.0, -0.0)
-0.0
```
"""
function add_up(a::T, b::T) where {T<:SysFloat}
    x, y = add12(a, b) # twosum
    if isinf(x)
        ifelse(x == typemin(x) && isfinite(a) && isfinite(b), -floatmax(x), x)
    else
        y > zero(y) ? nextfloat(x) : x
    end
end

"""
    add_down(a, b)

Computes `a + b` with the rounding mode
[`Base.Rounding.RoundDown`](https://docs.julialang.org/en/v1/base/math/#Base.Rounding.RoundDown).

```jldoctest
julia> add_down(0.1, 0.2)
0.3

julia> add_down(10.0^308, 10.0^308)
1.7976931348623157e308

julia> add_down(-10.0^308, -10.0^308)
-Inf

julia> add_down(-0.1, 0.1)
-0.0

julia> add_down(0.0, 0.0)
0.0

julia> add_down(0.0, -0.0)
-0.0

julia> add_down(-0.0, -0.0)
-0.0
```
"""
function add_down(a::T, b::T) where {T<:SysFloat}
    x, y = add12(a, b) # twosum
    if isinf(x)
        ifelse(x == typemax(x) && isfinite(a) && isfinite(b), floatmax(x), x)
    elseif y < zero(y)
        prevfloat(x)
    else
        ifelse(iszero(x) && (signbit(a) || signbit(b)), -zero(x), x)
    end
end

"""
    sub_up(a, b)

Computes `a - b` with the rounding mode
[`Base.Rounding.RoundUp`](https://docs.julialang.org/en/v1/base/math/#Base.Rounding.RoundUp).

```jldoctest
julia> sub_up(-0.1, 0.2)
-0.3

julia> sub_up(-10.0^308, 10.0^308)
-1.7976931348623157e308

julia> sub_up(10.0^308, -10.0^308)
Inf

julia> sub_up(0.1, 0.1)
0.0

julia> sub_up(0.0, 0.0)
0.0

julia> sub_up(0.0, -0.0)
0.0

julia> sub_up(-0.0, 0.0)
-0.0

julia> sub_up(-0.0, -0.0)
0.0
```
"""
sub_up(a::T, b::T) where {T<:SysFloat} = add_up(a, -b)

"""
    sub_down(a, b)

Computes `a - b` with the rounding mode
[`Base.Rounding.RoundDown`](https://docs.julialang.org/en/v1/base/math/#Base.Rounding.RoundDown).

```jldoctest
julia> sub_down(-0.1, 0.2)
-0.30000000000000004

julia> sub_down(-10.0^308, 10.0^308)
-Inf

julia> sub_down(10.0^308, -10.0^308)
1.7976931348623157e308

julia> sub_down(0.1, 0.1)
-0.0

julia> sub_down(0.0, 0.0)
-0.0

julia> sub_down(0.0, -0.0)
0.0

julia> sub_down(-0.0, 0.0)
-0.0

julia> sub_down(-0.0, -0.0)
-0.0
```
"""
sub_down(a::T, b::T) where {T<:SysFloat} = add_down(a, -b)

"""
    mul_up(a, b)

Computes `a * b` with the rounding mode
[`Base.Rounding.RoundUp`](https://docs.julialang.org/en/v1/base/math/#Base.Rounding.RoundUp).

```jldoctest
julia> mul_up(0.1, 0.2)
0.020000000000000004

julia> mul_up(10.0^308, 10.0^308)
Inf

julia> mul_up(10.0^308, -10.0^308)
-1.7976931348623157e308

julia> mul_up(5.0e-324, 5.0e-324)
5.0e-324

julia> mul_up(-0.1, 0.1)
-0.01

julia> mul_up(0.0, 0.0)
0.0

julia> mul_up(0.0, -0.0)
-0.0

julia> mul_up(-0.0, -0.0)
0.0
```
"""
function mul_up(a::T, b::T) where {T<:SysFloat}
    x, y = mul12(a, b)
    if isinf(x)
        ifelse(x == typemin(x) && isfinite(a) && isfinite(b), -floatmax(x), x)
    elseif abs(x) > product_errorfree_threshold(T) # not zero(x): (a, b) = (-2.1634867667116802e-200, 1.6930929484402486e-119) fails
        y > zero(y) ? nextfloat(x) : x
    else
        mult = product_underflow_mult(T)
        s, s2 = mul12(a * mult, b * mult)
        t = (x * mult) * mult
        t < s || (t == s && s2 > zero(s2)) ? nextfloat(x) : x
    end
end

"""
    mul_down(a, b)

Computes `a * b` with the rounding mode
[`Base.Rounding.RoundDown`](https://docs.julialang.org/en/v1/base/math/#Base.Rounding.RoundDown).

```jldoctest
julia> mul_down(0.1, 0.2)
0.02

julia> mul_down(10.0^308, 10.0^308)
1.7976931348623157e308

julia> mul_down(10.0^308, -10.0^308)
-Inf

julia> mul_down(5.0e-324, 5.0e-324)
0.0

julia> mul_down(-0.1, 0.1)
-0.010000000000000002

julia> mul_down(0.0, 0.0)
0.0

julia> mul_down(0.0, -0.0)
-0.0

julia> mul_down(-0.0, -0.0)
0.0
```
"""
function mul_down(a::T, b::T) where {T<:SysFloat}
    x, y = mul12(a, b)
    if isinf(x)
        ifelse(x == typemax(x) && isfinite(a) && isfinite(b), floatmax(x), x)
    elseif abs(x) > product_errorfree_threshold(T) # not zero(x): (a, b) = (6.640350825165134e-116, -1.1053488936824272e-202) fails
        y < zero(y) ? prevfloat(x) : x
    else
        mult = product_underflow_mult(T)
        s, s2 = mul12(a * mult, b * mult)
        t = (x * mult) * mult 
        t > s || (t == s && s2 < zero(s2)) ? prevfloat(x) : x
    end
end

"""
    div_up(a, b)

Computes `a / b` with the rounding mode
[`Base.Rounding.RoundUp`](https://docs.julialang.org/en/v1/base/math/#Base.Rounding.RoundUp).

```jldoctest
julia> div_up(0.1, 0.3)
0.33333333333333337

julia> div_up(2.0^-100, 2.0^1000)
5.0e-324

julia> div_up(-0.0, 1.0)
-0.0
```
"""
function div_up(a::T, b::T) where {T<:SysFloat}
    if iszero(a) || iszero(b) || isinf(a) || isinf(b) || isnan(a) || isnan(b)
        a / b
    else
        # if b < 0, flip sign of a and b
        a = flipsign(a, b)
        b = abs(b)
        if abs(a) < product_errorfree_threshold(T) && abs(b) < quotient_errorfree_threshold(T)
            mult = quotient_underflow_mult(T)
            a *= mult
            b *= mult
        end
        d = a / b
        x, y = mul12(d, b)
        x < a || (x == a && y < zero(y)) ? nextfloat(d) : d
    end
end

"""
    div_down(a, b)

Computes `a / b` with the rounding mode
[`Base.Rounding.RoundDown`](https://docs.julialang.org/en/v1/base/math/#Base.Rounding.RoundDown).

```jldoctest
julia> div_down(0.1, 0.3)
0.3333333333333333

julia> div_down(2.0^-100, 2.0^1000)
0.0

julia> div_down(-0.0, 1.0)
-0.0
```
"""
function div_down(a::T, b::T) where {T<:SysFloat}
    if iszero(a) || iszero(b) || isinf(a) || isinf(b) || isnan(a) || isnan(b)
        a / b
    else
        # if b < 0, flip sign of a and b
        a = flipsign(a, b)
        b = abs(b)
        if abs(a) < product_errorfree_threshold(T) && abs(b) < quotient_errorfree_threshold(T)
                mult = quotient_underflow_mult(T)
                a *= mult
                b *= mult
        end
        d = a / b
        x, y = mul12(d, b)
        x > a || (x == a && y > zero(y)) ? prevfloat(d) : d
    end
end

"""
    sqrt_up(a)

Computes `sqrt(a)` with the rounding mode
[`Base.Rounding.RoundUp`](https://docs.julialang.org/en/v1/base/math/#Base.Rounding.RoundUp).

```jldoctest
julia> sqrt_up(2.0)
1.4142135623730951

julia> sqrt_up(0.0)
0.0

julia> sqrt_up(-0.0)
-0.0
```
"""
function sqrt_up(a::SysFloat)
    d = sqrt(a)
    if isinf(d)
        typemax(d)
    elseif a < product_errorfree_threshold(typeof(a))
        invn = inverse_smallest_normal(typeof(a))
        a2 = a * invn^2
        d2 = d * invn
        x, y = mul12(d2, d2)
        x < a2 || (x == a2 && y < zero(y)) ? nextfloat(d) : d
    else
        x, y = mul12(d, d)
        x < a || (x == a  && y < zero(y)) ? nextfloat(d) : d
    end
end

"""
    sqrt_down(a)

Computes `sqrt(a)` with the rounding mode
[`Base.Rounding.RoundDown`](https://docs.julialang.org/en/v1/base/math/#Base.Rounding.RoundDown).

```jldoctest
julia> sqrt_down(2.0)
1.414213562373095

julia> sqrt_down(0.0)
0.0

julia> sqrt_down(-0.0)
-0.0
```
"""
function sqrt_down(a::SysFloat)
    d = sqrt(a)
    if isinf(d)
        typemax(d)
    elseif a < product_errorfree_threshold(typeof(a))
        invn = inverse_smallest_normal(typeof(a))
        a2 = a * invn^2
        d2 = d * invn
        x, y = mul12(d2, d2)
        x > a2 || (x == a2 && y > zero(y)) ? prevfloat(d) : d
    else
        x, y = mul12(d, d)
        x > a || (x == a  && y > zero(y)) ? prevfloat(d) : d
    end
end
