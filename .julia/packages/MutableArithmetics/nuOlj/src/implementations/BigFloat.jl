# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

# This file contains methods to implement the MutableArithmetics API for
# Base.BigFloat.

mutability(::Type{BigFloat}) = IsMutable()

# These methods are copied from `deepcopy_internal` in `base/mpfr.jl`. We don't
# use `mutable_copy(x) = deepcopy(x)` because this creates an empty `IdDict()`
# which costs some extra allocations. We don't need the IdDict case because we
# never call `mutable_copy` recursively.
@static if VERSION >= v"1.12.0-DEV.1343"
    mutable_copy(x::BigFloat) = Base.MPFR._BigFloat(copy(getfield(x, :d)))
else
    function mutable_copy(x::BigFloat)
        d = x._d
        GC.@preserve d begin
            d′ = unsafe_string(pointer(d), sizeof(d))
            return Base.MPFR._BigFloat(x.prec, x.sign, x.exp, d′)
        end
    end
end

const _MPFRRoundingMode = Base.MPFR.MPFRRoundingMode

# copy

promote_operation(::typeof(copy), ::Type{BigFloat}) = BigFloat

function operate_to!(out::BigFloat, ::typeof(copy), in::BigFloat)
    ccall(
        (:mpfr_set, :libmpfr),
        Int32,
        (Ref{BigFloat}, Ref{BigFloat}, _MPFRRoundingMode),
        out,
        in,
        Base.MPFR.ROUNDING_MODE[],
    )
    return out
end

operate!(::typeof(copy), x::BigFloat) = x

# zero

promote_operation(::typeof(zero), ::Type{BigFloat}) = BigFloat

function _set_si!(x::BigFloat, value)
    ccall(
        (:mpfr_set_si, :libmpfr),
        Int32,
        (Ref{BigFloat}, Clong, _MPFRRoundingMode),
        x,
        value,
        Base.MPFR.ROUNDING_MODE[],
    )
    return x
end
operate!(::typeof(zero), x::BigFloat) = _set_si!(x, 0)

# one

promote_operation(::typeof(one), ::Type{BigFloat}) = BigFloat

operate!(::typeof(one), x::BigFloat) = _set_si!(x, 1)

# +

function promote_operation(::typeof(+), ::Type{BigFloat}, ::Type{BigFloat})
    return BigFloat
end

function operate_to!(output::BigFloat, ::typeof(+), a::BigFloat, b::BigFloat)
    ccall(
        (:mpfr_add, :libmpfr),
        Int32,
        (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, _MPFRRoundingMode),
        output,
        a,
        b,
        Base.MPFR.ROUNDING_MODE[],
    )
    return output
end

operate_to!(out::BigFloat, ::typeof(+), a::BigFloat) = operate_to!(out, copy, a)

operate!(::typeof(+), a::BigFloat) = a

# -

promote_operation(::typeof(-), ::Vararg{Type{BigFloat},N}) where {N} = BigFloat

function operate_to!(output::BigFloat, ::typeof(-), a::BigFloat, b::BigFloat)
    ccall(
        (:mpfr_sub, :libmpfr),
        Int32,
        (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, _MPFRRoundingMode),
        output,
        a,
        b,
        Base.MPFR.ROUNDING_MODE[],
    )
    return output
end

function operate!(::typeof(-), x::BigFloat)
    x.sign = -x.sign
    return x
end

function operate_to!(o::BigFloat, ::typeof(-), x::BigFloat)
    operate_to!(o, copy, x)
    return operate!(-, o)
end

# Base.abs

function operate!(::typeof(Base.abs), x::BigFloat)
    x.sign = abs(x.sign)
    return x
end

function operate_to!(o::BigFloat, ::typeof(abs), x::BigFloat)
    operate_to!(o, copy, x)
    return operate!(abs, o)
end

# *

promote_operation(::typeof(*), ::Type{BigFloat}, ::Type{BigFloat}) = BigFloat

function operate_to!(output::BigFloat, ::typeof(*), a::BigFloat, b::BigFloat)
    ccall(
        (:mpfr_mul, :libmpfr),
        Int32,
        (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, _MPFRRoundingMode),
        output,
        a,
        b,
        Base.MPFR.ROUNDING_MODE[],
    )
    return output
end

operate_to!(out::BigFloat, ::typeof(*), a::BigFloat) = operate_to!(out, copy, a)

operate!(::typeof(*), a::BigFloat) = a

# Base.fma

function promote_operation(
    ::typeof(Base.fma),
    ::Type{F},
    ::Type{F},
    ::Type{F},
) where {F<:BigFloat}
    return F
end

function operate_to!(
    output::F,
    ::typeof(Base.fma),
    x::F,
    y::F,
    z::F,
) where {F<:BigFloat}
    ccall(
        (:mpfr_fma, :libmpfr),
        Int32,
        (Ref{F}, Ref{F}, Ref{F}, Ref{F}, _MPFRRoundingMode),
        output,
        x,
        y,
        z,
        Base.MPFR.ROUNDING_MODE[],
    )
    return output
end

function operate!(::typeof(Base.fma), x::F, y::F, z::F) where {F<:BigFloat}
    return operate_to!(x, Base.fma, x, y, z)
end

# Base.muladd

function promote_operation(
    ::typeof(Base.muladd),
    ::Type{F},
    ::Type{F},
    ::Type{F},
) where {F<:BigFloat}
    return F
end

function operate_to!(
    output::F,
    ::typeof(Base.muladd),
    x::F,
    y::F,
    z::F,
) where {F<:BigFloat}
    return operate_to!(output, Base.fma, x, y, z)
end

function operate!(::typeof(Base.muladd), x::F, y::F, z::F) where {F<:BigFloat}
    return operate!(Base.fma, x, y, z)
end

function operate_to!(
    output::BigFloat,
    op::Union{typeof(+),typeof(-),typeof(*)},
    a::BigFloat,
    b::BigFloat,
    c::Vararg{BigFloat,N},
) where {N}
    operate_to!(output, op, a, b)
    return operate!(op, output, c...)
end

function operate!(op::Function, x::BigFloat, args::Vararg{Any,N}) where {N}
    return operate_to!(x, op, x, args...)
end

# add_mul and sub_mul

# Buffer to hold the product
buffer_for(::AddSubMul, args::Vararg{Type{BigFloat},N}) where {N} = BigFloat()

function operate_to!(
    output::BigFloat,
    op::AddSubMul,
    x::BigFloat,
    y::BigFloat,
    z::BigFloat,
    args::Vararg{BigFloat,N},
) where {N}
    return buffered_operate_to!(BigFloat(), output, op, x, y, z, args...)
end

function buffered_operate_to!(
    buffer::BigFloat,
    output::BigFloat,
    op::AddSubMul,
    a::BigFloat,
    x::BigFloat,
    y::BigFloat,
    args::Vararg{BigFloat,N},
) where {N}
    operate_to!(buffer, *, x, y, args...)
    return operate_to!(output, add_sub_op(op), a, buffer)
end

function buffered_operate!(
    buffer::BigFloat,
    op::AddSubMul,
    x::BigFloat,
    args::Vararg{Any,N},
) where {N}
    return buffered_operate_to!(buffer, x, op, x, args...)
end

_scaling_to_bigfloat(x) = _scaling_to(BigFloat, x)

function operate_to!(
    output::BigFloat,
    op::Union{typeof(+),typeof(-),typeof(*)},
    args::Vararg{Scaling,N},
) where {N}
    return operate_to!(output, op, _scaling_to_bigfloat.(args)...)
end

function operate_to!(
    output::BigFloat,
    op::AddSubMul,
    x::Scaling,
    y::Scaling,
    z::Scaling,
    args::Vararg{Scaling,N},
) where {N}
    return operate_to!(
        output,
        op,
        _scaling_to_bigfloat(x),
        _scaling_to_bigfloat(y),
        _scaling_to_bigfloat(z),
        _scaling_to_bigfloat.(args)...,
    )
end

# Called for instance if `args` is `(v', v)` for a vector `v`.
function operate_to!(
    output::BigFloat,
    op::AddSubMul,
    x,
    y,
    z,
    args::Vararg{Any,N},
) where {N}
    return operate_to!(output, add_sub_op(op), x, *(y, z, args...))
end

struct DotBuffer{F<:Real}
    c::F
    t::F
    input::F
    tmp::F

    DotBuffer{F}() where {F<:Real} = new{F}(zero(F), zero(F), zero(F), zero(F))
end

function buffer_for(
    ::typeof(LinearAlgebra.dot),
    ::Type{V},
    ::Type{V},
) where {V<:AbstractVector{BigFloat}}
    return DotBuffer{BigFloat}()
end

# Dot product using Kahan-Babuška-Neumaier compensated summation.
#
# Currently restricted to BigFloat, but easily extendable to other
# similar types, if any such types appear in the ecosystem.
#
# Neumaier's 1974 paper is in German and not currently available from
# Google Scholar, but a scanned version is here:
# https://www.mat.univie.ac.at/~neum/scan/01.pdf
#
# Neumaier 1974:
#   * Rundungsfehleranalyse einiger Verfahren zur Summation endlicher
#     Summen
#   * English name: Rounding Error Analysis of Some Methods for Summing
#     Finite Sums
#   * DOI: https://doi.org/10.1002%2Fzamm.19740540106
#
# The paper has the pseudocode on page two of the PDF, in section 2,
# subsection "IV. Verbessertes Kahan-Babuška-Verfahren".
#
# Wikipedia page and section:
# https://en.wikipedia.org/w/index.php?title=Kahan_summation_algorithm&oldid=1114844162#Further_enhancements
#
# TODO: further improvement: Ogita, Rump and Oishi have a compensated
#       algorithm specifically tailored for dot product, so it would
#       be better to use that (I found out only after writing the code
#       here)
#
# Pseudocode as in Neumaier's 1974 paper:
#
#   s_0 := 0
#
#   # A running compensation for lost low-order bits.
#   w_0 := 0
#
#   for m ∈ 1:n
#     s_m := a_m + s_(m-1)
#
#     if abs(a_m) ≤ abs(s_(m-1))
#       w_m := w_(m-1) + (a_m + (s_(m-1) - s_m))
#     else
#       w_m := w_(m-1) + (s_(m-1) + (a_m - s_m))
#     end
#   end
#
#   # The result, with the correction only applied once in the very
#   # end.
#   s := s_n + w_n
#
# Pseudocode as in the Wikipedia page on Kahan summation (equivalent to
# Neumaier, but closer to the implementation):
#
#   function KahanBabushkaNeumaierSum(input)
#       sum = 0.0
#       # A running compensation for lost low-order bits.
#       c = 0.0
#       for i ∈ eachindex(input)
#           t = sum + input[i]
#           if abs(input[i]) ≤ abs(sum)
#               tmp = (sum - t) + input[i]
#           else
#               tmp = (input[i] - t) + sum
#           end
#           c += tmp
#           sum = t
#       end
#       # The result, with the correction only applied once in the very end.
#       sum + c
#   end

# Returns abs(x) <= abs(y) without allocating.
function _abs_lte_abs(x::BigFloat, y::BigFloat)
    x_is_neg, y_is_neg = signbit(x), signbit(y)
    if x_is_neg != y_is_neg
        operate!(-, x)
    end
    ret = y_is_neg ? y <= x : x <= y
    if x_is_neg != y_is_neg
        operate!(-, x)
    end
    return ret
end

function buffered_operate_to!(
    buf::DotBuffer{BigFloat},
    sum::BigFloat,
    ::typeof(LinearAlgebra.dot),
    x::AbstractVector{BigFloat},
    y::AbstractVector{BigFloat},
)                                                   # See pseudocode description
    operate!(zero, sum)                             # sum = 0
    operate!(zero, buf.c)                           # c = 0
    for (xi, yi) in zip(x, y)                       # for i in eachindex(input)
        operate_to!(buf.input, copy, xi)            # input = x[i]
        operate!(*, buf.input, yi)                  # input = x[i] * y[i]
        operate_to!(buf.t, +, sum, buf.input)       # t = sum + input
        if _abs_lte_abs(buf.input, sum)             # if |input| < |sum|
            operate_to!(buf.tmp, copy, sum)         # tmp = sum
            operate!(-, buf.tmp, buf.t)             # tmp = sum - t
            operate!(+, buf.tmp, buf.input)         # tmp = (sum - t) + input
        else
            operate_to!(buf.tmp, copy, buf.input)   # tmp = input
            operate!(-, buf.tmp, buf.t)             # tmp = input - t
            operate!(+, buf.tmp, sum)               # tmp = (input - t) + sum
        end
        operate!(+, buf.c, buf.tmp)                 # c += tmp
        operate_to!(sum, copy, buf.t)               # sum = t
    end
    operate!(+, sum, buf.c)                         # sum += c
    return sum
end

# Base.evalpoly

function operate!(
    op::typeof(Base.evalpoly),
    out::F,
    coefs::Union{Tuple,AbstractVector},
) where {F<:BigFloat}
    return operate_to!(out, op, mutable_copy(out), coefs)
end

function operate(
    op::typeof(Base.evalpoly),
    x::F,
    coefs::Union{Tuple,AbstractVector},
) where {F<:BigFloat}
    return operate_to!(zero(x), op, x, coefs)
end
