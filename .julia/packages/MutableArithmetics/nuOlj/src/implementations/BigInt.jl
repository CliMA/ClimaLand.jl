# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

# This file contains methods to implement the MutableArithmetics API for
# Base.BigInt.

mutability(::Type{BigInt}) = IsMutable()

# Copied from `deepcopy_internal` implementation in Julia:
# https://github.com/JuliaLang/julia/blob/7d41d1eb610cad490cbaece8887f9bbd2a775021/base/gmp.jl#L772
mutable_copy(x::BigInt) = Base.GMP.MPZ.set(x)

# copy

promote_operation(::typeof(copy), ::Type{BigInt}) = BigInt

function operate_to!(out::BigInt, ::typeof(copy), in::BigInt)
    Base.GMP.MPZ.set!(out, in)
    return out
end

operate!(::typeof(copy), x::BigInt) = x

# zero

promote_operation(::typeof(zero), ::Type{BigInt}) = BigInt

operate!(::typeof(zero), x::BigInt) = Base.GMP.MPZ.set_si!(x, 0)

# one

promote_operation(::typeof(one), ::Type{BigInt}) = BigInt

operate!(::typeof(one), x::BigInt) = Base.GMP.MPZ.set_si!(x, 1)

# +

promote_operation(::typeof(+), ::Type{BigInt}, ::Type{BigInt}) = BigInt

function operate_to!(output::BigInt, ::typeof(+), a::BigInt, b::BigInt)
    return Base.GMP.MPZ.add!(output, a, b)
end

operate_to!(out::BigInt, ::typeof(+), a::BigInt) = operate_to!(out, copy, a)

operate!(::typeof(+), a::BigInt) = a

# -

promote_operation(::typeof(-), ::Vararg{Type{BigInt},N}) where {N} = BigInt

function operate_to!(output::BigInt, ::typeof(-), a::BigInt, b::BigInt)
    return Base.GMP.MPZ.sub!(output, a, b)
end

operate_to!(out::BigInt, ::typeof(-), a::BigInt) = Base.GMP.MPZ.neg!(out, a)

operate!(::typeof(-), a::BigInt) = operate_to!(a, -, a)

# abs

promote_operation(::typeof(abs), ::Type{BigInt}) = BigInt

function operate!(::typeof(abs), n::BigInt)
    n.size = abs(n.size)
    return n
end

function operate_to!(o::BigInt, ::typeof(abs), n::BigInt)
    operate_to!(o, copy, n)
    return operate!(abs, o)
end

# *

promote_operation(::typeof(*), ::Type{BigInt}, ::Type{BigInt}) = BigInt

function operate_to!(output::BigInt, ::typeof(*), a::BigInt, b::BigInt)
    return Base.GMP.MPZ.mul!(output, a, b)
end

operate_to!(out::BigInt, ::typeof(*), a::BigInt) = operate_to!(out, copy, a)

operate!(::typeof(*), a::BigInt) = a

promote_operation(::typeof(div), ::Type{BigInt}, ::Type{BigInt}) = BigInt

function operate_to!(output::BigInt, ::typeof(div), a::BigInt, b::BigInt)
    return Base.GMP.MPZ.tdiv_q!(output, a, b)
end

# gcd

function promote_operation(
    ::Union{typeof(gcd),typeof(lcm)},
    ::Type{BigInt},
    ::Type{BigInt},
)
    return BigInt
end

function operate_to!(output::BigInt, ::typeof(gcd), a::BigInt, b::BigInt)
    return Base.GMP.MPZ.gcd!(output, a, b)
end

function operate_to!(output::BigInt, ::typeof(lcm), a::BigInt, b::BigInt)
    return Base.GMP.MPZ.lcm!(output, a, b)
end

function operate_to!(
    output::BigInt,
    op::Union{typeof(gcd),typeof(lcm)},
    a::BigInt,
    b::BigInt,
    c::Vararg{BigInt,N},
) where {N}
    operate_to!(output, op, a, b)
    return operate!(op, output, c...)
end

function operate_to!(
    output::BigInt,
    op::Union{typeof(+),typeof(-),typeof(*)},
    a::BigInt,
    b::BigInt,
    c::Vararg{BigInt,N},
) where {N}
    operate_to!(output, op, a, b)
    return operate!(op, output, c...)
end

function operate!(op::Function, x::BigInt, args::Vararg{Any,N}) where {N}
    return operate_to!(x, op, x, args...)
end

# add_mul and sub_mul

# Buffer to hold the product
buffer_for(::AddSubMul, args::Vararg{Type{BigInt},N}) where {N} = BigInt()

function operate_to!(
    output::BigInt,
    op::AddSubMul,
    x::BigInt,
    y::BigInt,
    z::BigInt,
    args::Vararg{BigInt,N},
) where {N}
    return buffered_operate_to!(BigInt(), output, op, x, y, z, args...)
end

function buffered_operate_to!(
    buffer::BigInt,
    output::BigInt,
    op::AddSubMul,
    a::BigInt,
    x::BigInt,
    y::BigInt,
    args::Vararg{BigInt,N},
) where {N}
    operate_to!(buffer, *, x, y, args...)
    return operate_to!(output, add_sub_op(op), a, buffer)
end

function buffered_operate!(
    buffer::BigInt,
    op::AddSubMul,
    x::BigInt,
    args::Vararg{Any,N},
) where {N}
    return buffered_operate_to!(buffer, x, op, x, args...)
end

function _scaling_to(::Type{T}, x) where {T}
    return convert(T, scaling_to_number(x))
end

_scaling_to_bigint(x) = _scaling_to(BigInt, x)

function operate_to!(
    output::BigInt,
    op::Union{
        typeof(+),
        typeof(-),
        typeof(*),
        typeof(div),
        typeof(gcd),
        typeof(lcm),
    },
    args::Vararg{Scaling,N},
) where {N}
    return operate_to!(output, op, _scaling_to_bigint.(args)...)
end

function operate_to!(
    output::BigInt,
    op::AddSubMul,
    x::Scaling,
    y::Scaling,
    z::Scaling,
    args::Vararg{Scaling,N},
) where {N}
    return operate_to!(
        output,
        op,
        _scaling_to_bigint(x),
        _scaling_to_bigint(y),
        _scaling_to_bigint(z),
        _scaling_to_bigint.(args)...,
    )
end

# Called for instance if `args` is `(v', v)` for a vector `v`.
function operate_to!(
    output::BigInt,
    op::AddSubMul,
    x,
    y,
    z,
    args::Vararg{Any,N},
) where {N}
    return operate_to!(output, add_sub_op(op), x, *(y, z, args...))
end
