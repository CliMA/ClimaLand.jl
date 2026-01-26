# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

# This file contains methods to implement the MutableArithmetics API for
# Base.Rational{T}.

# `Rational` is a `struct`, not a `mutable struct`, so if `T` is not
# mutable, we cannot mutate the rational.
mutability(::Type{Rational{T}}) where {T} = mutability(T)

mutable_copy(x::Rational) = Rational(mutable_copy(x.num), mutable_copy(x.den))

# copy

promote_operation(::typeof(copy), ::Type{Q}) where {Q<:Rational} = Q

function operate_to!(out::Q, ::typeof(copy), in::Q) where {Q<:Rational}
    operate_to!(out.num, copy, in.num)
    operate_to!(out.den, copy, in.den)
    return out
end

operate!(::typeof(copy), x::Rational) = x

# zero

promote_operation(::typeof(zero), ::Type{Rational{T}}) where {T} = Rational{T}

function operate!(::typeof(zero), x::Rational)
    operate!(zero, x.num)
    operate!(one, x.den)
    return x
end

# one

promote_operation(::typeof(one), ::Type{Rational{T}}) where {T} = Rational{T}

function operate!(::typeof(one), x::Rational)
    operate!(one, x.num)
    operate!(one, x.den)
    return x
end

function _buffered_divgcd(buffer, x, x_output, y, y_output)
    operate_to!(buffer, gcd, x, y)
    operate_to!(x_output, div, x, buffer)
    operate_to!(y_output, div, y, buffer)
    return
end

function _buffered_divgcd(buffer, x, y)
    operate_to!(buffer, gcd, x, y)
    if !isone(buffer)
        operate!(div, x, buffer)
        operate!(div, y, buffer)
    end
    return
end

_buffered_simplify(buffer, x::Rational) = _buffered_divgcd(buffer, x.num, x.den)

# + and -

function promote_operation(
    ::Union{typeof(+),typeof(-)},
    ::Type{Rational{S}},
    ::Type{Rational{T}},
) where {S,T}
    return Rational{promote_sum_mul(S, T)}
end

function buffer_for(
    ::Union{typeof(+),typeof(-)},
    ::Type{Rational{S}},
    ::Type{Rational{T}},
) where {S,T}
    U = promote_operation(gcd, S, T)
    return zero(U), zero(U), zero(U)
end

function buffered_operate_to!(
    buffer::Tuple,
    output::Rational,
    op::Union{typeof(-),typeof(+)},
    x::Rational,
    y::Rational,
)
    _buffered_divgcd(buffer[1], x.den, buffer[2], y.den, buffer[3])
    # TODO: Use `checked_mul` and `checked_sub` like in Base
    operate_to!(output.num, *, x.num, buffer[3])
    buffered_operate!(
        buffer[2],
        add_sub_mul_op(op),
        output.num,
        y.num,
        buffer[2],
    )
    operate_to!(output.den, *, x.den, buffer[3])
    _buffered_simplify(buffer[1], output)
    return output
end

function operate_to!(
    output::Rational,
    op::Union{typeof(+),typeof(-)},
    x::Rational,
    y::Rational,
)
    return buffered_operate_to!(
        buffer_for(op, typeof(x), typeof(y)),
        output,
        op,
        x,
        y,
    )
end

function operate_to!(output::Rational, ::typeof(+), x::Rational, y::Rational)
    xd, yd = Base.divgcd(promote(x.den, y.den)...)
    # TODO: Use `checked_mul` and `checked_add` like in Base
    operate_to!(output.num, *, x.num, yd)
    operate!(add_mul, output.num, y.num, xd)
    operate_to!(output.den, *, x.den, yd)
    # Reuse `xd` as it is a local copy created by this method
    _buffered_simplify(xd, output)
    return output
end

function operate_to!(out::Q, ::typeof(+), x::Q) where {Q<:Rational}
    return operate_to!(out, copy, x)
end

operate!(::typeof(+), x::Rational) = x

function operate_to!(out::Q, ::typeof(-), x::Q) where {Q<:Rational}
    operate_to!(out.num, -, x.num)
    operate_to!(out.den, copy, x.den)
    return out
end

function operate!(::typeof(-), x::Rational)
    operate!(-, x.num)
    return x
end

# abs

promote_operation(::typeof(abs), ::Type{Q}) where {Q<:Rational} = Q

function operate!(::typeof(abs), x::Rational)
    operate!(abs, x.num)
    return x
end

function operate_to!(o::Q, ::typeof(abs), x::Q) where {Q<:Rational}
    operate_to!(o.num, abs, x.num)
    operate_to!(o.den, copy, x.den)
    return o
end

# *

function promote_operation(
    ::typeof(*),
    ::Type{Rational{S}},
    ::Type{Rational{T}},
) where {S,T}
    return Rational{promote_operation(*, S, T)}
end

function buffer_for(
    ::typeof(*),
    ::Type{Rational{S}},
    ::Type{Rational{T}},
) where {S,T}
    U = promote_operation(gcd, S, T)
    return zero(Rational{U}), zero(Rational{U}), zero(U)
end

function buffered_operate_to!(
    buffer::Tuple,
    output::Rational,
    ::typeof(*),
    x::Rational,
    y::Rational,
)
    # Cannot use `output.num` and `output.den` as buffer as `output` might be an alias for `x`
    _buffered_divgcd(buffer[3], x.num, buffer[1].num, y.den, buffer[2].den)
    _buffered_divgcd(buffer[3], x.den, buffer[1].den, y.num, buffer[2].num)
    operate_to!(output.num, *, buffer[1].num, buffer[2].num)
    operate_to!(output.den, *, buffer[1].den, buffer[2].den)
    return output
end

function operate_to!(output::Rational, ::typeof(*), x::Rational, y::Rational)
    return buffered_operate_to!(
        buffer_for(*, typeof(x), typeof(y)),
        output,
        *,
        x,
        y,
    )
end

function operate_to!(out::Q, ::typeof(*), x::Q) where {Q<:Rational}
    return operate_to!(out, copy, x)
end

operate!(::typeof(*), x::Rational) = x

# gcd

function promote_operation(
    ::Union{typeof(gcd),typeof(lcm)},
    ::Type{Rational{S}},
    ::Type{Rational{T}},
) where {S,T}
    return Rational{promote_operation(gcd, S, T)}
end

function operate_to!(output::Rational, ::typeof(gcd), a::Rational, b::Rational)
    operate_to!(output.num, gcd, a.num, b.num)
    operate_to!(output.den, lcm, a.den, b.den)
    return output
end

function operate_to!(output::Rational, ::typeof(lcm), a::Rational, b::Rational)
    operate_to!(output.num, lcm, a.num, b.num)
    operate_to!(output.den, gcd, a.den, b.den)
    return output
end

function operate_to!(
    output::Rational,
    op::Union{typeof(+),typeof(-),typeof(*)},
    a::Rational,
    b::Rational,
    c::Vararg{Rational,N},
) where {N}
    operate_to!(output, op, a, b)
    return operate!(op, output, c...)
end

function operate!(op::Function, x::Rational, args::Vararg{Any,N}) where {N}
    return operate_to!(x, op, x, args...)
end

# add_mul and sub_mul

# Buffer to hold the product
function buffer_for(op::AddSubMul, args::Vararg{Type{<:Rational},N}) where {N}
    U = promote_operation(*, args...)
    return buffer_for(*, Base.tail(args)...),
    zero(U),
    buffer_for(add_sub_op(op), args[1], U)
end

function operate_to!(
    output::Rational,
    op::AddSubMul,
    x::Rational,
    y::Rational,
    z::Rational,
    args::Vararg{Rational,N},
) where {N}
    buffer = buffer_for(op, typeof(x), typeof(y), typeof(z), typeof.(args)...)
    return buffered_operate_to!(buffer, output, op, x, y, z, args...)
end

function buffered_operate_to!(
    buffer::Tuple,
    output::Rational,
    op::AddSubMul,
    a::Rational,
    x::Rational,
    y::Rational,
    args::Vararg{Rational,N},
) where {N}
    buffered_operate_to!(buffer[1], buffer[2], *, x, y, args...)
    return buffered_operate_to!(buffer[3], output, add_sub_op(op), a, buffer[2])
end

function buffered_operate!(
    buffer,
    op::AddSubMul,
    x::Rational,
    args::Vararg{Any,N},
) where {N}
    return buffered_operate_to!(buffer, x, op, x, args...)
end

function operate_to!(
    output::Rational,
    op::Union{typeof(+),typeof(-),typeof(*)},
    args::Vararg{Scaling,N},
) where {N}
    return operate_to!(output, op, _scaling_to.(typeof(output), args)...)
end

function operate_to!(
    output::Rational,
    op::AddSubMul,
    x::Scaling,
    y::Scaling,
    z::Scaling,
    args::Vararg{Scaling,N},
) where {N}
    return operate_to!(
        output,
        op,
        _scaling_to(typeof(output), x),
        _scaling_to(typeof(output), y),
        _scaling_to(typeof(output), z),
        _scaling_to.(typeof(output), args)...,
    )
end

# Called for instance if `args` is `(v', v)` for a vector `v`.
function operate_to!(
    output::Rational,
    op::AddSubMul,
    x,
    y,
    z,
    args::Vararg{Any,N},
) where {N}
    return operate_to!(output, add_sub_op(op), x, *(y, z, args...))
end
