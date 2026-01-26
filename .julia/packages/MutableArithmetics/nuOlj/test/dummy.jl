# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

using LinearAlgebra
import MutableArithmetics as MA

# It does not support operation with floats on purpose to test that
# MutableArithmetics does not convert to float when it shouldn't.
struct DummyBigInt <: MA.AbstractMutable
    data::BigInt
end
DummyBigInt(J::UniformScaling) = DummyBigInt(J.λ)

# Broadcast
Base.broadcastable(x::DummyBigInt) = Ref(x)
# The version with `DummyBigInt` without `Type` is needed in LinearAlgebra for
# Julia v1.6+.
Base.ndims(::Union{Type{DummyBigInt},DummyBigInt}) = 0

function Base.promote_rule(
    ::Type{DummyBigInt},
    ::Type{<:Union{Integer,UniformScaling{<:Integer}}},
)
    return DummyBigInt
end
# `copy` on BigInt returns the same instance anyway
Base.copy(x::DummyBigInt) = x
MA.mutable_copy(x::DummyBigInt) = DummyBigInt(MA.mutable_copy(x.data))
LinearAlgebra.symmetric_type(::Type{DummyBigInt}) = DummyBigInt
LinearAlgebra.symmetric(x::DummyBigInt, ::Symbol) = x
LinearAlgebra.issymmetric(::DummyBigInt) = true
LinearAlgebra.hermitian_type(::Type{DummyBigInt}) = DummyBigInt
LinearAlgebra.hermitian(x::DummyBigInt, ::Symbol) = x
LinearAlgebra.dot(x::DummyBigInt, y::DummyBigInt) = x * y
function LinearAlgebra.dot(
    x::DummyBigInt,
    y::Union{Integer,UniformScaling{<:Integer}},
)
    return x * y
end
function LinearAlgebra.dot(
    x::Union{Integer,UniformScaling{<:Integer}},
    y::DummyBigInt,
)
    return x * y
end
LinearAlgebra.transpose(x::DummyBigInt) = x
LinearAlgebra.adjoint(x::DummyBigInt) = x
MA.mutability(::Type{DummyBigInt}) = MA.IsMutable()
MA.promote_operation(::typeof(zero), ::Type{DummyBigInt}) = DummyBigInt
MA.promote_operation(::typeof(one), ::Type{DummyBigInt}) = DummyBigInt
function MA.promote_operation(
    ::typeof(+),
    ::Type{DummyBigInt},
    ::Type{DummyBigInt},
)
    return DummyBigInt
end
_data(x) = x
_data(x::DummyBigInt) = x.data
MA.scaling(x::DummyBigInt) = x

function MA.operate_to!(
    x::DummyBigInt,
    op::Function,
    args::Union{MA.Scaling,DummyBigInt}...,
)
    return DummyBigInt(MA.operate_to!(x.data, op, _data.(args)...))
end

function MA.buffer_for(
    ::MA.AddSubMul,
    args::Vararg{Type{DummyBigInt},N},
) where {N}
    return DummyBigInt(BigInt())
end
_undummy(x) = x
_undummy(x::DummyBigInt) = x.data
function MA.buffered_operate_to!(
    buffer::DummyBigInt,
    output::DummyBigInt,
    op::Function,
    args::Vararg{Union{MA.Scaling,DummyBigInt},N},
) where {N}
    MA.buffered_operate_to!(
        buffer.data,
        output.data,
        output,
        op,
        _undummy.(args)...,
    )
    return output
end
function MA.buffered_operate!(
    buffer::DummyBigInt,
    op::Function,
    args::Vararg{Union{MA.Scaling,DummyBigInt},N},
) where {N}
    MA.buffered_operate!(buffer.data, op, _undummy.(args)...)
    return args[1]
end

# Called for instance if `args` is `(v', v)` for a vector `v`.
function MA.operate_to!(
    output::DummyBigInt,
    op::MA.AddSubMul,
    x::Union{MA.Scaling,DummyBigInt},
    y::Union{MA.Scaling,DummyBigInt},
    z::Union{MA.Scaling,DummyBigInt},
    args::Union{MA.Scaling,DummyBigInt}...,
)
    return MA.operate_to!(output, MA.add_sub_op(op), x, *(y, z, args...))
end
function MA.operate_to!(output::DummyBigInt, op::MA.AddSubMul, x, y, z, args...)
    return MA.operate_to!(output, MA.add_sub_op(op), x, *(y, z, args...))
end
function MA.operate!(
    op::Function,
    x::DummyBigInt,
    args::Vararg{Any,N},
) where {N}
    return MA.operate_to!(x, op, x, args...)
end

function MA.operate!(op::Union{typeof(zero),typeof(one)}, x::DummyBigInt)
    return DummyBigInt(MA.operate!(op, x.data))
end

function MA.promote_operation(
    ::typeof(*),
    ::Type{DummyBigInt},
    ::Type{DummyBigInt},
)
    return DummyBigInt
end
Base.convert(::Type{DummyBigInt}, x::Int) = DummyBigInt(x)
MA.isequal_canonical(x::DummyBigInt, y::DummyBigInt) = x.data == y.data
Base.iszero(x::DummyBigInt) = iszero(x.data)
Base.isone(x::DummyBigInt) = isone(x.data)
# We don't define == to tests that implementation of MA can pass the tests without defining ==.
# This is the case for MOI functions for instance.
# For th same reason, we only define `zero` and `one` for `Type{DummyBigInt}`, not for `DummyBigInt`.
Base.zero(::Type{DummyBigInt}) = DummyBigInt(zero(BigInt))
Base.one(::Type{DummyBigInt}) = DummyBigInt(one(BigInt))
Base.:+(x::DummyBigInt) = DummyBigInt(+x.data)
Base.:+(x::DummyBigInt, y::DummyBigInt) = DummyBigInt(x.data + y.data)
function Base.:+(x::DummyBigInt, y::Union{Integer,UniformScaling{<:Integer}})
    return DummyBigInt(x.data + y)
end
function Base.:+(x::Union{Integer,UniformScaling{<:Integer}}, y::DummyBigInt)
    return DummyBigInt(x + y.data)
end
Base.:-(x::DummyBigInt) = DummyBigInt(-x.data)
Base.:-(x::DummyBigInt, y::DummyBigInt) = DummyBigInt(x.data - y.data)
function Base.:-(x::Union{Integer,UniformScaling{<:Integer}}, y::DummyBigInt)
    return DummyBigInt(x - y.data)
end
function Base.:-(x::DummyBigInt, y::Union{Integer,UniformScaling{<:Integer}})
    return DummyBigInt(x.data - y)
end
Base.:*(x::DummyBigInt) = x
Base.:*(x::DummyBigInt, y::DummyBigInt) = DummyBigInt(x.data * y.data)
function Base.:*(x::Union{Integer,UniformScaling{<:Integer}}, y::DummyBigInt)
    return DummyBigInt(x * y.data)
end
function Base.:*(x::DummyBigInt, y::Union{Integer,UniformScaling{<:Integer}})
    return DummyBigInt(x.data * y)
end
function Base.:^(x::DummyBigInt, y::Union{Integer,UniformScaling{<:Integer}})
    return DummyBigInt(x.data^y)
end

Base.:*(x::Complex, y::DummyBigInt) = x * y.data
