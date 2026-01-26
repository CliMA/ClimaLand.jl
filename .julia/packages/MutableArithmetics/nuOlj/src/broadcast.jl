# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

# This file provides a mutable broadcast interface for types that support
# mutation.

function _broadcasted_type(
    ::Broadcast.DefaultArrayStyle{N},
    ::Base.HasShape{N},
    ::Type{Eltype},
) where {N,Eltype}
    return Array{Eltype,N}
end

function _broadcasted_type(
    ::Broadcast.DefaultArrayStyle{N},
    ::Base.HasShape{N},
    ::Type{Bool},
) where {N}
    return BitArray{N}
end

"""
This method is a generic fallback for array types that are not
`DefaultArrayStyle`. Because we can't tell the container from a generic
broadcast style, we fallback to `Any`, which is always a valid super type (just
not a helpful one).

In MutableArithmetics, `_broadcasted_type` appears only in `promote_broadcast`,
which itself appears only in `broadcast_mutability`, and so types hitting this
method will fallback to the `IsNotMutable()` branch, which is the expected
outcome.
"""
_broadcasted_type(::Broadcast.BroadcastStyle, ::Base.HasShape, ::Type) = Any

# Same as `Base.Broadcast._combine_styles` but with types as argument.
_combine_styles() = Broadcast.DefaultArrayStyle{0}()

_combine_styles(c::Type) = Broadcast.result_style(Broadcast.BroadcastStyle(c))

function _combine_styles(c1::Type, c2::Type)
    return Broadcast.result_style(_combine_styles(c1), _combine_styles(c2))
end

@inline function _combine_styles(
    c1::Type,
    c2::Type,
    cs::Vararg{Type,N},
) where {N}
    return Broadcast.result_style(
        _combine_styles(c1),
        _combine_styles(c2, cs...),
    )
end

_combine_shapes(s) = s

function _combine_2_shapes(::Base.HasShape{N}, ::Base.HasShape{M}) where {N,M}
    return Base.HasShape{max(N, M)}()
end

function _combine_shapes(s1, s2, args::Vararg{Any,N}) where {N}
    return _combine_shapes(_combine_2_shapes(s1, s2), args...)
end

_shape(T) = Base.HasShape{ndims(T)}()

function _combine_sizes(args::Vararg{Any,N}) where {N}
    return _combine_shapes(_shape.(args)...)
end

function promote_broadcast(op::F, args::Vararg{Any,N}) where {F<:Function,N}
    # TODO: we could use `promote_operation` instead as `combine_eltypes` uses
    # `return_type` hence it may return a non-concrete type and we do not handle
    # that case.
    T = Base.Broadcast.combine_eltypes(op, args)
    return _broadcasted_type(
        _combine_styles(args...),
        _combine_sizes(args...),
        T,
    )
end

"""
    broadcast_mutability(T::Type, ::typeof(op), args::Type...)::MutableTrait

Return `IsMutable` to indicate an object of type `T` can be modified to be
equal to `broadcast(op, args...)`.
"""
function broadcast_mutability(T::Type, op, args::Vararg{Type,N}) where {N}
    if mutability(T) isa IsMutable && promote_broadcast(op, args...) == T
        return IsMutable()
    else
        return IsNotMutable()
    end
end

function broadcast_mutability(x, op, args::Vararg{Any,N}) where {N}
    return broadcast_mutability(typeof(x), op, typeof.(args)...)
end

# Some AbstractArray, like JuMP.Containers.SparseAxisArray, do not support
# Base.size. In such cases, we default to returning `false`, since we cannot
# safely decide whether a broadcast can be stored in `x` unless we know the
# sizes of the two entries.
function _try_size(x::AbstractArray)
    try
        return size(x)
    catch
        return missing
    end
end
_try_size(x::Array) = size(x)

_checked_size(x_size::Any, y::AbstractArray) = x_size == _try_size(y)
_checked_size(::Any, ::Any) = true
_checked_size(::Any, ::Tuple{}) = true
function _checked_size(x_size::Any, y::Tuple)
    return _checked_size(x_size, y[1]) && _checked_size(x_size, Base.tail(y))
end
_checked_size(::Missing, ::Tuple) = false
_checked_size(::Missing, ::Tuple{}) = false

# This method is a slightly tricky one:
#
# If the elements in the broadcast are different sized arrays, weird things can
# happen during broadcasting since we'll either need to return a different size
# to `x`, or multiple copies of an argument will be used for different parts of
# `x`. To simplify, let's just return `IsNotMutable` if the sizes are different,
# which will be slower but correct.
function broadcast_mutability(
    x::AbstractArray,
    op,
    args::Vararg{Any,N},
) where {N}
    if !_checked_size(_try_size(x), args)::Bool
        return IsNotMutable()
    end
    return broadcast_mutability(typeof(x), op, typeof.(args)...)
end

broadcast_mutability(::Type) = IsNotMutable()

"""
    broadcast!(op::Function, args...)

Modify the value of `args[1]` to be equal to the value of
`broadcast(op, args...)`.

This method can only be called if `mutability(args[1], op, args...)` returns
[`IsMutable`](@ref).
"""
function broadcast! end

function mutable_broadcasted(broadcasted::Broadcast.Broadcasted{S}) where {S}
    function f(args::Vararg{Any,N}) where {N}
        return operate!!(broadcasted.f, args...)
    end
    return Broadcast.Broadcasted{S}(f, broadcasted.args, broadcasted.axes)
end

# If A is `Symmetric`, we cannot do a normal broadcast because we might modify
# the same entry twice. See https://github.com/jump-dev/JuMP.jl/issues/2102
function broadcast!(
    op::F,
    A::Union{Array,SparseArrays.AbstractSparseArray},
    args::Vararg{Any,N},
) where {F<:Function,N}
    bc = Broadcast.broadcasted(op, A, args...)
    instantiated = Broadcast.instantiate(bc)
    return copyto!(A, mutable_broadcasted(instantiated))
end

_any_uniform_scaling() = false

function _any_uniform_scaling(
    ::LinearAlgebra.UniformScaling,
    args::Vararg{Any,N},
) where {N}
    return true
end

function _any_uniform_scaling(::Any, args::Vararg{Any,N}) where {N}
    return _any_uniform_scaling(args...)
end

"""
    broadcast!!(op::Function, args...)

Returns the value of `broadcast(op, args...)`, possibly modifying `args[1]`.
"""
function broadcast!!(op::F, args::Vararg{Any,N}) where {F<:Function,N}
    # `any(x -> x isa LinearAlgebra.UniformScaling, args)` produces
    # `(1 allocation: 32 bytes)` on Julia v1.6.1 so we use
    # `_any_uniform_scaling` instead.
    if _any_uniform_scaling(args...)
        return _broadcast_with_uniform_scaling!(op, args...)
    else
        return broadcast_fallback!(
            broadcast_mutability(args[1], op, args...),
            op,
            args...,
        )
    end
end

function _broadcast_with_uniform_scaling!(
    op::F,
    args::Vararg{Any,N},
) where {F<:Function,N}
    return op(args...)
end

function broadcast_fallback!(
    ::IsNotMutable,
    op::F,
    args::Vararg{Any,N},
) where {F<:Function,N}
    return broadcast(op, args...)
end

function broadcast_fallback!(
    ::IsMutable,
    op::F,
    args::Vararg{Any,N},
) where {F<:Function,N}
    return broadcast!(op, args...)
end
