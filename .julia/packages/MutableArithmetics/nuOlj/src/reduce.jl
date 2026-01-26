# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

_check_same_length(::Any) = nothing

function _check_same_length(a, b, c::Vararg{Any,N}) where {N}
    if length(a) != length(b)
        throw(
            DimensionMismatch(
                "one array has length $(length(a)) which does not match the " *
                "length of the next one, $(length(b)).",
            ),
        )
    end
    return _check_same_length(b, c...)
end

reduce_op(op::AddSubMul) = add_sub_op(op)

reduce_op(::typeof(add_dot)) = +

neutral_element(::typeof(+), ::Type{T}) where {T} = zero(T)
neutral_element(::typeof(+), ::Type{T}) where {T<:AbstractArray} = Zero()

map_op(::AddSubMul) = *

map_op(::typeof(add_dot)) = LinearAlgebra.dot

function promote_map_reduce(op::Function, args::Vararg{Any,N}) where {N}
    return promote_operation(
        op,
        promote_operation(map_op(op), args...),
        args...,
    )
end

_concrete_eltype(x) = isempty(x) ? eltype(x) : typeof(first(x))

function fused_map_reduce(op::F, args::Vararg{Any,N}) where {F<:Function,N}
    _check_same_length(args...)
    T = promote_map_reduce(op, _concrete_eltype.(args)...)
    accumulator = neutral_element(reduce_op(op), T)
    buffer = buffer_for(op, T, eltype.(args)...)
    for I in zip(eachindex.(args)...)
        accumulator =
            buffered_operate!!(buffer, op, accumulator, getindex.(args, I)...)
    end
    return accumulator
end

function operate(
    ::typeof(sum),
    a::AbstractArray;
    init = zero(promote_operation(+, eltype(a), eltype(a))),
)
    return mapreduce(identity, add!!, a; init)
end
