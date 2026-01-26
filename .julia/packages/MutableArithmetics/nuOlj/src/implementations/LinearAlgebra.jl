# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

# This file contains methods to implement the MutableArithmetics API for types
# in the LinearAlgebra stdlib.

mutability(::Type{<:Array}) = IsMutable()

mutable_copy(A::Array) = copy_if_mutable.(A)

# Sum

# By default, we assume the return value is an `Array` because having a
# different method for all combinations of cases would be cumbersome. A more
# specific method could be implemented for specific cases.
function promote_operation(
    op::Union{typeof(+),typeof(-)},
    ::Type{<:AbstractArray{S,N}},
    ::Type{<:AbstractArray{T,M}},
) where {S,T,N,M}
    # If `N != M`, we need the axes between `min(N,M)+1` and `max(N,M)` to be
    # `Base.OneTo(1)`. In any cases, the axes from `1` to `min(N,M)` must also
    # match.
    return Array{promote_operation(op, S, T),max(N, M)}
end

function promote_operation(
    op::Union{typeof(+),typeof(-)},
    ::Type{LinearAlgebra.UniformScaling{S}},
    ::Type{Matrix{T}},
) where {S,T}
    return Matrix{promote_operation(op, S, T)}
end

function promote_operation(
    op::Union{typeof(+),typeof(-)},
    ::Type{Matrix{T}},
    ::Type{LinearAlgebra.UniformScaling{S}},
) where {S,T}
    return Matrix{promote_operation(op, S, T)}
end

# Only `Scaling`
function operate!(
    op::Union{typeof(+),typeof(-)},
    A::Matrix,
    B::LinearAlgebra.UniformScaling,
)
    n = LinearAlgebra.checksquare(A)
    for i in 1:n
        A[i, i] = operate!!(op, A[i, i], B)
    end
    return A
end

function operate!(
    op::AddSubMul,
    A::Matrix,
    B::Scaling,
    C::Scaling,
    D::Vararg{Scaling,N},
) where {N}
    return operate!(add_sub_op(op), A, *(B, C, D...))
end

# TODO(odow): these are the only cases that appear in all of JuliaHub. They
# should become private.
mul_rhs(::typeof(+)) = add_mul
mul_rhs(::typeof(-)) = sub_mul

# We redirect the mutable `A + B` into `A .+ B`.
# To be consistent with Julia Base, we first call `promote_shape`
# which throws an error if the broadcasted dimension are not singleton
# and we check that the axes of `A` are indeed the axes of the array
# that would be returned in Julia Base (maybe we could relax this ?).
function _check_dims(A, B)
    if axes(A) != promote_shape(A, B)
        throw(
            DimensionMismatch(
                "Cannot sum or substract a matrix of axes `$(axes(B))` into" *
                " matrix of axes `$(axes(A))`, expected axes" *
                " `$(promote_shape(A, B))`.",
            ),
        )
    end
    return
end

function operate!(
    op::Union{typeof(+),typeof(-)},
    A::AbstractArray,
    B::AbstractArray,
)
    _check_dims(A, B)
    return broadcast!(op, A, B)
end

function operate_to!(
    output::AbstractArray,
    op::Union{typeof(+),typeof(-)},
    A::AbstractArray,
)
    if axes(output) != axes(A)
        throw(
            DimensionMismatch(
                "Cannot sum or substract a matrix of axes `$(axes(A))`" *
                " into a matrix of axes `$(axes(output))`, expected" *
                " axes `$(axes(A))`.",
            ),
        )
    end
    # We don't have `MA.broadcast_to!` as it would be exactly `Base.broadcast!`.
    return Base.broadcast!(op, output, A)
end

function operate_to!(
    output::AbstractArray,
    op::Union{typeof(+),typeof(-)},
    A::AbstractArray,
    B::AbstractArray,
)
    if axes(output) != promote_shape(A, B)
        throw(
            DimensionMismatch(
                "Cannot sum or substract matrices of axes `$(axes(A))` and" *
                " `$(axes(B))` into a matrix of axes `$(axes(output))`," *
                " expected axes `$(promote_shape(A, B))`.",
            ),
        )
    end
    # We don't have `MA.broadcast_to!` as it would be exactly `Base.broadcast!`.
    return Base.broadcast!(op, output, A, B)
end

# We call `scaling_to_number` as `UniformScaling` do not support broadcasting
function operate!(
    op::AddSubMul,
    A::AbstractArray,
    B::AbstractArray,
    α::Vararg{Scaling,M},
) where {M}
    _check_dims(A, B)
    return broadcast!(op, A, B, scaling_to_number.(α)...)
end

function operate!(
    op::AddSubMul,
    A::AbstractArray,
    α::Scaling,
    B::AbstractArray,
    β::Vararg{Scaling,M},
) where {M}
    _check_dims(A, B)
    return broadcast!(op, A, scaling_to_number(α), B, scaling_to_number.(β)...)
end

function operate!(
    op::AddSubMul,
    A::AbstractArray,
    α1::Scaling,
    α2::Scaling,
    B::AbstractArray,
    β::Vararg{Scaling,M},
) where {M}
    _check_dims(A, B)
    return broadcast!(
        op,
        A,
        scaling_to_number(α1),
        scaling_to_number(α2),
        B,
        scaling_to_number.(β)...,
    )
end

# Fallback, we may be able to be more efficient in more cases by adding more
# specialized methods.
function operate!(op::AddSubMul, A::AbstractArray, x, y)
    return operate!(op, A, x * y)
end

function operate!(
    op::AddSubMul,
    A::AbstractArray,
    x,
    y,
    args::Vararg{Any,N},
) where {N}
    @assert N > 0
    return operate!(op, A, x, *(y, args...))
end

# Product

function similar_array_type(
    ::Type{LinearAlgebra.Symmetric{T,MT}},
    ::Type{S},
) where {S,T,MT}
    return LinearAlgebra.Symmetric{S,similar_array_type(MT, S)}
end

function similar_array_type(
    ::Type{LinearAlgebra.Diagonal{T,VT}},
    ::Type{S},
) where {S,T,VT<:AbstractVector{T}}
    return LinearAlgebra.Diagonal{S,similar_array_type(VT, S)}
end

similar_array_type(::Type{<:AbstractVector}, ::Type{T}) where {T} = Vector{T}
similar_array_type(::Type{Array{T,N}}, ::Type{S}) where {S,T,N} = Array{S,N}

similar_array_type(::Type{BitArray{N}}, ::Type{S}) where {S,N} = Array{S,N}
similar_array_type(::Type{BitArray{N}}, ::Type{Bool}) where {N} = BitArray{N}

function similar_array_type(
    ::Type{<:SubArray{T,N,<:Array{T}}},
    ::Type{S},
) where {S,T,N}
    return Array{S,N}
end

function promote_operation(
    op::typeof(*),
    A::Type{<:AbstractArray{T}},
    ::Type{S},
) where {S,T}
    return similar_array_type(A, promote_operation(op, T, S))
end

function promote_operation(
    op::typeof(*),
    ::Type{S},
    A::Type{<:AbstractArray{T}},
) where {S,T}
    return similar_array_type(A, promote_operation(op, S, T))
end

# `{S}` and `{T}` are used to avoid ambiguity with above methods.
function promote_operation(
    ::typeof(*),
    A::Type{<:AbstractArray{S}},
    B::Type{<:AbstractArray{T}},
) where {S,T}
    return promote_array_mul(A, B)
end

function promote_sum_mul(T::Type, S::Type)
    U = promote_operation(*, T, S)
    return promote_operation(+, U, U)
end

function promote_array_mul(::Type{Matrix{S}}, ::Type{Vector{T}}) where {S,T}
    return Vector{promote_sum_mul(S, T)}
end

function promote_array_mul(
    ::Type{<:AbstractMatrix{S}},
    ::Type{<:AbstractMatrix{T}},
) where {S,T}
    return Matrix{promote_sum_mul(S, T)}
end

function promote_array_mul(
    ::Type{<:AbstractVector{S}},
    ::Type{<:LinearAlgebra.Adjoint{T,<:AbstractVector{T}}},
) where {S,T}
    return Matrix{promote_operation(*, S, T)}
end

function promote_array_mul(
    ::Type{<:AbstractVector{S}},
    ::Type{<:LinearAlgebra.Transpose{T,<:AbstractVector{T}}},
) where {S,T}
    return Matrix{promote_operation(*, S, T)}
end

function promote_array_mul(
    ::Type{<:AbstractMatrix{S}},
    ::Type{<:AbstractVector{T}},
) where {S,T}
    return Vector{promote_sum_mul(S, T)}
end

function _dim_check(C::AbstractVector, A::AbstractMatrix, B::AbstractVector)
    mB = length(B)
    mA, nA = size(A)
    if mB != nA
        throw(
            DimensionMismatch(
                "matrix A has dimensions ($mA,$nA), vector B has length $mB",
            ),
        )
    end
    if mA != length(C)
        throw(
            DimensionMismatch(
                "result C has length $(length(C)), needs length $mA",
            ),
        )
    end
    return
end

function _dim_check(C::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix)
    mB, nB = size(B)
    mA, nA = size(A)
    if mB != nA
        throw(
            DimensionMismatch(
                "matrix A has dimensions ($mA,$nA), matrix B has dimensions ($mB,$nB)",
            ),
        )
    end
    if size(C, 1) != mA || size(C, 2) != nB
        throw(
            DimensionMismatch(
                "result C has dimensions $(size(C)), needs ($mA,$nB)",
            ),
        )
    end
    return
end

function buffered_operate!(
    buffer,
    ::typeof(add_mul),
    C::Vector,
    A::AbstractMatrix,
    B::AbstractVector,
)
    _dim_check(C, A, B)
    for (ci, ai) in zip(axes(C, 1), axes(A, 1))
        for (aj, bj) in zip(axes(A, 2), axes(B, 1))
            C[ci] = buffered_operate!!(buffer, add_mul, C[ci], A[ai, aj], B[bj])
        end
    end
    return C
end

function buffered_operate!(
    buffer,
    ::typeof(add_mul),
    C::Matrix,
    A::AbstractMatrix,
    B::AbstractMatrix,
)
    _dim_check(C, A, B)
    for (ci, ai) in zip(axes(C, 1), axes(A, 1))
        for (cj, bj) in zip(axes(C, 2), axes(B, 2))
            for (aj, bi) in zip(axes(A, 2), axes(B, 1))
                C[ci, cj] = buffered_operate!!(
                    buffer,
                    add_mul,
                    C[ci, cj],
                    A[ai, aj],
                    B[bi, bj],
                )
            end
        end
    end
    return C
end

function buffer_for(
    ::typeof(add_mul),
    ::Type{<:VecOrMat{S}},
    ::Type{<:AbstractMatrix{T}},
    ::Type{<:AbstractVecOrMat{U}},
) where {S,T,U}
    return buffer_for(add_mul, S, T, U)
end

function operate!(
    ::typeof(add_mul),
    C::VecOrMat,
    A::AbstractMatrix,
    B::AbstractVecOrMat,
)
    buffer = buffer_for(add_mul, typeof(C), typeof(A), typeof(B))
    return buffered_operate!(buffer, add_mul, C, A, B)
end

function operate!(::typeof(zero), C::Union{Vector,Matrix})
    # C may contain undefined values so we cannot call `zero!`
    for i in eachindex(C)
        @inbounds C[i] = zero(eltype(C))
    end
    return
end

function operate_to!(
    C::AbstractArray,
    ::typeof(*),
    A::AbstractArray,
    B::AbstractArray,
)
    operate!(zero, C)
    return operate!(add_mul, C, A, B)
end

function undef_array(::Type{Array{T,N}}, axes::Vararg{Base.OneTo,N}) where {T,N}
    return Array{T,N}(undef, length.(axes))
end

# This method is for things like StaticArrays which return something other than
# Base.OneTo for their axes. It isn't typed because there can be a mix of axes
# in the call.
function undef_array(::Type{T}, axes...) where {T}
    return undef_array(T, convert.(Base.OneTo, axes)...)
end

# Does what `LinearAlgebra/src/matmul.jl` does for abstract matrices and
# vectors: estimate the resulting element type, allocate the resulting array but
# it redirects to `mul_to!` instead of `LinearAlgebra.mul!`.
function operate(
    ::typeof(*),
    A::AbstractMatrix{S},
    B::AbstractVector{T},
) where {T,S}
    # Only use the efficient in-place operate_to! if both arrays are
    # concrete. Bad things can happen if S or T is abstract and we pick the
    # wrong type for C.
    if !(isconcretetype(S) && isconcretetype(T))
        if S <: AbstractMutable || T <: AbstractMutable
            # We can't use A * B here, because this will StackOverflow via:
            #   *(A, B) -> mul(A, B) -> operate(*, A, B) -> *(A, B)
            # See MutableArithmetics.jl#336
            if axes(A, 2) != axes(B, 1)
                throw(DimensionMismatch())
            end
            return [sum(A[i, j] * B[j] for j in axes(B, 1)) for i in axes(A, 1)]
        end
        return A * B
    end
    C = undef_array(promote_array_mul(typeof(A), typeof(B)), axes(A, 1))
    return operate_to!(C, *, A, B)
end

function operate(
    ::typeof(*),
    A::AbstractMatrix{S},
    B::AbstractMatrix{T},
) where {T,S}
    # Only use the efficient in-place operate_to! if both arrays are
    # concrete. Bad things can happen if S or T is abstract and we pick the
    # wrong type for C.
    if !(isconcretetype(S) && isconcretetype(T))
        # It's safe to use A * B here, because there is no fallback for
        #   *(A, B) -> mul(A, B) -> operate(*, A, B)
        # See MutableArithmetics.jl#336
        return A * B
    end
    C = undef_array(
        promote_array_mul(typeof(A), typeof(B)),
        axes(A, 1),
        axes(B, 2),
    )
    return operate_to!(C, *, A, B)
end

const _TransposeOrAdjoint{T,MT} =
    Union{LinearAlgebra.Transpose{T,MT},LinearAlgebra.Adjoint{T,MT}}

function _mirror_transpose_or_adjoint(x, ::LinearAlgebra.Transpose)
    return LinearAlgebra.transpose(x)
end

function _mirror_transpose_or_adjoint(x, ::LinearAlgebra.Adjoint)
    return LinearAlgebra.adjoint(x)
end

function _mirror_transpose_or_adjoint(
    A::Type{<:AbstractArray{T}},
    ::Type{<:LinearAlgebra.Transpose},
) where {T}
    return LinearAlgebra.Transpose{T,A}
end

function _mirror_transpose_or_adjoint(
    A::Type{<:AbstractArray{T}},
    ::Type{<:LinearAlgebra.Adjoint},
) where {T}
    return LinearAlgebra.Adjoint{T,A}
end

function similar_array_type(
    TA::Type{<:_TransposeOrAdjoint{T,A}},
    ::Type{S},
) where {S,T,A}
    return _mirror_transpose_or_adjoint(similar_array_type(A, S), TA)
end

# dot product

function promote_array_mul(
    ::Type{<:_TransposeOrAdjoint{S,<:AbstractVector}},
    ::Type{<:AbstractVector{T}},
) where {S,T}
    return promote_sum_mul(S, T)
end

function promote_array_mul(
    A::Type{<:_TransposeOrAdjoint{S,V}},
    M::Type{<:AbstractMatrix{T}},
) where {S,T,V<:AbstractVector}
    B = promote_array_mul(_mirror_transpose_or_adjoint(M, A), V)
    return _mirror_transpose_or_adjoint(B, A)
end

function operate(
    ::typeof(*),
    x::LinearAlgebra.Adjoint{<:Any,<:AbstractVector},
    y::AbstractVector,
)
    return operate(LinearAlgebra.dot, parent(x), y)
end

function operate(
    ::typeof(*),
    x::_TransposeOrAdjoint{<:Any,<:AbstractVector},
    y::AbstractMatrix,
)
    return _mirror_transpose_or_adjoint(
        operate(*, _mirror_transpose_or_adjoint(y, x), parent(x)),
        x,
    )
end

function operate(
    ::typeof(*),
    x::_TransposeOrAdjoint{<:Any,<:AbstractVector},
    y::AbstractVector,
)
    return fused_map_reduce(add_mul, x, y)
end

function operate(
    ::typeof(LinearAlgebra.dot),
    x::AbstractArray,
    y::AbstractArray,
)
    return fused_map_reduce(add_dot, x, y)
end
