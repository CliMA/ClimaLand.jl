@inline function Base.view(
  A::AbstractPtrArray{T,N},
  i::Vararg{Union{Integer,AbstractRange,Colon},N}
) where {T,N}
  j = Base.to_indices(A, i)
  @boundscheck checkbounds(A, j...)
  PtrArray(SubArray(A, j))
end
@inline function Base.view(
  A::AbstractPtrArray{T,N},
  i::AbstractUnitRange
) where {T,N}
  V = vec(A)
  @boundscheck checkbounds(V, i)
  view(V, i)
end
@inline function Base.view(
  A::AbstractPtrArray{T,1},
  i::AbstractUnitRange
) where {T}
  @boundscheck checkbounds(A, i)
  sx = only(static_strides(A))
  if sx === static(1)
    p = pointer(A) + (first(i) - first(offsets(A))) * sizeof(T)
    PtrArray(p, (static_length(i),), (nothing,))
  else
    p = pointer(A) + (first(i) - first(offsets(A))) * sizeof(T) * sx
    PtrArray(p, (static_length(i),), (sx,))
  end
end

@inline function _view(
  B::BitPtrArray{N},
  i::Vararg{Union{Integer,AbstractRange,Colon},M}
) where {N,M}
  A = SubArray(B, Base.to_indices(B, i))
  p = _offset_ptr(stridedpointer(B), i)
  sz = static_size(A)
  sx = _sparse_strides(dense_dims(A), strides(A))
  R = map(Int, stride_rank(A))
  PtrArray(p, sz, sx, offsets(A), _compact_rank(Val(R)))
end
@inline function Base.view(
  B::BitPtrArray{N},
  i::Vararg{Union{Integer,AbstractRange,Colon},M}
) where {N,M}
  _view(B, i...)
end
@inline function Base.view(
  B::BitPtrArray{N},
  i::Vararg{Union{Integer,AbstractRange,Colon},N}
) where {N}
  _view(B, i...)
end

@inline function _bview_unitrange(A, i)
  sx = stride(A, static(1))
  @assert sx === static(1)
  off = (first(i) - first(offsets(A)))
  @assert off % 8 == 0
  p = pointer(A) + (off >>> 3)
  PtrArray(p, (static_length(i),), (sx,))
end
@inline function Base.view(A::BitPtrArray{N}, i::AbstractUnitRange) where {N}
  @assert static_size(A, static(1)) == static_strides(A, static(2))
  _bview_unitrange(A, i)
end
@inline function Base.view(A::BitPtrArray{1}, i::AbstractUnitRange)
  _bview_unitrange(A, i)
end
@inline Base.view(A::BitPtrArray{1}, ::Colon) = A
@inline Base.view(A::BitPtrArray{N}, ::Colon) where {N} = A

@inline function zview(
  A::AbstractPtrArray{T,N,R,S,X,O,P},
  i::Vararg{Union{Integer,AbstractRange,Colon},M}
) where {T,N,R,S,X,O,P,M}
  zero_offsets(view(A, i...))
end

@inline Base.view(A::AbstractPtrArray, ::Colon) = vec(A)
@inline zview(A::AbstractPtrArray, ::Colon) = vec(A)

@inline Base.view(A::AbstractPtrArray{<:Any,N}, ::Vararg{Colon,N}) where {N} = A
@inline zview(A::AbstractPtrArray{<:Any,N}, ::Vararg{Colon,N}) where {N} = A

@inline Base.view(A::AbstractPtrVector, ::Colon) = A
@inline zview(A::AbstractPtrVector, ::Colon) = A

"""
    rank_to_sortperm(::NTuple{N,Int}) -> NTuple{N,Int}

Returns the `sortperm` of the stride ranks.
"""
function rank_to_sortperm(R::NTuple{N,Int}) where {N}
  sp = ntuple(zero, Val{N}())
  r = ntuple(n -> sum(R[n] .≥ R), Val{N}())
  @inbounds for n = 1:N
    sp = Base.setindex(sp, n, r[n])
  end
  sp
end
rank_to_sortperm(R) = sortperm(R)
_one_to_nothing(x::I) where {I} = StrideReset(x)
_one_to_nothing(::One) = nothing
Base.@propagate_inbounds function square_view(A::PtrMatrix{T,R,S,X,O}, i::I) where {T,R,S,X,O,I}
  # sizes = size(A)
  # @boundscheck i <= min(sizes[1], sizes[2]) || throw(BoundsError(A, (i, i)))
  x = map(_one_to_nothing, static_strides(A))
  SquarePtrMatrix{T,R,I,typeof(x),O}(pointer(A), i, x, offsets(A))
end
# Base.@propagate_inbounds function square_view(A::AbstractMatrix, i)
#   StrideArray(square_view(PtrArray(A), i), preserve_buffer(A))
# end
Base.@propagate_inbounds function square_view(A::AbstractMatrix, i)
  @view(A[begin:begin-1+i, begin:begin-1+i])
end
