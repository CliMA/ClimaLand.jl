@inline Base.vec(A::PtrArray{T,1}) where {T} = A
@inline function Base.vec(
  A::PtrArray{T,N,R,S,NTuple{N,Nothing}}
) where {T,N,R,S}
  PtrArray(pointer(A), (static_length(A),), (nothing,), (One(),), Val((1,)))
end
@inline Base.vec(A::StaticStrideArray{T,1}) where {T} = A
@inline Base.vec(A::StaticStrideArray{T,N}) where {T,N} =
  StrideArray(vec(PtrArray(A)), A)
@inline Base.vec(A::StrideArray) =
  StrideArray(vec(PtrArray(A)), preserve_buffer(A))

@inline function Base.reshape(
  A::PtrArray{T,N,R,S,NTuple{N,Nothing}},
  dims::Tuple{Vararg{Integer}}
) where {T,N,R,S}
  PtrArray(pointer(A), dims)
end

@inline function Base.reshape(
  A::PtrArray{T,N,R,S,NTuple{N,Nothing}},
  dims::Tuple{Vararg{Base.Integer}}
) where {T,N,R,S}
  PtrArray(pointer(A), dims)
end

@inline function Base.reshape(
  A::PtrArray{T,N,R,S,NTuple{N,Nothing}},
  dims::Tuple{Int,Vararg{Int}}
) where {T,N,R,S}
  PtrArray(pointer(A), dims)
end
@inline function Base.reshape(
  A::PtrArray{T,N,R,S,NTuple{N,Nothing}},
  dims::Tuple{Integer,Vararg{Integer}}
) where {T,N,R,S}
  PtrArray(pointer(A), dims)
end
@inline function Base.reshape(
  A::PtrArray{T,N,R,S,NTuple{N,Nothing}},
  dims::Tuple{Base.Integer,Vararg{Base.Integer}}
) where {T,N,R,S}
  PtrArray(pointer(A), dims)
end

@inline function Base.reshape(
  A::PtrArray{T,N,R,S,NTuple{N,Nothing}},
  dims::Tuple{Vararg{Int}}
) where {T,N,R,S}
  PtrArray(pointer(A), dims)
end

@inline function Base.reshape(
  A::PtrArray{T,N,R,S,NTuple{N,Nothing}},
  ::Tuple{}
) where {T,N,R,S}
  PtrArray(pointer(A), ())
end

@inline function Base.reshape(A::StrideArray, dims::Tuple{Int})
  StrideArray(PtrArray(pointer(A), dims), preserve_buffer(A))
end
@inline function Base.reshape(A::StrideArray, dims::Tuple{Vararg{Int}})
  StrideArray(PtrArray(pointer(A), dims), preserve_buffer(A))
end
@inline function Base.reshape(A::StrideArray, dims::Tuple{Int,Vararg{Int}})
  StrideArray(PtrArray(pointer(A), dims), preserve_buffer(A))
end
@inline function Base.reshape(A::StrideArray, dims::Tuple{Int,Int,Vararg{Int}})
  StrideArray(PtrArray(pointer(A), dims), preserve_buffer(A))
end
@inline function Base.reshape(A::StrideArray, dims::Tuple{Integer})
  StrideArray(PtrArray(pointer(A), dims), preserve_buffer(A))
end
@inline function Base.reshape(
  A::StrideArray,
  dims::Tuple{Integer,Integer,Vararg{Integer}}
)
  StrideArray(PtrArray(pointer(A), dims), preserve_buffer(A))
end
@inline function Base.reshape(A::StrideArray, dims::Tuple{Vararg{Base.Integer}})
  StrideArray(PtrArray(pointer(A), dims), preserve_buffer(A))
end
@inline function Base.reshape(
  A::StrideArray,
  dims::Tuple{Base.Integer,Vararg{Base.Integer}}
)
  StrideArray(PtrArray(pointer(A), dims), preserve_buffer(A))
end
