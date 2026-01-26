module StaticArrayInterfaceOffsetArraysExt

using StaticArrayInterface
using StaticArrayInterface.Static

if isdefined(Base, :get_extension) 
    using OffsetArrays
else 
    using ..OffsetArrays
end

relative_offsets(r::OffsetArrays.IdOffsetRange) = (getfield(r, :offset),)
relative_offsets(A::OffsetArrays.OffsetArray) = getfield(A, :offsets)
function relative_offsets(A::OffsetArrays.OffsetArray, ::StaticInt{dim}) where {dim}
    if dim > ndims(A)
        return static(0)
    else
        return getfield(relative_offsets(A), dim)
    end
end
function relative_offsets(A::OffsetArrays.OffsetArray, dim::Int)
    if dim > ndims(A)
        return 0
    else
        return getfield(relative_offsets(A), dim)
    end
end
StaticArrayInterface.parent_type(::Type{<:OffsetArrays.OffsetArray{T,N,A}}) where {T,N,A} = A
function _offset_axis_type(::Type{T}, dim::StaticInt{D}) where {T,D}
    OffsetArrays.IdOffsetRange{Int,StaticArrayInterface.axes_types(T, dim)}
end
function StaticArrayInterface.axes_types(::Type{T}) where {T<:OffsetArrays.OffsetArray}
    Static.eachop_tuple(
        _offset_axis_type,
        ntuple(static, StaticInt(ndims(T))),
        StaticArrayInterface.parent_type(T)
    )
end
StaticArrayInterface.static_strides(A::OffsetArray) = StaticArrayInterface.static_strides(parent(A))
function StaticArrayInterface.known_offsets(::Type{A}) where {A<:OffsetArrays.OffsetArray}
    ntuple(identity -> nothing, Val(ndims(A)))
end
function StaticArrayInterface.offsets(A::OffsetArrays.OffsetArray)
    map(+, StaticArrayInterface.offsets(parent(A)), relative_offsets(A))
end
@inline function StaticArrayInterface.offsets(A::OffsetArrays.OffsetArray, dim)
    d = StaticArrayInterface.to_dims(A, dim)
    StaticArrayInterface.offsets(parent(A), d) + relative_offsets(A, d)
end
@inline function StaticArrayInterface.static_axes(A::OffsetArrays.OffsetArray)
    map(OffsetArrays.IdOffsetRange, StaticArrayInterface.static_axes(parent(A)), relative_offsets(A))
end
@inline function StaticArrayInterface.static_axes(A::OffsetArrays.OffsetArray, dim)
    d = StaticArrayInterface.to_dims(A, dim)
    OffsetArrays.IdOffsetRange(StaticArrayInterface.static_axes(parent(A), d), relative_offsets(A, d))
end
function StaticArrayInterface.stride_rank(T::Type{<:OffsetArray})
  StaticArrayInterface.stride_rank(StaticArrayInterface.parent_type(T))
end
function StaticArrayInterface.dense_dims(T::Type{<:OffsetArray})
    StaticArrayInterface.dense_dims(StaticArrayInterface.parent_type(T))
end
function StaticArrayInterface.contiguous_axis(T::Type{<:OffsetArray})
  StaticArrayInterface.contiguous_axis(StaticArrayInterface.parent_type(T))
end

end # module
