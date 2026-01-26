module IndirectArrays

export IndirectArray

"""
    IndirectArray(index, values)

creates an array `A` where the values are looked up in the value table,
`values`, using the `index`.  Concretely, `A[i...] =
values[index[i...]]`.

`values` can be an `AbstractVector` (useful when the values in `index` are
contiguous) or an `AbstractDict` (useful when they are not).
When there are not very many such distinct values, "frozen" LittleDicts from
[OrderedCollections](https://github.com/JuliaCollections/OrderedCollections.jl)
are recommended for performance.
"""
struct IndirectArray{T,N,I,A<:AbstractArray{I,N},V<:Union{AbstractVector{T},AbstractDict{I,T}}} <: AbstractArray{T,N}
    index::A
    values::V

    @inline function IndirectArray{T,N,I,A,V}(index, values) where {T,N,I,A,V}
        # The typical logic for testing bounds and then using
        # @inbounds will not check whether index is inbounds for
        # values. So we had better check this on construction.
        if isa(values, AbstractVector)
            I <: Integer || error("with a vector `values`, the index must have integer eltype")
            @boundscheck checkbounds(values, index)
        else
            @boundscheck begin
                uindex = unique(index)
                all(idx -> haskey(values, idx), index) || Base.throw_boundserror(index, values)
            end
        end
        new{T,N,I,A,V}(index, values)
    end
end
const IndirectArrayVec{T,N,I,A,V<:AbstractVector} = IndirectArray{T,N,I,A,V}

Base.@propagate_inbounds IndirectArray(index::AbstractArray{<:Integer,N}, values::AbstractVector{T}) where {T,N} =
    IndirectArray{T,N,eltype(index),typeof(index),typeof(values)}(index, values)
Base.@propagate_inbounds IndirectArray(index::AbstractArray{I,N}, values::AbstractDict{I,T}) where {I,T,N} =
    IndirectArray{T,N,I,typeof(index),typeof(values)}(index, values)

function IndirectArray{T}(A::AbstractArray) where {T}
    values = unique(A)
    index = convert(Array{T}, indexin(A, values))
    return IndirectArray(index, values)
end
IndirectArray(A::AbstractArray) = IndirectArray{UInt8}(A)

Base.size(A::IndirectArray) = size(A.index)
Base.axes(A::IndirectArray) = axes(A.index)
Base.IndexStyle(::Type{IndirectArray{T,N,I,A,V}}) where {T,N,I,A,V} = IndexStyle(A)

function Base.copy(A::IndirectArray)
    @inbounds ret = IndirectArray(copy(A.index), copy(A.values))
    return ret
end

@inline function Base.getindex(A::IndirectArray, i::Int)
    @boundscheck checkbounds(A.index, i)
    @inbounds idx = A.index[i]
    @inbounds ret = A.values[idx]
    ret
end

@inline function Base.getindex(A::IndirectArray{T,N}, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A.index, I...)
    @inbounds idx = A.index[I...]
    @inbounds ret = A.values[idx]
    ret
end

@inline function Base.setindex!(A::IndirectArrayVec, x, i::Int)
    @boundscheck checkbounds(A.index, i)
    idx = findfirst(isequal(x), A.values)
    if idx === nothing
        push!(A.values, x)
        A.index[i] = length(A.values)
    else
        A.index[i] = idx
    end
    return A
end

@inline function Base.push!(A::IndirectArrayVec{T,1} where T, x)
    idx = findfirst(isequal(x), A.values)
    if idx === nothing
        push!(A.values, x)
        push!(A.index, length(A.values))
    else
        push!(A.index, idx)
    end
    return A
end

function Base.append!(A::IndirectArrayVec{T,1}, B::IndirectArray{T,1}) where T
    if A.values == B.values
        append!(A.index, B.index)
    else # pretty inefficient but let's get something going
        for b in B
            push!(A, b)
        end
    end
    return A
end

function Base.append!(A::IndirectArrayVec{<:Any,1}, B::AbstractVector)
    for b in B
        push!(A, b)
    end
    return A
end

end # module
