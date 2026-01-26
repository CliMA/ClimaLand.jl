"""
    RepeatedRangeMatrix(r::AbstractRange{T}, at::AbstractVector{T}) where T

A RepeatedRange is like a RangeMatrix, except that it only stores one range and
a vector of offsets, at which the range repeats. For now, both the range and
vector of offsets must have the same element type.
"""
struct RepeatedRangeMatrix{T,R,A} <: AbstractMatrix{T}
    r::R # <: Range{T}
    at::A #<: AbstractVector{T}
end
RepeatedRangeMatrix(r::AbstractRange{T}, at::AbstractVector{T}) where T =
    RepeatedRangeMatrix{T, typeof(r), typeof(at)}(r, at)

Base.size(R::RepeatedRangeMatrix) = (length(R.r), length(R.at))
Base.IndexStyle(::Type{<:RepeatedRangeMatrix}) = IndexCartesian()

if VERSION < v"0.7.0-DEV.5126"
# This coupled iteration over the two fields is 10-20x faster than Cartesian iteration
# TODO: re-implement in the new iteration protocol
@inline function Base.start(R::RepeatedRangeMatrix)
    is = start(R.r)
    idone = done(R.r, is)
    js = start(R.at)
    jdone = done(R.at, js)
    return (idone | jdone) ? ((one(eltype(R.r)), is), (one(eltype(R.at)), js), true) :
                             (next(R.r, is), next(R.at, js), false)
end
@inline function Base.next(R::RepeatedRangeMatrix, state)
    (i, is), (j, js), _ = state
    val = i + j
    if done(R.r, is)
        if done(R.at, js)
            return (val, ((i, is), (j, js), true))
        end
        is = start(R.r)
        j, js = next(R.at, js)
    end
    return (val, (next(R.r, is), (j, js), false))
end
@inline Base.done(R::RepeatedRangeMatrix, state) = state[end]
end

# Scalar indexing
@inline function Base.getindex(R::RepeatedRangeMatrix, i::Int, j::Int)
    @boundscheck checkbounds(R, i, j)
    @inbounds return R.r[i] + R.at[j]
end

# For non-scalar indexing, only specialize with inner Ranges and Colons to
# return Ranges or RangeMatrixes. For everything else, we can use the fallbacks.
if VERSION < v"0.7.0-alpha"
    @inline function Base.getindex(R::RepeatedRangeMatrix, I::Union{AbstractRange, Colon}, j::Real)
        @boundscheck checkbounds(R, I, j)
        @inbounds return R.r[I] + R.at[j]
    end
else
    @inline function Base.getindex(R::RepeatedRangeMatrix, I::Union{AbstractRange, Colon}, j::Real)
        @boundscheck checkbounds(R, I, j)
        @inbounds return R.r[I] .+ R.at[j]
    end
end
@inline function Base.getindex(R::RepeatedRangeMatrix, I::Union{AbstractRange, Colon}, J)
    @boundscheck checkbounds(R, I, J)
    @inbounds return RepeatedRangeMatrix(R.r[I], R.at[J])
end

# We can also optimize bounds checks to only look at the range's endpoints
function Base.checkindex(::Type{Bool}, inds::AbstractUnitRange, R::RepeatedRangeMatrix)
    b = true
    @inbounds for a in R.at
        b &= checkindex(Bool, inds, R.r[1] + a)
        b &= checkindex(Bool, inds, R.r[end] + a)
    end
    b
end
