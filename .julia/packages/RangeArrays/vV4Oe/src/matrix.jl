"""
    RangeMatrix(rs::AbstractVector{T}) where T <: AbstractRange

A RangeMatrix is a simple matrix representation of a vector of ranges, with
each range representing one column. Construct a RangeMatrix with a vector of
ranges; the ranges must all have the same length.

Note that it is only efficient when all the ranges are of the same type and in
a concretely typed Vector.
"""
struct RangeMatrix{T,A} <: AbstractArray{T,2}
    rs::A # A <: AbstractVector{_<:Range{T}}
    dims::Tuple{Int,Int}
end
RangeMatrix(rs::AbstractRange...) = RangeMatrix(collect(rs)) # TODO: use tuple storage?
function RangeMatrix(rs::AbstractVector{T}) where T <: AbstractRange
    n = length(rs)
    n == 0 && return RangeMatrix{T}(rs, (0, 0))
    m = length(rs[1])
    for j=2:n
        m == length(rs[j]) || throw(ArgumentError("all ranges must have the same length; expected $m, got $(length(rs[j]))"))
    end
    RangeMatrix{eltype(T), typeof(rs)}(rs, (m, n))
end

Base.size(R::RangeMatrix) = R.dims
Base.IndexStyle(::Type{<:RangeMatrix}) = IndexCartesian()

# Scalar indexing
@inline function Base.getindex(R::RangeMatrix, i::Int, j::Int)
    @boundscheck checkbounds(R, i, j);
    @inbounds return R.rs[j][i]
end

# For non-scalar indexing, only specialize with inner Ranges and Colons to
# return Ranges or RangeMatrixes. For everything else, we can use the fallbacks.
@inline function Base.getindex(R::RangeMatrix, I::Union{AbstractRange, Colon}, J)
    @boundscheck checkbounds(R, I, J)
    unsafe_getindex(R, I, J)
end
@inline unsafe_getindex(R::RangeMatrix, I::Union{AbstractRange, Colon}, j::Real) =
    @inbounds return R.rs[j][I]
@inline unsafe_getindex(R::RangeMatrix, I::Union{AbstractRange, Colon}, ::Colon) =
    @inbounds return RangeMatrix([R.rs[j][I] for j=1:length(R.rs)])
@inline unsafe_getindex(R::RangeMatrix, I::Union{AbstractRange, Colon}, J) =
    @inbounds return RangeMatrix([R.rs[j][I] for j in J])

# We can also optimize bounds checks to only look at each range's endpoints
function Base.checkindex(::Type{Bool}, inds::AbstractUnitRange, R::RangeMatrix)
    b = true
    @inbounds for r in R.rs
        b &= checkindex(Bool, inds, r[1])
        b &= checkindex(Bool, inds, r[end])
    end
    b
end
