module StackViews

export StackView

using OffsetArrays

const SlicesType = Union{Tuple, AbstractArray{<:AbstractArray}}

"""
    StackView(slices...; [dims])
    StackView(slices, [dims])
    StackView{T}(args...; kwargs...)

Stack/concatenate a list of arrays `slices` along dimension `dims` without copying the data.

If not specified, `dims` is defined as `ndims(first(slices))+1`, i.e., a new dimension in the tail.
This works better for normal Julia arrays when their memory layout is column-major.

# Example

`StackView` works very similar to `cat` and its `vcat`/`hcat`/`hvcat` variants.

```jldoctest; setup=:(using StackViews)
julia> A = reshape(collect(1:6), 2, 3);

julia> B = reshape(collect(7:12), 2, 3);

julia> StackView(A, B, dims=3) # mostly equivalent to `cat(A, B, dims=3)`
2×3×2 StackView{$Int, 3, 3, $(Tuple{Matrix{Int}, Matrix{Int}})}:
[:, :, 1] =
 1  3  5
 2  4  6

[:, :, 2] =
 7   9  11
 8  10  12
```

The difference comes to when `dim <= ndims(first(slices))`, which is `dim <= 2` for this example.
`StackView` always creats new dimension while `cat`s do not.

```jldoctest; setup=:(using StackViews; A = reshape(collect(1:6), 2, 3); B = reshape(collect(7:12), 2, 3);)
julia> StackView(A, B, dims=1) # `cat(A, B, dims=1)` outputs 4×3 Matrix
2×2×3 StackView{$Int, 3, 1, $(Tuple{Matrix{Int}, Matrix{Int}})}:
[:, :, 1] =
 1  2
 7  8

[:, :, 2] =
 3   4
 9  10

[:, :, 3] =
  5   6
 11  12

julia> StackView(A, B, dims=2) # `cat(A, B, dims=2)` outputs 2×6 Matrix
2×2×3 StackView{$Int, 3, 2, $(Tuple{Matrix{Int}, Matrix{Int}})}:
[:, :, 1] =
 1  7
 2  8

[:, :, 2] =
 3   9
 4  10

[:, :, 3] =
 5  11
 6  12
```

!!! tip
    For type-stability, you can pass `Val` as dim, e.g., `StackView(slices, Val(3))`.
"""
struct StackView{T, N, D, A} <: AbstractArray{T, N}
    slices::A
    function StackView{T, N, D}(slices) where {T, N, D}
        all(A->axes(A) == axes(first(slices)), slices) || throw(ArgumentError("all slices should be of the same axes."))
        new{T, N, D, typeof(slices)}(slices)
    end
    StackView{T, N, 0}(::A) where {T, N, A} = throw(ArgumentError("`dims=Val(0)` is not supported, do you mean `dims=Val(1)`?"))
end
StackView(slices::AbstractArray...; dims=Val(_default_dims(slices))) = StackView(slices, dims)
StackView{T}(slices::AbstractArray...; dims=Val(_default_dims(slices))) where T = StackView{T}(slices, dims)
StackView(slices, dims::Int) = StackView(slices, Val(dims)) # type-unstable

function StackView(slices::SlicesType, dims::Val = Val(_default_dims(slices)))
   return StackView{_default_eltype(slices)}(slices, dims)
end
function StackView{T}(slices::SlicesType, dims::Val = Val(_default_dims(slices))) where T
    N = _max(dims, Val(_default_dims(slices)))
    # unify all the axes to 1-based ranges
    slices = map(OffsetArrays.no_offset_view, slices)
    return StackView{T, _value(N), _value(dims)}(slices)
end

@inline _default_dims(slices) = ndims(first(slices)) + 1
@inline function _default_eltype(slices)
    T = mapreduce(eltype, promote_type, slices)
    _isconcretetype(T) || throw(ArgumentError("Input arrays should be homogenous."))
    return T
end
function _isconcretetype(T)
    # We relax the restriction and allow `Union`
    # This is particularily useful for arrays with `missing` and `nothing`
    isconcretetype(T) && return true
    isa(T, Union) && return isconcretetype(T.a) && _isconcretetype(T.b)
    return false
end

function Base.size(A::StackView{T,N,D}) where {T,N,D}
    frame_size = size(first(A.slices))
    prev, post = Base.IteratorsMD.split(frame_size, Val(D-1))
    return (_append_tuple(prev, Val(D-1))..., length(A.slices), post...)
end

function Base.axes(A::StackView{T,N,D}) where {T,N,D}
    frame_axes = axes(first(A.slices))
    prev, post = Base.IteratorsMD.split(frame_axes, Val(D-1))

    # use homogenous range to make _append_tuple happy
    fill_range = _convert(eltype(prev), Base.OneTo(1))
    return (_append_tuple(prev, Val(D-1), fill_range)...,
            Base.OneTo(length(A.slices)),
            post...)
end

@inline function Base.getindex(A::StackView{T,N,D}, inds::Vararg{Int,N}) where {T,N,D}
    @boundscheck checkbounds(A, inds...)
    prev, post = Base.IteratorsMD.split(inds, Val(D-1))
    idx, post = first(post), Base.tail(post)
    return @inbounds A.slices[idx][prev..., post...]
end

@inline function Base.setindex!(A::StackView{T,N,D}, x, inds::Vararg{Int, N}) where {T,N,D}
    @boundscheck checkbounds(A, inds...)
    prev, post = Base.IteratorsMD.split(inds, Val(D-1))
    idx, post = first(post), Base.tail(post)
    @inbounds A.slices[idx][prev..., post...] = x
end

# utils
_convert(::Type{Base.Bottom}, idx)=idx
_convert(T::Type, idx)=convert(T, idx)

# For type stability
@inline _max(::Val{x}, ::Val{y}) where {x, y} = Val(max(x, y))
@inline _value(::Val{N}) where N = N
@inline _append_tuple(t::NTuple{N1}, ::Val{N1}, x=1) where N1 = t
@inline _append_tuple(t::NTuple{N1}, ::Val{N2}, x=1) where {N1, N2} = _append_tuple((t..., x), Val(N2), x)
end
