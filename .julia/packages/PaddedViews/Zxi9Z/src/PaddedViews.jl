module PaddedViews
using Base: OneTo, tail
using OffsetArrays
using OffsetArrays: no_offset_view
@static if !isdefined(Base, :IdentityUnitRange)
    const IdentityUnitRange = Base.Slice
else
    using Base: IdentityUnitRange
end

export PaddedView, paddedviews, sym_paddedviews

"""
    datapadded = PaddedView(fillvalue, data, padded_axes)
    datapadded = PaddedView(fillvalue, data, padded_axes, data_axes)
    datapadded = PaddedView(fillvalue, data, sz)
    datapadded = PaddedView(fillvalue, data, sz, first_datum)
    datapadded = PaddedView{T}(args...)

Create a padded version of the array `data`, where any elements within
the span of `padded_axes` not assigned in `data` will have value
`fillvalue`.

Supply `data_axes` to specify an alterate set of axes for `data`, effectively
relocating `data` to a different set of indices.
This is shorthand for

    offsetdata = OffsetArray(data, data_axes)
    datapadded = PaddedView(fillvalue, offsetdata, padded_axes)

using the [OffsetArrays](https://github.com/JuliaArrays/OffsetArrays.jl) package.

Alternately, the padded array size `sz` can be specified, in which case `datapadded`
starts indexing at 1.
One may optionally specify the location of the `[1, 1, ...]` element of `data` with
`first_datum`.
Specifically, `datapadded[first_datum...]` corresponds to `data[1, 1, ...]`.
`first_datum` defaults to all-1s.

The view eltype `T` is optional. If not specified, then in most cases, `T` is inferred to be
`eltype(data)`. In cases when `fillvalue` can't be converted to `eltype(data)`, `T` will be
promoted the one that does. For example, when `fillvalue == nothing` and `eltype(data) == Float32`,
the inferred eltype `T` will be `Union{Nothing, Float32}`.

# Example

```jldoctest
julia> using PaddedViews

julia> a = collect(reshape(1:9, 3, 3))
3×3 $(Array{Int,2}):
 1  4  7
 2  5  8
 3  6  9

julia> PaddedView(-1, a, (4, 5))
4×5 PaddedView(-1, ::$(Array{Int,2}), (Base.OneTo(4), Base.OneTo(5))) with eltype $Int:
  1   4   7  -1  -1
  2   5   8  -1  -1
  3   6   9  -1  -1
 -1  -1  -1  -1  -1

julia> PaddedView(-1, a, (1:5,1:5), (2:4,2:4))
5×5 PaddedView(-1, OffsetArray(::$(Array{Int,2}), 2:4, 2:4), (1:5, 1:5)) with eltype $Int with indices 1:5×1:5:
 -1  -1  -1  -1  -1
 -1   1   4   7  -1
 -1   2   5   8  -1
 -1   3   6   9  -1
 -1  -1  -1  -1  -1

julia> PaddedView(-1, a, (0:4, 0:4))
5×5 PaddedView(-1, ::$(Array{Int,2}), (0:4, 0:4)) with eltype $Int with indices 0:4×0:4:
 -1  -1  -1  -1  -1
 -1   1   4   7  -1
 -1   2   5   8  -1
 -1   3   6   9  -1
 -1  -1  -1  -1  -1

julia> PaddedView(-1, a, (5,5), (2,2))
5×5 PaddedView(-1, OffsetArray(::$(Array{Int,2}), 2:4, 2:4), (Base.OneTo(5), Base.OneTo(5))) with eltype $Int:
 -1  -1  -1  -1  -1
 -1   1   4   7  -1
 -1   2   5   8  -1
 -1   3   6   9  -1
 -1  -1  -1  -1  -1
```
"""
struct PaddedView{T,N,I,A} <: AbstractArray{T,N}
    fillvalue::T
    data::A
    indices::I

    function PaddedView{T,N,I,A}(fillvalue,
                                 data,
                                 indices::NTuple{N,AbstractUnitRange}) where {T,N,I,A}
        ndims(data) == N || throw(DimensionMismatch("data and indices should have the same dimension, instead they're $(ndims(data)) and $N."))
        new{T,N,I,A}(fillvalue, data, indices)
    end
end

function PaddedView(fillvalue::FT, data::AbstractArray{T}, args...) where {FT, T}
    PaddedView{filltype(FT, T)}(fillvalue, data, args...)
end

function PaddedView(fillvalue::FT,
                    data::AbstractArray{T,N},
                    padded_inds::NTuple{N,AbstractUnitRange},
                    data_inds::NTuple{N,AbstractUnitRange}) where {FT,T,N}
    PaddedView{filltype(FT, T)}(fillvalue, data, padded_inds, data_inds)
end

_to_axis(x::Union{OneTo, IdentityUnitRange}) = x
_to_axis(r::AbstractUnitRange) = IdentityUnitRange(r)

function PaddedView{FT}(fillvalue,
                        data::AbstractArray{T,N},
                        indices) where {FT,T,N}
    indsoffset = map(_to_axis, indices)
    PaddedView{FT,N,typeof(indsoffset),typeof(data)}(convert(FT, fillvalue), data, indsoffset)
end

function PaddedView{FT}(fillvalue,
                        data::AbstractArray{T,N},
                        sz::Tuple{Integer,Vararg{Integer}}) where {FT,T,N}
    inds = map(OneTo, sz)
    PaddedView{FT,N,typeof(inds),typeof(data)}(convert(FT, fillvalue), data, inds)
end

# No need to export this
#
# **Regardless of accuracy**, when `data` cannot represent the exact meaning of fillvalue
# object, the type of `fillvalue` should be lifted to a common type.
# Examples of this are: `Nothing`, `Missing`, `ColorTypes.Colorant`
#
# Since eltype of `data` will be lazily promoted to the filltype, it's likely to hit a
# performance issue if we abuse `filltype` by also considering storage type (e.g,
# `Float32`/`Float64`). Thus we shouldn't add lines like:
#     `filltype(::Type{FT}, ::Type{T}) where {FT<:Real, T<:Real} = promote_type(FT, T)`
# Ref: https://github.com/JuliaArrays/PaddedViews.jl/pull/25#issuecomment-610039569
filltype(::Type, ::Type{T}) where T = T
filltype(::Type{FT}, ::Type{T}) where {FT<:Union{Nothing, Missing}, T} = Union{FT, T}
filltype(::Type{FT}, ::Type{T}) where {FT, T<:Union{Nothing, Missing}} = Union{FT, T}
# ambiguity patch
filltype(::Type{FT}, ::Type{T}) where {FT<:Union{Nothing, Missing}, T<:Union{Nothing, Missing}} = Union{FT, T}

# This method eliminates an ambiguity between the two below it
function PaddedView{FT}(fillvalue,
                        data::AbstractArray{T,0},
                        ::Tuple{},
                        ::Tuple{}) where {FT,T}
    return PaddedView{FT}(fillvalue, data, ())
end

function PaddedView{FT}(fillvalue,
                        data::AbstractArray{T,N},
                        padded_inds::NTuple{N,AbstractUnitRange},
                        data_inds::NTuple{N,AbstractUnitRange}) where {FT,T,N}
    off_data = OffsetArray(data, data_inds...)
    return PaddedView{FT}(fillvalue, off_data, padded_inds)
end

function PaddedView{FT}(fillvalue,
                        data::AbstractArray{T,N},
                        sz::NTuple{N,Integer},
                        first_datum::NTuple{N,Integer}) where {FT,T,N}
    padded_inds = map(OneTo, sz)
    data_inds   = map((ax, o)->ax.+o, axes(data), first_datum .- 1)
    return PaddedView{FT}(fillvalue, data, padded_inds, data_inds)
end

Base.axes(A::PaddedView) = A.indices
@inline Base.axes(A::PaddedView, d::Integer) = d <= ndims(A) ? A.indices[d] : default_axes(A.indices)
default_axes(::NTuple{N,I}) where {N,I<:AbstractUnitRange} = convert(I, OneTo(1))
default_axes(::Tuple{}) = OneTo(1)
default_axes(::Any) = OneTo(1)

Base.size(A::PaddedView) = map(length, axes(A))

Base.parent(A::PaddedView) = A.data

Base.@propagate_inbounds function Base.getindex(A::PaddedView{T,N}, i::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, i...)
    if Base.checkbounds(Bool, A.data, i...)
        return convert(T, A.data[i...])
    end
    return A.fillvalue
end

"""
    Aspad = paddedviews(fillvalue, A1, A2, ...; [dims])

Pad the arrays `A1`, `A2`, ..., to a common size or set of axes,
chosen as the span of axes enclosing all of the input arrays.

The padding is applied to one direction in dimensions `dims`. For example, values are filled to
bottom-right part of the new array in two-dimensional case. Use [`sym_paddedviews`](@ref) if
_both_ directions need to be padded.

The axes of original array `A` will be preserved in the padded result `Ap`, hence it's true
that `Ap[CartesianIndices(A)] == A`.

# Example:
```jldoctest
julia> using PaddedViews

julia> a1 = reshape([1, 2, 3], 3, 1)
3×1 $(Array{Int,2}):
 1
 2
 3

julia> a2 = [4 5 6]
1×3 $(Array{Int,2}):
 4  5  6

julia> a1p, a2p = paddedviews(-1, a1, a2);

julia> a1p
3×3 PaddedView(-1, ::$(Array{Int,2}), (Base.OneTo(3), Base.OneTo(3))) with eltype $Int:
 1  -1  -1
 2  -1  -1
 3  -1  -1

julia> a2p
3×3 PaddedView(-1, ::$(Array{Int,2}), (Base.OneTo(3), Base.OneTo(3))) with eltype $Int:
  4   5   6
 -1  -1  -1
 -1  -1  -1

julia> a1p[CartesianIndices(a1)]
3×1 $(Array{Int,2}):
 1
 2
 3
```

`dims` keyword allows padding only for specified dimensions.

```jldoctest; setup=:(using PaddedViews)
julia> a1 = reshape(collect(1:9), 3, 3)
3×3 $(Matrix{Int}):
 1  4  7
 2  5  8
 3  6  9

julia> a2 = [4 5;6 7]
2×2 $(Matrix{Int}):
 4  5
 6  7

julia> a1f, a2f = paddedviews(-1, a1, a2; dims=1);

julia> a2f
3×2 PaddedView(-1, ::$(Matrix{Int}), (Base.OneTo(3), Base.OneTo(2))) with eltype $(Int):
  4   5
  6   7
 -1  -1

julia> a1f, a2f = paddedviews(-1, a1, a2; dims=(1,2));

julia> a2f
3×3 PaddedView(-1, ::$(Matrix{Int}), (Base.OneTo(3), Base.OneTo(3))) with eltype $(Int):
  4   5  -1
  6   7  -1
 -1  -1  -1
```
"""
function paddedviews(fillvalue,  As::AbstractArray...; dims=1:ndims(first(As)))
    inds = outerinds(As...)
    map(A->PaddedView(fillvalue, A, _extended_axes(A, inds, dims)), As)
end
_extended_axes(A, inds, ::Nothing) = inds
function _extended_axes(A, inds, dims::Int)
    map((r1, r2, d)->d ? r2 : r1, axes(A), inds, ntuple(i->i==dims, ndims(A)))
end
function _extended_axes(A, inds, dims::Tuple)
    map((r1, r2, d)->d ? r2 : r1, axes(A), inds, ntuple(i->i in collect(dims), ndims(A)))
end
function _extended_axes(A, inds, dims::UnitRange{Int})
    map((r1, r2, d)->d ? r2 : r1, axes(A), inds, ntuple(i->i in collect(dims), ndims(A)))
end

# Zero, one, and two arrays are common, improve inferrability
paddedviews(fillvalue) = ()
paddedviews(fillvalue, A1::AbstractArray) = (PaddedView(fillvalue, A1, outerinds(A1)),)
function paddedviews(fillvalue, A1::AbstractArray, A2::AbstractArray)
    inds = outerinds(A1, A2)
    PaddedView(fillvalue, A1, inds), PaddedView(fillvalue, A2, inds)
end

# This is an (unexported) optimization if you're supplying the arrays from a vector.
# MosaicViews uses this.
function paddedviews_itr(fillvalue, itr)
    inds = outerinds(itr...)
    [PaddedView(fillvalue, A, inds) for A in itr]
end

@inline outerinds(A::AbstractArray, Bs...) = _outerinds(axes(A), Bs...)
@inline _outerinds(inds, A::AbstractArray, Bs...) =
    _outerinds(__outerinds(inds, axes(A)), Bs...)
_outerinds(inds) = inds
@inline function __outerinds(inds1::NTuple{N,I}, inds2) where {N,I<:AbstractUnitRange}
    map((i1, i2) -> padrange(i1, i2), inds1, inds2)
end

# This shouldn't be reached due to the `paddedviews(fillvalue)` method above,
# but it's here in case anyone extends `paddedviews` for types other than AbstractArrays.
# See https://github.com/JuliaImages/ImageCore.jl/pull/32#discussion_r111545756
outerinds() = error("must supply at least one array with concrete axes")

padrange(i1::T, i2::T) where T<: OneTo = OneTo(max(last(i1), last(i2)))
padrange(i1::T, i2::T) where T<: AbstractUnitRange = convert(T, min(first(i1),first(i2)):max(last(i1), last(i2)))
padrange(i1::AbstractRange, i2::AbstractRange) = padrange(convert(UnitRange{Int}, i1), convert(UnitRange{Int}, i2))



"""
    Aspad = sym_paddedviews(fillvalue, A1, A2, ...; [dims])

Pad the arrays `A1`, `A2`, ..., to a common size or set of axes, chosen as the span of axes
enclosing all of the input arrays.

The padding is applied to both directions in dimensions `dims`, which means original array
located at the center the padded result. Use [`paddedviews`](@ref) if only one direction
need to be padded.

The axes of original array `A` will be preserved in the padded result `Ap`, hence it's true
that `Ap[CartesianIndices(A)] == A`.

```jldoctest
julia> using PaddedViews

julia> a1 = reshape([1, 2, 3], 3, 1)
3×1 $(Array{Int,2}):
 1
 2
 3

julia> a2 = [4 5 6]
1×3 $(Array{Int,2}):
 4  5  6

julia> a1p, a2p = sym_paddedviews(-1, a1, a2);

julia> a1p
3×3 PaddedView(-1, ::$(Array{Int,2}), (1:3, 0:2)) with eltype $Int with indices 1:3×0:2:
 -1  1  -1
 -1  2  -1
 -1  3  -1

julia> a2p
3×3 PaddedView(-1, ::$(Array{Int,2}), (0:2, 1:3)) with eltype $Int with indices 0:2×1:3:
 -1  -1  -1
  4   5   6
 -1  -1  -1

julia> a1p[CartesianIndices(a1)]
3×1 $(Array{Int,2}):
 1
 2
 3
```

`dims` keyword allows padding only for specified dimensions.

```jldoctest; setup=:(using PaddedViews)
julia> a1 = reshape(collect(1:9), 3, 3)
 3×3 $(Matrix{Int}):
  1  4  7
  2  5  8
  3  6  9

julia> a2 = reshape([5, 6], 2, 1)
 2×1 $(Matrix{Int}):
  5
  6

julia> a1f, a2f = sym_paddedviews(-1, a1, a2; dims=1);

julia> a2f
 3×1 PaddedView(-1, ::$(Matrix{Int}), (1:3, 1:1)) with eltype $(Int) with indices 1:3×1:1:
   5
   6
  -1

julia> a1f, a2f = sym_paddedviews(-1, a1, a2; dims=(1,2));

julia> a2f
 3×3 PaddedView(-1, ::$(Matrix{Int}), (1:3, 0:2)) with eltype $(Int) with indices 1:3×0:2:
  -1   5  -1
  -1   6  -1
  -1  -1  -1
 ```
"""
function sym_paddedviews(fillvalue, As::AbstractArray...; dims=1:ndims(first(As)))
    inds = outerinds(As...)
    map(As) do A
        PaddedView(fillvalue, A, _sym_pad_inds(A, inds, dims))
    end
end
sym_paddedviews(fillvalue) = ()
sym_paddedviews(fillvalue, A::AbstractArray) = (PaddedView(fillvalue, A, _sym_pad_inds(A, outerinds(A))),)
function sym_paddedviews(fillvalue, A1::AbstractArray, A2::AbstractArray)
    inds = outerinds(A1, A2)
    PaddedView(fillvalue, A1, _sym_pad_inds(A1, inds)), PaddedView(fillvalue, A2, _sym_pad_inds(A2, inds))
end

function sym_paddedviews_itr(fillvalue, itr)
    inds = outerinds(itr...)
    [PaddedView(fillvalue, A, _sym_pad_inds(A, inds)) for A in itr]
end

function _sym_pad_inds(A, inds, dims=1:ndims(A))
    inds = _extended_axes(A, inds, dims)
    A_axes = axes(A)
    map(A_axes, inds) do ax, i
        pad_sz = length(i) - length(ax)
        offset = pad_sz ÷ 2
        first(ax)-offset:last(ax)+pad_sz-offset
    end
end


function Base.showarg(io::IO, A::PaddedView, toplevel)
    print(io, "PaddedView(", A.fillvalue, ", ")
    Base.showarg(io, parent(A), false)
    print(io, ", (", join(map(no_offset_view, A.indices), ", "))
    print(io, ndims(A) == 1 ? ",))" : "))")
    toplevel && print(io, " with eltype ", eltype(A))
    return nothing
end

end # module
