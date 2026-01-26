module ImageAxes

using Base: @pure, tail
using Reexport, SimpleTraits

# maybe return to "@reexport AxisArrays" if AxisArrays is fixed

import AxisArrays
using AxisArrays: AxisArray, Axis, axisnames, axisvalues, axisdim, atindex, atvalue, collapse
export AxisArray, Axis, axisnames, axisvalues, axisdim, atindex, atvalue, collapse

@reexport using ImageCore
using ImageCore.MappedArrays
using ImageCore.OffsetArrays
import ImageBase: restrict


export # types
    HasTimeAxis,
    IndexAny,
    IndexIncremental,
    StreamingContainer,
    TimeAxis,
    StreamIndexStyle,
    # functions
    colordim,
    arraydata,
    getindex!,
    istimeaxis,
    timeaxis,
    timedim

"""
    TimeAxis{Ax}

A trait (from SimpleTraits) indicating whether axis `Ax` corresponds
to time. This decision is based on the symbol-name given to `Ax`. For
example, the following declares that all `Axis{:time}` objects
correspond to time:

    @traitimpl TimeAxis{Axis{:time}}

This definition has already been made in ImageAxes, but you can add
new names as well.
"""
@traitdef TimeAxis{X}

@traitimpl TimeAxis{Axis{:time}}

# Note: any axis not marked as a TimeAxis is assumed to correspond to
# space. It might be useful to allow users to add other possibilities,
# but this is not currently possible with SimpleTraits.

"""
    timeaxis(A)

Return the time axis, if present, of the array `A`, and `nothing` otherwise.
"""
@inline timeaxis(A::AxisArray) = _timeaxis(A.axes...)
timeaxis(A::AbstractArray) = nothing
timeaxis(A::AbstractMappedArray) = timeaxis(parent(A))
@traitfn _timeaxis(ax::Ax, axes...) where {Ax<:Axis; !TimeAxis{Ax}} = _timeaxis(axes...)
@traitfn _timeaxis(ax::Ax, axes...) where {Ax<:Axis;  TimeAxis{Ax}} = ax
_timeaxis() = nothing


"""
    istimeaxis(ax)

Test whether the axis `ax` corresponds to time.
"""
istimeaxis(ax::Axis) = istimeaxis(typeof(ax))
@traitfn istimeaxis(::Type{Ax}) where {Ax<:Axis; !TimeAxis{Ax}} = false
@traitfn istimeaxis(::Type{Ax}) where {Ax<:Axis;  TimeAxis{Ax}} = true

@traitdef HasTimeAxis{X}
"""
    HasTimeAxis{AA}

A trait for testing whether type `AA` has a time axis. Time axes must
be declared before use.

# Examples

```julia
using ImageAxes, SimpleTraits

# Declare that all axes named `:time` are time axes
@traitimpl TimeAxis{Axis{:time}}

# Define functions that dispatch on AxisArrays that may or may not have time axes
@traitfn got_time{AA<:AxisArray;  HasTimeAxis{AA}}(img::AA) = "yep, I've got time"
@traitfn got_time{AA<:AxisArray; !HasTimeAxis{AA}}(img::AA) = "no, I'm too busy"

julia> A = AxisArray(1:5, Axis{:time}(1:5));

julia> got_time(A)
"yep, I've got time"

julia> A = AxisArray(1:5, Axis{:x}(1:5));

julia> got_time(A)
"no, I'm too busy"
```
"""
HasTimeAxis

axtype(::Type{AxisArray{T,N,D,Ax}}) where {T,N,D,Ax} = Ax
axtype(A::AxisArray) = axtype(typeof(A))

Base.@pure function SimpleTraits.trait(t::Type{HasTimeAxis{AA}}) where AA<:AxisArray
    axscan = map(S->istimeaxis(S), axtype(AA).parameters)
    any(axscan) ? HasTimeAxis{AA} : Not{HasTimeAxis{AA}}
end

# Specializations to preserve the AxisArray wrapper
function Base.PermutedDimsArray(A::AxisArray, perm)
    axs = AxisArrays.axes(A)
    AxisArray(PermutedDimsArray(A.data, perm), axs[[perm...]]...)
end
function ImageCore.channelview(A::AxisArray)
    Ac = channelview(A.data)
    _channelview(A, Ac)
end
# without extra dimension:
_channelview(A::AxisArray{C,N}, Ac::AbstractArray{T,N}) where {C,T,N} = AxisArray(Ac, AxisArrays.axes(A)...)
# with extra dimension:
_channelview(A::AxisArray, Ac::AbstractArray) = AxisArray(Ac, Axis{:color}(axes(Ac,1)), AxisArrays.axes(A)...)


### Image properties based on traits ###

"""
    timedim(img) -> d::Int

Return the dimension of the array used for encoding time, or 0 if not
using an axis for this purpose.

Note: if you want to recover information about the time axis, it is
generally better to use `timeaxis`.
"""
timedim(img::AxisArray{T,N}) where {T,N} = _timedim(filter_time_axis(AxisArrays.axes(img), ntuple(identity, Val(N))))
_timedim(dim::Tuple{Int}) = dim[1]
_timedim(::Tuple{}) = 0

ImageCore.nimages(img::AxisArray) = _nimages(timeaxis(img))
_nimages(::Nothing) = 1
_nimages(ax::Axis) = length(ax)

function colordim(img::AxisArray)
    d = _colordim(1, AxisArrays.axes(img))
    d > ndims(img) ? 0 : d
end
_colordim(d, ax::Tuple{Ax,Vararg{Any}}) where {Ax<:Axis{:color}} = d
_colordim(d, ax::Tuple{Any,Vararg{Any}}) = _colordim(d+1, tail(ax))
_colordim(d, ax::Tuple{}) = d+1

ImageCore.pixelspacing(img::AxisArray) = map(step, filter_space_axes(AxisArrays.axes(img), axisvalues(img)))

ImageCore.spacedirections(img::AxisArray) = ImageCore._spacedirections(pixelspacing(img))

ImageCore.coords_spatial(img::AxisArray{T,N}) where {T,N} = filter_space_axes(AxisArrays.axes(img), ntuple(identity, Val(N)))

ImageCore.spatialorder(img::AxisArray) = filter_space_axes(AxisArrays.axes(img), axisnames(img))

ImageCore.size_spatial(img::AxisArray)    = filter_space_axes(AxisArrays.axes(img), size(img))
ImageCore.indices_spatial(img::AxisArray) = filter_space_axes(AxisArrays.axes(img), axes(img))

arraydata(img::AxisArray) = img.data

### Utilities for writing "simple algorithms" safely ###

# Check that the time dimension, if present, is last
@traitfn function ImageCore.assert_timedim_last(img::AA) where {AA<:AxisArray; HasTimeAxis{AA}}
    ax = AxisArrays.axes(img)[end]
    istimeaxis(ax) || error("time dimension is not last")
    nothing
end
@traitfn ImageCore.assert_timedim_last(img::AA) where {AA<:AxisArray; !HasTimeAxis{AA}} = nothing

### Convert ###
function Base.convert(::Type{Array{C,n}},
                      img::AxisArray{C,n}) where {C<:Colorant,n}
    copyto!(Array{C}(undef, size(img)), img)
end
function Base.convert(::Type{Array{Cdest,n}},
                      img::AxisArray{Csrc,n}) where {Cdest<:Colorant,n,Csrc<:Colorant}
    copyto!(Array{ccolor(Cdest, Csrc)}(undef, size(img)), img)
end

### StreamingContainer ###

checknames(axnames, ::Type{P}) where {P} = checknames(axnames, axisnames(P))
@noinline function checknames(axnames, parentnames::Tuple{Symbol,Vararg{Symbol}})
    mapreduce(x->in(x, parentnames), &, axnames) || throw(DimensionMismatch("names $axnames are not included among $parentnames"))
    nothing
end

"""
    A = StreamingContainer{T}(f!, parent, streamingaxes::Axis...)

An array-like object possessing one or more axes for which changing "slices" may
be expensive or subject to restrictions. A canonical example would be
an AVI stream, where addressing pixels within the same frame is fast
but jumping between frames might be slow.

Here's a simple example of dividing by the mean of each slice of an image before returning values.

    A = AxisArrays.AxisArray(reshape(1:36, 3, 3, 4))

    function f!(buffer, slice)
        meanslice = mean(slice)
        buffer .= slice./meanslice
    end

    B = StreamingContainer{Float64}(f!, A, AxisArrays.axes(A)[3])

    julia> A[:,:,1]
    3×3 AxisArray{Int64,2,Array{Int64,2},Tuple{Axis{:row,Base.OneTo{Int64}},Axis{:col,Base.OneTo{Int64}}}}:
     1  4  7
     2  5  8
     3  6  9

    julia> B[:,:,1]
    3×3 Array{Float64,2}:
     0.2  0.8  1.4
     0.4  1.0  1.6
     0.6  1.2  1.8

The user-provided `f!` function should take arguments:

    f!(buffer, slice)

Where `buffer` will be an empty array that can hold a slice of your series, and `slice` will hold the current input slice.

It's worth noting that `StreamingContainer` is *not* a subtype of
`AbstractArray`, but that much of the array interface (`eltype`,
`ndims`, `axes`, `size`, `getindex`, and `IndexStyle`) is
supported. A StreamingContainer `A` can be built from an AxisArray,
but it may also be constructed from other "parent" objects, even
non-arrays, as long as they support the same functions. In either
case, the parent should also support the standard AxisArray functions
`axes`, `axisnames`, `axisvalues`, and `axisdim`; this support will be
extended to the `StreamingContainer`.

Additionally, a StreamingContainer `A` supports

    getindex!(dest, A, axt::Axis{:time}, ...)

to obtain slices along the streamed axes (here it is assumed that
`:time` is a streamed axis of `A`). You can implement this directly
(dispatching on the parameters of `A`), or (if the parent is an
`AbstractArray`) rely on the fallback

    A.getindex!(dest, view(parent, axs...))

where `A.getindex! = f!` as passed as an argument at construction. `dest` should
have dimensionality `ndims(parent)-length(streamingaxes)`.

Optionally, define [`StreamIndexStyle(typeof(parent),typeof(f!))`](@ref).
"""
struct StreamingContainer{T,N,streamingaxisnames,P,GetIndex}
    getindex!::GetIndex
    parent::P
end
function StreamingContainer{T}(f!::Function, parent, axs::Axis...) where T
    N = ndims(parent)
    axnames = axisnames(axs...)
    checknames(axnames, typeof(parent))
    StreamingContainer{T,N,axnames,typeof(parent),typeof(f!)}(f!, parent)
end

Base.parent(S::StreamingContainer) = S.parent
Base.axes(S::StreamingContainer) = axes(S.parent)
Base.size(S::StreamingContainer)    = size(S.parent)
Base.axes(S::StreamingContainer, d) = axes(S.parent, d)
Base.size(S::StreamingContainer, d)    = size(S.parent, d)

AxisArrays.axes(S::StreamingContainer) = AxisArrays.axes(parent(S))
AxisArrays.axisnames(S::StreamingContainer)  = axisnames(AxisArrays.axes(S)...)
AxisArrays.axisvalues(S::StreamingContainer) = axisvalues(AxisArrays.axes(S)...)
function AxisArrays.axisdim(S::StreamingContainer, ::Type{Axis{name}}) where name
    isa(name, Int) && return name <= ndims(S) ? name : error("axis $name greater than array dimensionality $(ndims(S))")
    names = axisnames(S)
    idx = findfirst(isequal(name), names)
    idx isa Nothing && error("axis $name not found in array axes $names")
    idx
end
AxisArrays.axisdim(S::StreamingContainer, ax::Axis) = axisdim(S, typeof(ax))
AxisArrays.axisdim(S::StreamingContainer, ::Type{Axis{name,T}}) where {name,T} = axisdim(S, Axis{name})

ImageCore.nimages(S::StreamingContainer) = _nimages(timeaxis(S))
ImageCore.coords_spatial(S::StreamingContainer{T,N}) where {T,N} =
    filter_space_axes(AxisArrays.axes(S), ntuple(identity, Val(N)))
ImageCore.spatialorder(S::StreamingContainer) = filter_space_axes(AxisArrays.axes(S), axisnames(S))
ImageCore.size_spatial(img::StreamingContainer)    = filter_space_axes(AxisArrays.axes(img), size(img))
ImageCore.indices_spatial(img::StreamingContainer) = filter_space_axes(AxisArrays.axes(img), axes(img))
function ImageCore.assert_timedim_last(S::StreamingContainer)
    istimeaxis(AxisArrays.axes(S)[end]) || error("time dimension is not last")
    nothing
end

Base.eltype(::Type{StreamingContainer{T,N,names,P,GetIndex}}) where {T,N,names,P,GetIndex} = T
Base.eltype(S::StreamingContainer) = eltype(typeof(S))
Base.ndims(::Type{StreamingContainer{T,N,names,P,GetIndex}}) where {T,N,names,P,GetIndex} = N
Base.ndims(S::StreamingContainer) = ndims(typeof(S))
Base.length(S::StreamingContainer) = prod(size(S))

streamingaxisnames(::Type{StreamingContainer{T,N,names,P,GetIndex}}) where {T,N,names,P,GetIndex} =
    names
streamingaxisnames(S::StreamingContainer) = streamingaxisnames(typeof(S))

timeaxis(S::StreamingContainer) = _timeaxis(AxisArrays.axes(S)...)

@inline function getindex!(dest, S::StreamingContainer, axs::Axis...)
    all(ax->isstreamedaxis(ax,S), axs) || throw(ArgumentError("$axs do not coincide with the streaming axes $(streamingaxisnames(S))"))
    _getindex!(dest, S.getindex!, parent(S), axs...)
end

@inline function getindex!(dest, S::StreamingContainer, I...)
    axs = getslicedindices(S, I)
    all(x->isa(x, Colon), filter_notstreamed(I, S)) || throw(ArgumentError("positional indices must use `:` for any non-streaming axes"))
    _getindex!(dest, S.getindex!, parent(S), axs...)
end

# _getindex! just makes it easy to specialize for particular parent or function types
@inline function _getindex!(dest, f!, P::AbstractArray, axs::Axis...)
    v = view(P, axs...)
    f!(dest, v)
end

function isstreamedaxis(ax::Axis{name},
                        S::StreamingContainer{T,N,saxnames}) where {name,T,N,saxnames}
    in(name, saxnames)
end

sliceindices(S::StreamingContainer) = filter_notstreamed(axes(S), S)
sliceaxes(S::StreamingContainer) = filter_notstreamed(AxisArrays.axes(S), S)
getslicedindices(S::StreamingContainer, I) = filter_streamed(map((ax, i) -> ax(i), AxisArrays.axes(S), I), S)

filter_streamed(inds, S)    = _filter_streamed(inds, AxisArrays.axes(S), S)
filter_notstreamed(inds, S) = _filter_notstreamed(inds, AxisArrays.axes(S), S)
filter_streamed(inds::Tuple{Axis,Vararg{Axis}}, S)    = _filter_streamed(inds, inds, S)
filter_notstreamed(inds::Tuple{Axis,Vararg{Axis}}, S) = _filter_notstreamed(inds, inds, S)

@generated function _filter_streamed(a, axs::NTuple{N,Axis}, S::StreamingContainer) where N
    inds = findall(x->x in streamingaxisnames(S), axisnames(axs.parameters...))
    Expr(:tuple, Expr[:(a[$i]) for i in inds]...)
end
@generated function _filter_notstreamed(a, axs::NTuple{N,Axis}, S::StreamingContainer) where N
    inds = findall(x->x in streamingaxisnames(S), axisnames(axs.parameters...))
    inds = setdiff(1:N, inds)
    Expr(:tuple, Expr[:(a[$i]) for i in inds]...)
end

@inline function Base.getindex(S::StreamingContainer{T,N}, inds::Vararg{Union{Colon,Base.ViewIndex},N}) where {T,N}
    tmp = similar(Array{T}, sliceindices(S))
    getindex!(tmp, S, getslicedindices(S, inds)...)
    tmp[filter_notstreamed(inds, S)...]
end
@inline function Base.getindex(S::StreamingContainer, ind1::Axis, inds_rest::Axis...)
    axs = sliceaxes(S)
    tmp = AxisArray(Array{eltype(S)}(undef, map(length, axs)), axs)
    inds = (ind1, inds_rest...)
    getindex!(tmp, S, _filter_streamed(inds, inds, S)...)
    getindex_rest(tmp, _filter_notstreamed(inds, inds, S))
end
getindex_rest(tmp, ::Tuple{}) = tmp
getindex_rest(tmp, inds) = tmp[inds...]

"""
    style = StreamIndexStyle(A)

A trait that indicates the degree of support for indexing the streaming axes
of `A`. Choices are [`IndexAny()`](@ref) and
[`IndexIncremental()`](@ref) (for arrays that only permit advancing
the time axis, e.g., a video stream from a webcam). The default value
is `IndexAny()`.

This should be specialized for the type rather than the instance. For
a StreamingContainer `S`, you can define this trait via

```julia
StreamIndexStyle(::Type{P}, ::typeof(f!)) = IndexIncremental()
```

where `P = typeof(parent(S))`.
"""
abstract type StreamIndexStyle end

"""
    IndexAny()

Indicates that an axis supports full random-access indexing.
"""
struct IndexAny <: StreamIndexStyle end
"""
    IndexIncremental()

Indicates that an axis supports only incremental indexing, i.e., from `i` to `i+1`.
This is commonly used for the temporal axis with media streams.
"""
struct IndexIncremental <: StreamIndexStyle end

StreamIndexStyle(::Type{A}) where {A<:AbstractArray} = IndexAny()
StreamIndexStyle(A::AbstractArray) = StreamIndexStyle(typeof(A))

StreamIndexStyle(::Type{StreamingContainer{T,N,axnames,P,GetIndex}}) where {T,N,axnames,P,GetIndex} = StreamIndexStyle(P, GetIndex)
StreamIndexStyle(::Type{P},::Type{GetIndex}) where {P,GetIndex} = IndexAny()

StreamIndexStyle(S::StreamingContainer) = StreamIndexStyle(typeof(S))

### Low level utilities ###

filter_space_axes(axes::NTuple{N,Axis}, items::NTuple{N,Any}) where {N} =
    _filter_space_axes(axes, items)
@inline @traitfn _filter_space_axes(axes::Tuple{Ax,Vararg{Any}}, items) where {Ax<:Axis;  TimeAxis{Ax}} =
    _filter_space_axes(tail(axes), tail(items))
@inline @traitfn _filter_space_axes(axes::Tuple{Ax,Vararg{Any}}, items) where {Ax<:Axis; !TimeAxis{Ax}} =
    (items[1], _filter_space_axes(tail(axes), tail(items))...)
_filter_space_axes(::Tuple{}, ::Tuple{}) = ()
@inline _filter_space_axes(axes::Tuple{Ax,Vararg{Any}}, items) where {Ax<:Axis{:color}} =
    _filter_space_axes(tail(axes), tail(items))

filter_time_axis(axes::NTuple{N,Axis}, items::NTuple{N}) where {N} =
    _filter_time_axis(axes, items)
@inline @traitfn _filter_time_axis(axes::Tuple{Ax,Vararg{Any}}, items) where {Ax<:Axis; !TimeAxis{Ax}} =
    _filter_time_axis(tail(axes), tail(items))
@inline @traitfn _filter_time_axis(axes::Tuple{Ax,Vararg{Any}}, items) where {Ax<:Axis;  TimeAxis{Ax}} =
    (items[1], _filter_time_axis(tail(axes), tail(items))...)
_filter_time_axis(::Tuple{}, ::Tuple{}) = ()

# summary: print color types & fixed-point types compactly
function AxisArrays._summary(io, A::AxisArray{T,N}) where {T<:Union{Fractional,Colorant},N}
    print(io, "$N-dimensional AxisArray{")
    if T<:Colorant
        ColorTypes.colorant_string_with_eltype(io, T)
    else
        ColorTypes.showcoloranttype(io, T)
    end
    println(io, ",$N,...} with axes:")
end

include("deprecated.jl")

# glue codes
include("offsetarrays.jl")
include("restrict.jl")

end # module
