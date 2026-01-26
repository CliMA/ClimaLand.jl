module ImageCore

using Reexport
@reexport using FixedPointNumbers
@reexport using Colors
@reexport using ColorVectorSpace

@reexport using MosaicViews
@reexport using PaddedViews
using MappedArrays
using OffsetArrays # for show.jl
using .ColorTypes: colorant_string
using Colors: Fractional
using MappedArrays: AbstractMultiMappedArray

using Base: tail, @pure, Indices
import Base: float

# TODO: just use .+
# See https://github.com/JuliaLang/julia/pull/22932#issuecomment-330711997
plus(r::AbstractUnitRange, i::Integer) = broadcast(+, r, i)
plus(a::AbstractArray, i::Integer) = a .+ i

using .ColorTypes: AbstractGray, TransparentGray, Color3, Transparent3
const Color1{T} = Colorant{T,1}
const Color2{T} = Colorant{T,2}
const Color4{T} = Colorant{T,4}
const AColor{N,C,T} = AlphaColor{C,T,N}
const ColorA{N,C,T} = ColorAlpha{C,T,N}
const NonparametricColors = Union{RGB24,ARGB32,Gray24,AGray32}
const Color1Array{C<:Color1,N} = AbstractArray{C,N}
# Type that arises from reshape(reinterpret(To, A), sz):
const RRArray{To,From,M,P} = Base.ReinterpretArray{To,M,From,P,true}
const RGArray = Union{Base.ReinterpretArray{<:AbstractGray,M,<:Number,P}, Base.ReinterpretArray{<:Number,M,<:AbstractGray,P}} where {M,P}

# Deliberately not export these constants to enable extensibility for downstream packages
const NumberLike = Union{Number,AbstractGray}
const Pixel = Union{Number,Colorant}
const GenericGrayImage{T<:NumberLike,N} = AbstractArray{T,N}
const GenericImage{T<:Pixel,N} = AbstractArray{T,N}

export
    ## Types
    HasDimNames,
    HasProperties,
    StackedView,
    ## constants
    zeroarray,
    ## functions
    # views
    channelview,
    colorview,
    permuteddimsview,
    rawview,
    normedview,
    reinterpretc,
    # conversions
#    float16,
    float32,
    float64,
    n0f8,
    n6f10,
    n4f12,
    n2f14,
    n0f16,
    # mapping values
    clamp01,
    clamp01!,
    clamp01nan,
    clamp01nan!,
    colorsigned,
    scaleminmax,
    scalesigned,
    takemap,
    # traits
    assert_timedim_last,
    coords_spatial,
    indices_spatial,
    namedaxes,
    nimages,
    pixelspacing,
    sdims,
    size_spatial,
    spacedirections,
    spatialorder

include("colorchannels.jl")
include("stackedviews.jl")
include("convert_reinterpret.jl")
include("traits.jl")
include("map.jl")
include("show.jl")
include("deprecations.jl")

"""
    rawview(img::AbstractArray{FixedPoint})

returns a "view" of `img` where the values are interpreted in terms of
their raw underlying storage. For example, if `img` is an `Array{N0f8}`,
the view will act like an `Array{UInt8}`.

See also: [`normedview`](@ref)
"""
rawview(a::AbstractArray{T}) where {T<:FixedPoint} = mappedarray(reinterpret, y->T(y,0), a)
rawview(a::Array{T}) where {T<:FixedPoint} = reinterpret(FixedPointNumbers.rawtype(T), a)
rawview(a::AbstractArray{T}) where {T<:Real} = a

"""
    normedview([T], img::AbstractArray{Unsigned})

returns a "view" of `img` where the values are interpreted in terms of
`Normed` number types. For example, if `img` is an `Array{UInt8}`, the
view will act like an `Array{N0f8}`.  Supply `T` if the element
type of `img` is `UInt16`, to specify whether you want a `N6f10`,
`N4f12`, `N2f14`, or `N0f16` result.

See also: [`rawview`](@ref)
"""
normedview(::Type{T}, a::AbstractArray{S}) where {T<:FixedPoint,S<:Unsigned} = mappedarray(y->T(y,0),reinterpret, a)
normedview(::Type{T}, a::Array{S}) where {T<:FixedPoint,S<:Unsigned} = reinterpret(T, a)
normedview(::Type{T}, a::AbstractArray{T}) where {T<:Normed} = a
normedview(a::AbstractArray{UInt8}) = normedview(N0f8, a)
normedview(a::AbstractArray{T}) where {T<:Normed} = a

# PaddedViews support
# This make sure Colorants as `fillvalue` are correctly filled, for example, let
# `PaddedView(ARGB(0, 0, 0, 0), img)` correctly filled with transparent color even when
# `img` is of eltype `RGB`
function PaddedViews.filltype(::Type{FC}, ::Type{C}) where {FC<:Colorant, C<:Colorant}
    # rand(RGB, 4, 4) has eltype RGB{Any} but it isn't a concrete type
    # although the consensus[1] is to not make a concrete eltype, this op is needed to make a
    # type-stable colorant construction in _filltype without error; there's no RGB{Any} thing
    # [1]: https://github.com/JuliaLang/julia/pull/34948
    T = eltype(C) === Any ? eltype(FC) : eltype(C)
    _filltype(FC, base_colorant_type(C){T})
end
PaddedViews.filltype(::Type{FC}, ::Type{C}) where {FC<:Colorant, C<:Number} = base_color_type(FC){C}
PaddedViews.filltype(::Type{FC}, ::Type{<:Number}) where {FC<:Union{Gray24, AGray32, RGB24, ARGB32}} = FC
_filltype(::Type{<:Colorant}, ::Type{C}) where {C<:Colorant} = C
_filltype(::Type{FC}, ::Type{C}) where {FC<:Color3, C<:AbstractGray} =
    base_colorant_type(FC){promote_type(eltype(FC), eltype(C))}
_filltype(::Type{FC}, ::Type{C}) where {FC<:TransparentColor, C<:AbstractGray} =
    alphacolor(FC){promote_type(eltype(FC), eltype(C))}
_filltype(::Type{FC}, ::Type{C}) where {FC<:TransparentColor, C<:Color3} =
    alphacolor(C){promote_type(eltype(FC), eltype(C))}


# MosaicViews support (require MosaicViews >= v0.3.0)
MosaicViews.promote_wrapped_type(::Type{T}, ::Type{S}) where {T<:Colorant,S<:Colorant} = promote_type(T, S)
MosaicViews.promote_wrapped_type(::Type{T}, ::Type{S}) where {T,S<:Colorant} = S === Union{} ? T : base_colorant_type(S){MosaicViews.promote_wrapped_type(T, eltype(S))}
MosaicViews.promote_wrapped_type(::Type{T}, ::Type{S}) where {T<:Colorant,S} = MosaicViews.promote_wrapped_type(S, T)

# Support transpose
Base.transpose(a::AbstractMatrix{C}) where {C<:Colorant} = permutedims(a, (2,1))
function Base.transpose(a::AbstractVector{C}) where C<:Colorant
    ind = axes(a, 1)
    out = similar(Array{C}, (oftype(ind, Base.OneTo(1)), ind))
    outr = reshape(out, ind)
    copy!(outr, a)
    out
end

# It's better not to define fft on Colorant arrays, because keeping
# track of the color dimension and the fft-dims is prone to omissions
# or problems due to later operations. So we put the bookkeeping on
# the user, but we try to give helpful suggestions.
function throw_ffterror(io, @nospecialize(f), x, dims=1:ndims(x))
    newdims = plus(dims, channelview_dims_offset(x))
    print(io, '\n', f, " not defined for eltype $(eltype(x)). Use channelview, and likely $newdims for the dims in the fft.")
end

function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        f = exc.f
        if isa(f, Function) && nameof(f) ∈ (:fft, :rfft, :plan_fft, :plan_rfft, :realfloat) && argtypes[1] <: AbstractArray{<:Colorant}
            throw_ffterror(io, exc.f, exc.args...)
        end
    end
end

include("precompile.jl")

end ## module
