abstract type AbstractTIFF{T,N} <: AbstractArray{T,N} end

abstract type AbstractDenseTIFF{T,N} <: AbstractTIFF{T,N} end

abstract type AbstractStridedTIFF{T,N} <: AbstractTIFF{T,N} end

const ColorOrTuple = Union{WidePixel,Colorant}

function Base.getproperty(img::T, sym::Symbol) where {T<:AbstractTIFF}
    if sym === :ifds
        Base.depwarn("Directly accessing the ifds property is deprecated" *
                     " and will be removed in TiffImages v0.7.0+. Please use `ifds(img)` instead.",
            :AbstractTIFF)
    end
    getfield(img, sym)
end

"""
    ifds(img)

Get the Image File Directories of `img`, see [`TiffImages.IFD`] for details.
Returns a single IFD for a 2D TIFF. Otherwise, returns a list of IFDs with a
length equal to the third dimension of a 3D TIFF.
```
"""
ifds(img::I) where {I<:AbstractTIFF{T,2} where {T}} = first(getfield(img, :ifds))
ifds(img::I) where {I<:AbstractTIFF} = getfield(img, :ifds)

interpretation(img::AbstractArray) = interpretation(eltype(img))
interpretation(::IndirectArray{T,N,I}) where {T,N,I} = PHOTOMETRIC_PALETTE
interpretation(::SubArray{T,N,P}) where {T,N,P<:IndirectArray{X,Y,I}} where {X,Y,I} = PHOTOMETRIC_PALETTE
interpretation(::Type{WidePixel{C,X}}) where {C,X} = interpretation(C)
interpretation(::Type{T}) where {T<:Gray} = PHOTOMETRIC_MINISBLACK
interpretation(::Type{T}) where {T<:AbstractRGB} = PHOTOMETRIC_RGB
interpretation(::Type{<:TransparentColor{C,T,N}}) where {C,T,N} = interpretation(C)
interpretation(::Type) = PHOTOMETRIC_MINISBLACK

samplesperpixel(img::AbstractArray) = samplesperpixel(eltype(img))
samplesperpixel(::IndirectArray{T,N,I}) where {T,N,I} = 1
samplesperpixel(::SubArray{T,N,P}) where {T,N,P<:IndirectArray{X,Y,I}} where {X,Y,I} = 1
samplesperpixel(::Type{<:Colorant{T,N}}) where {T,N} = N
samplesperpixel(t::Type{<:WidePixel{C,X}}) where {C,X} = samplesperpixel(C) + length(fieldnames(X))

bitspersample(img::AbstractArray) = bitspersample(eltype(img))
bitspersample(::IndirectArray{T,N,I}) where {T,N,I} = I <: Integer ? bitspersample(I) : bitspersample(eltype(I)) # I is an integer after IndirectArray v1
bitspersample(::SubArray{T,N,P}) where {T,N,P<:IndirectArray{X,Y,I}} where {X,Y,I} = bitspersample(I)
bitspersample(T::Type{<:WidePixel}) = collect(flatten(bitspersample.(fieldtypes(T))))
bitspersample(T::Type{<:Colorant}) = collect(bitspersample.(fieldtypes(T)))
bitspersample(T::Type{<:Tuple}) = collect(bitspersample.(fieldtypes(T)))
bitspersample(::Type{<:Normed{T,S}}) where {T,S} = S
bitspersample(::Type{<:Fixed{T,S}}) where {T,S} = S + 1
bitspersample(::Type{T}) where T = sizeof(T) * 8

sampleformat(img::AbstractArray) = sampleformat(eltype(img))
sampleformat(T::Type{<:WidePixel}) = collect(flatten(sampleformat.(fieldtypes(T))))
sampleformat(T::Type{<:Colorant}) = collect(sampleformat.(fieldtypes(T)))
sampleformat(T::Type{<:Tuple}) = collect(sampleformat.(fieldtypes(T)))
sampleformat(::Type{<:Normed{T,S}}) where {T,S} = 1
sampleformat(::Type{<:Fixed{T,S}}) where {T,S} = 2
sampleformat(::Type{T}) where T<:AbstractFloat = 3
sampleformat(::Type{Complex{T}}) where T<:Signed = 5
sampleformat(::Type{Complex{T}}) where T<:AbstractFloat = 6

extrasamples(img::AbstractArray) = extrasamples(eltype(img))
extrasamples(T::Type{WidePixel{C,X}}) where {C<:TransparentColor,X} = vcat([1], fill(0, length(fieldtypes(X))))
extrasamples(T::Type{WidePixel{C,X}}) where {C,X} = fill(0, length(fieldtypes(X)))
extrasamples(::Type) = nothing
extrasamples(img::Type{<:TransparentColor}) = 1

"""
    channel(img, i)

Get the `i`'th channel of data from `img`

For example, `channel(img, 2)` would yield the green channel
of an RGB image, or the (unlabeled) second channel of a multispectral
image
"""
channel(img::AbstractTIFF, i::Int) = channel.(img, i)
channel(x::Colorant, i::Int) = getfield(x, i)
channel(x::WidePixel, i::Int) = i <= length(x.color) ? getfield(x.color, i) : x.extra[i-length(x.color)]

_length(x) = length(x)
_length(T::Type{<:Tuple}) = length(fieldnames(T))

"""

    nchannels(img::AbstractTIFF)

Return the number of channels in each pixel

For example, an image with RGB pixels has 3 channels
"""
nchannels(img::AbstractTIFF) = nchannels(eltype(img))
nchannels(x) = _length(x)
nchannels(::Type{WidePixel{C,X}}) where {C,X} = _length(C) + _length(X)
nchannels(::WidePixel{C,X}) where {C,X} = _length(C) + _length(X)

"""
    color(img::AbstractTIFF{ <: WidePixel}; alpha)

Extract the color channels from multispectral image `img`

The optional `alpha` parameter allows an arbitrary channel to be used as the alpha channel

```jldoctest; setup=:(import TiffImages: DenseTaggedImage, WidePixel; import ColorTypes: RGB, RGBA)
julia> img = DenseTaggedImage(WidePixel.([rand(RGB{Float32}) for x in 1:256, y in 1:256], [(rand(Float32),) for x in 1:256, y in 1:256]));

julia> eltype(img)
WidePixel{RGB{Float32}, Tuple{Float32}}

julia> nchannels(img)
4

julia> eltype(color(img; alpha=4)) # use the fourth channel as an alpha channel
RGBA{Float32}
```
"""
color(img::AbstractTIFF{<:WidePixel}; alpha::Union{Nothing,Integer}=nothing) = color.(img; alpha)

"""
    color(img::AbstractTIFF)

This is an identity relation that just returns `img`
"""
color(img::AbstractTIFF) = img

"""
    color(x::WidePixel; alpha)

Return the value of the `color` field of `x`

The optional `alpha` parameter allows an arbitrary channel to be used as the alpha channel
"""
function color(x::WidePixel{C,X}; alpha::Union{Nothing,Integer}=nothing) where {C,X}
    alpha == nothing ? x.color : ColorTypes.coloralpha(C)(ColorTypes.color(x.color), channel(x, alpha))
end

"""
    color(x::Colorant)

This is an identity relation that just returns `x`
"""
color(x::Colorant) = x
