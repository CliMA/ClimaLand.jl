import .Iterators.drop

istiled(ifd::IFD) = TILEWIDTH in ifd
isplanar(ifd::IFD) = Int(getdata(ifd, PLANARCONFIG, 1)) == 2
tilecols(ifd::IFD) = Int(ifd[TILEWIDTH].data)::Int
tilerows(ifd::IFD) = Int(ifd[TILELENGTH].data)::Int
nrows(ifd::IFD) = Int(ifd[IMAGELENGTH].data)::Int
ncols(ifd::IFD) = Int(ifd[IMAGEWIDTH].data)::Int
nsamples(ifd::IFD) = Int(getdata(ifd, SAMPLESPERPIXEL, 1))
predictor(ifd::IFD) = Int(getdata(ifd, PREDICTOR, 0))
bitspersample(ifd::IFD) = Int(first(ifd[BITSPERSAMPLE].data))::Int
ispalette(ifd::IFD) = Int(getdata(ifd, PHOTOMETRIC, 0)) == 3
compression(ifd::IFD) = getdata(CompressionType, ifd, COMPRESSION, COMPRESSION_NONE)
colortype(ifd::IFD) = interpretation(ifd)

is_irregular_bps(ifd::IFD) = bitspersample(ifd) != sizeof(rawtype(ifd)) * 8
is_complicated(ifd::IFD) = !iscontiguous(ifd) || compression(ifd) != COMPRESSION_NONE || is_irregular_bps(ifd) == true || predictor(ifd) > 1

# returns true if all slices have the same size and color type
function is_homogeneous(ifds::Vector{<:IFD})
    return all(map(==(nrows(first(ifds))), nrows.(ifds))) &&
        all(map(==(ncols(first(ifds))), ncols.(ifds))) &&
        all(map(==(colortype(first(ifds))), colortype.(ifds))) &&
        # tiled images are padded during encoding, so 'nrows' and 'ncols'
        # may not indicate the true space requirements for decoding
        all(map(==(istiled(first(ifds))), istiled.(ifds)))
end

"""
    interpretation(ifd)

For a given IFD, determine the proper colorimetric interpretation of the data.
It returns subtypes of `ColorOrTuple` depending on the values of the tiff tags and
whether there are extrasamples it doesn't know how to deal with.
"""
function interpretation(ifd::IFD)
    interp = PhotometricInterpretations(ifd[PHOTOMETRIC].data)
    bps = Int.(ifd[BITSPERSAMPLE].data)
    extras = EXTRASAMPLES in ifd ? ExtraSamples(first(ifd[EXTRASAMPLES].data)) : missing
    type = interpretation(interp, extras)

    if length(type) == nsamples(ifd)
        return type{_mappedtype(rawtype(ifd), first(bps))}
    else
        sf = SampleFormats.(SAMPLEFORMAT in ifd ? Int.(ifd[SAMPLEFORMAT].data) : fill(1, length(bps)))

        # color type, eg RGB{N0f8}
        ctype = type{_mappedtype(rawtype(ifd), first(bps))}

        # number of components, eg 3 for RGB
        n = length(ctype)

        # tuple type for "extra" channels, dervied from "sample format" and "bits per sample" tags
        extratype = Tuple{_mappedtype.(rawtype.(drop(sf, n), drop(bps, n)), drop(bps, n))...}

        WidePixel{ctype, extratype}
    end
end

interpretation(p::PhotometricInterpretations, x::Union{Missing, ExtraSamples}) = interpretation(Val(p), Val(x))
interpretation(::Val{PHOTOMETRIC_RGB}) = RGB
interpretation(::Val{PHOTOMETRIC_MINISBLACK}) = Gray
interpretation(::Val{PHOTOMETRIC_PALETTE}) = Palette
interpretation(::Val{PHOTOMETRIC_YCBCR}) = YCbCr
interpretation(::Val{PHOTOMETRIC_CIELAB}) = Lab
interpretation(p::Val, ::Val) = interpretation(p)
interpretation(p::Val, ::Val{EXTRASAMPLE_UNASSALPHA}) = interpretation(p)
interpretation(p::Val, ::Val{EXTRASAMPLE_ASSOCALPHA}) = coloralpha(interpretation(p))
interpretation(p::Val, ::Val{EXTRASAMPLE_ASSOCALPHA_NS}) = coloralpha(interpretation(p))

# dummy color type for palette colored images to dispatch on
struct Palette{T} <: Colorant{T, 1}
    i::T
end

Base.adjoint(x::WidePixel) = x

_mappedtype(::Type{T}, bps) where {T} = T
_mappedtype(::Type{T}, bps) where {T <: Unsigned} = Normed{T, bps}
_mappedtype(::Type{T}, bps) where {T <: Signed}   = Fixed{T, bps - 1}

function rawtype(format::SampleFormats, bits::Int)
    bits < 1 || bits > 64 && error("unsupported bit depth ($bits)")

    if format == SAMPLEFORMAT_IEEEFP
        if bits == 16
            rtype = Float16
        elseif bits == 32
            rtype = Float32
        elseif bits == 64
            rtype = Float64
        else
            error("unsupported sample format")
        end
    elseif format == SAMPLEFORMAT_UINT || format == SAMPLEFORMAT_INT
        m = trailing_zeros(nextpow(2, bits))
        if m <= 3
            rtype = format == SAMPLEFORMAT_UINT ? UInt8 : Int8
        elseif m == 4
            rtype = format == SAMPLEFORMAT_UINT ? UInt16 : Int16
        elseif m == 5
            rtype = format == SAMPLEFORMAT_UINT ? UInt32 : Int32
        elseif m == 6
            rtype = format == SAMPLEFORMAT_UINT ? UInt64 : Int64
        else
            error("unsupported sample format")
        end
    else
        error("unsupported sample format")
    end

    @debug "raw type for ($format, $bits) is $rtype"

    rtype
end

function rawtype(ifd::IFD)
    samplesperpixel = nsamples(ifd)
    bps = bitspersample(ifd)
    sampleformats = fill(UInt16(0x01), samplesperpixel)
    if SAMPLEFORMAT in ifd
        sampleformats = ifd[SAMPLEFORMAT].data
    end

    format = SampleFormats(first(sampleformats))
    n = first(bps)

    rawtype(format, n)
end

"""
    $(SIGNATURES)

Allocate a cache for this IFD with correct type and size.
"""
function getcache(ifd::IFD)
    pixeltype = interpretation(ifd)
    if istiled(ifd)
        tile_width = tilecols(ifd)
        tile_height = tilerows(ifd)
        return Array{pixeltype}(undef, cld(ncols(ifd), tile_width) * tile_width, cld(nrows(ifd), tile_height) * tile_height)
    else
        return Array{pixeltype}(undef, ncols(ifd), nrows(ifd))
    end
end

function uncompressed_size(ifd::IFD, columns::Int, rows::Int)
    bps = bitspersample(ifd)
    # each row is encoded separately
    cld(columns * bps, 8) * rows
end
