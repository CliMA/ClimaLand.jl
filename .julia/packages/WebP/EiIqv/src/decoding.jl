function decode(
    ::Type{TColor}, data::AbstractVector{UInt8}; transpose = false
)::Matrix{TColor} where {TColor <: Colorant}
    TDecodedColor = TColor
    if TColor == ARGB{N0f8}
        webp_decode_fn = Wrapper.WebPDecodeARGB
    elseif TColor == BGR{N0f8}
        webp_decode_fn = Wrapper.WebPDecodeBGR
    elseif TColor == BGRA{N0f8}
        webp_decode_fn = Wrapper.WebPDecodeBGRA
    elseif TColor == RGB{N0f8}
        webp_decode_fn = Wrapper.WebPDecodeRGB
    elseif TColor == RGBA{N0f8}
        webp_decode_fn = Wrapper.WebPDecodeRGBA
    elseif TColor == Gray{N0f8}
        webp_decode_fn = Wrapper.WebPDecodeRGB
        TDecodedColor = RGB{N0f8}
    else
        throw(ArgumentError("Unsupported color type: $TColor"))
    end
    width = Ref{Int32}(-1)
    height = Ref{Int32}(-1)
    decoded_data_ptr = webp_decode_fn(pointer(data), length(data), width, height)
    decoded_data_size = (sizeof(TDecodedColor), Int(width[]), Int(height[]))
    decoded_data = unsafe_wrap(Array{UInt8, 3}, decoded_data_ptr, decoded_data_size)
    image_view = colorview(TDecodedColor, normedview(decoded_data))
    if TDecodedColor == TColor
        image = transpose ? collect(image_view) : permutedims(image_view, (2, 1))
    else
        image = if transpose
            TColor.(image_view)
        else
            TColor.(PermutedDimsArray(image_view, (2, 1)))
        end
    end
    Wrapper.WebPFree(decoded_data_ptr)
    return image
end

function decode(
    data::AbstractVector{UInt8}; kwargs...
)::Union{Matrix{RGB{N0f8}}, Matrix{RGBA{N0f8}}}
    bitstream_features = Ref{Wrapper.WebPBitstreamFeatures}()
    # WebPGetFeatures is not available in libwebp dynamic library, but WebPGetFeaturesInternal is equivalent: https://github.com/webmproject/libwebp/blob/v1.4.0/src/webp/decode.h#L441
    Wrapper.WebPGetFeaturesInternal(
        pointer(data), length(data), bitstream_features, Wrapper.WEBP_DECODER_ABI_VERSION
    )
    has_alpha = bitstream_features[].has_alpha != 0
    TColor = has_alpha ? RGBA{N0f8} : RGB{N0f8}
    return decode(TColor, data; kwargs...)
end

function read_webp(
    ::Type{CT}, f::Union{AbstractString, IO}; kwargs...
)::Matrix{CT} where {CT <: Colorant}
    return decode(CT, Base.read(f); kwargs...)
end

function read_webp(
    f::Union{AbstractString, IO}; kwargs...
)::Union{Matrix{RGB{N0f8}}, Matrix{RGBA{N0f8}}}
    return decode(Base.read(f); kwargs...)
end
