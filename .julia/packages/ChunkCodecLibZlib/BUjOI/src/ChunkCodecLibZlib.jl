module ChunkCodecLibZlib

using Zlib_jll: libz

using ChunkCodecCore:
    Codec,
    EncodeOptions,
    DecodeOptions,
    check_contiguous,
    check_in_range,
    grow_dst!,
    DecodingError,
    MaybeSize,
    NOT_SIZE
import ChunkCodecCore:
    decode_options,
    can_concatenate,
    try_decode!,
    try_resize_decode!,
    try_encode!,
    encode_bound,
    is_thread_safe,
    try_find_decoded_size,
    decoded_size_range

export ZlibCodec,
    DeflateCodec,
    GzipCodec,
    ZlibEncodeOptions,
    DeflateEncodeOptions,
    GzipEncodeOptions,
    ZlibDecodeOptions,
    DeflateDecodeOptions,
    GzipDecodeOptions,
    LibzDecodingError

# reexport ChunkCodecCore
using ChunkCodecCore: ChunkCodecCore, encode, decode
export ChunkCodecCore, encode, decode


include("libz.jl")

"""
    struct ZlibCodec <: Codec
    ZlibCodec()

zlib compression using libzlib: https://www.zlib.net/

This is the zlib format described in RFC 1950

See also [`ZlibEncodeOptions`](@ref) and [`ZlibDecodeOptions`](@ref)
"""
struct ZlibCodec <: Codec
end
decode_options(::ZlibCodec) = ZlibDecodeOptions()

# windowBits setting for the codec
_windowBits(::ZlibCodec) = Cint(15)

"""
    struct DeflateCodec <: Codec
    DeflateCodec()

deflate compression using libzlib: https://www.zlib.net/

This is the deflate format described in RFC 1951

See also [`DeflateEncodeOptions`](@ref) and [`DeflateDecodeOptions`](@ref)
"""
struct DeflateCodec <: Codec
end
decode_options(::DeflateCodec) = DeflateDecodeOptions()

# windowBits setting for the codec
_windowBits(::DeflateCodec) = Cint(-15)

"""
    struct GzipCodec <: Codec
    GzipCodec()

gzip compression using libzlib: https://www.zlib.net/

This is the gzip (.gz) format described in RFC 1952

See also [`GzipEncodeOptions`](@ref) and [`GzipDecodeOptions`](@ref)
"""
struct GzipCodec <: Codec
end
decode_options(::GzipCodec) = GzipDecodeOptions()

# windowBits setting for the codec
_windowBits(::GzipCodec) = Cint(15+16)

include("encode.jl")
include("decode.jl")

end # module ChunkCodecLibZlib
