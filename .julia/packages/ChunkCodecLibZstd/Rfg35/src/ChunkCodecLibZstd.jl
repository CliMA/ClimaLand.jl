module ChunkCodecLibZstd

using Zstd_jll: libzstd

using ChunkCodecCore:
    Codec,
    EncodeOptions,
    DecodeOptions,
    check_contiguous,
    check_in_range,
    DecodingError,
    MaybeSize,
    NOT_SIZE
import ChunkCodecCore:
    can_concatenate,
    try_decode!,
    try_encode!,
    encode_bound,
    is_thread_safe,
    try_find_decoded_size,
    decoded_size_range,
    decode_options

export ZstdCodec,
    ZstdEncodeOptions,
    ZstdDecodeOptions,
    ZstdDecodingError


if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("""
        public
            ZSTD_minCLevel,
            ZSTD_maxCLevel,
            ZSTD_defaultCLevel,
            ZSTD_versionNumber,
            ZSTD_isError,
            ZSTD_bounds,
            ZSTD_cParam_getBounds,
            ZSTD_dParam_getBounds
    """))
end

# reexport ChunkCodecCore
using ChunkCodecCore: ChunkCodecCore, encode, decode
export ChunkCodecCore, encode, decode


include("libzstd.jl")

"""
    struct ZstdCodec <: Codec
    ZstdCodec()

Zstandard compression using libzstd: www.zstd.net

Zstandard's format is documented in [RFC8878](https://datatracker.ietf.org/doc/html/rfc8878)

Like libzstd's simple API, encode compresses data as a single frame with saved
decompressed size. Decoding will succeed even if the decompressed size is unknown.
Also like libzstd's simple API, decoding accepts concatenated frames 
and will error if there is invalid data appended.

See also [`ZstdEncodeOptions`](@ref) and [`ZstdDecodeOptions`](@ref)
"""
struct ZstdCodec <: Codec
end
decode_options(::ZstdCodec) = ZstdDecodeOptions()

include("encode.jl")
include("decode.jl")

end # module ChunkCodecLibZstd
