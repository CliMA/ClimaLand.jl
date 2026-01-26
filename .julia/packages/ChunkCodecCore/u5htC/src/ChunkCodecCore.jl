module ChunkCodecCore

export decode, encode

if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("""
        public
            Codec,
            EncodeOptions,
            DecodeOptions,

            MaybeSize,
            is_size,
            NOT_SIZE,
            DecodingError,
            DecodedSizeError,
            decode!,

            decode_options,

            decoded_size_range,
            encode_bound,
            try_encode!,

            try_find_decoded_size,
            try_decode!,

            check_in_range,
            check_contiguous,
            grow_dst!,

            can_concatenate,
            is_thread_safe,
            try_resize_decode!,
            is_lossless,

            NoopCodec,
            NoopEncodeOptions,
            NoopDecodeOptions,

            ShuffleCodec,
            ShuffleEncodeOptions,
            ShuffleDecodeOptions
    """))
end

include("types.jl")
include("errors.jl")
include("interface.jl")
include("noop.jl")
include("shuffle.jl")

end
