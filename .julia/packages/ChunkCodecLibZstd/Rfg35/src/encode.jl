"""
    struct ZstdEncodeOptions <: EncodeOptions
    ZstdEncodeOptions(; kwargs...)

Zstandard compression using libzstd: www.zstd.net

Zstandard's format is documented in [RFC8878](https://datatracker.ietf.org/doc/html/rfc8878)

# Keyword Arguments

- `codec::ZstdCodec=ZstdCodec()`
- `compressionLevel::Integer=0`: Compression level, regular levels are 1-22.

  Levels ≥ 20 should be used with caution, as they require more memory.
  The library also offers negative compression levels,
  which extend the range of speed vs. ratio preferences.
  The lower the level, the faster the speed (at the cost of compression).
  0 is a special value for `ZSTD_defaultCLevel()`.
  The level will be clamped to the range `ZSTD_minCLevel()` to `ZSTD_maxCLevel()`.
- `checksum::Bool=false`: A 32-bits checksum of content is written at end of frame.
- `advanced_parameters::Vector{Pair{Int32, Int32}}=[]`: Pairs of `param => value`.

  Warning, some parameters can result in encodings that are incompatible with default decoders.
  Some parameters are experimental and may change in new versions of libzstd,
  so you may need to check [`ZSTD_versionNumber`](@ref) and [`ZSTD_cParam_getBounds`](@ref).
  See comments in zstd.h.
  Additional parameters are set with `ZSTD_CCtx_setParameter`. These parameters
  are set after the compression level, and checksum options are set, 
  so they can override those values.
  The vector must not be mutated after calling the constructor.
"""
struct ZstdEncodeOptions <: EncodeOptions
    codec::ZstdCodec
    compressionLevel::Int32
    checksum::Bool
    advanced_parameters::Vector{Pair{Int32, Int32}}
end
function ZstdEncodeOptions(;
        codec::ZstdCodec=ZstdCodec(),
        compressionLevel::Integer=0,
        checksum::Bool=false,
        advanced_parameters::Vector{Pair{Int32, Int32}}=Pair{Int32, Int32}[],
        kwargs...
    )
    _clamped_compression_level = clamp(compressionLevel, ZSTD_minCLevel(), ZSTD_maxCLevel())
    ZstdEncodeOptions(codec, _clamped_compression_level, checksum, advanced_parameters)
end

is_thread_safe(::ZstdEncodeOptions) = true

function decoded_size_range(::ZstdEncodeOptions)
    # prevent overflow of encode_bound
    # like ZSTD_MAX_INPUT_SIZE for Int64
    # From ChunkCodecTests.find_max_decoded_size(ZstdEncodeOptions())
    Int64(0):Int64(1):Int64(0x7F807F807F807F7F)
end

function encode_bound(::ZstdEncodeOptions, src_size::Int64)::Int64
    # ZSTD_COMPRESSBOUND ported to Julia
    # This also works when streaming
    # assuming no extra flushes
    # https://github.com/facebook/zstd/issues/3935
    # From zstd.h
    # #define ZSTD_COMPRESSBOUND(srcSize) (((size_t)(srcSize) >= ZSTD_MAX_INPUT_SIZE) ? 0 : (srcSize) + ((srcSize)>>8) + (((srcSize) < (128<<10)) ? (((128<<10) - (srcSize)) >> 11) /* margin, from 64 to 0 */ : 0)) 
    # /* this formula ensures that bound(A) + bound(B) <= bound(A+B) as long as A and B >= 128 KB */

    # Here we use Int64 instead of size_t
    margin = if src_size < (Int64(128)<<10)
        (((Int64(128)<<10) - src_size) >> 11)
    else 
        Int64(0)
    end::Int64
    clamp(widen(src_size) + widen(src_size>>8 + margin), Int64)
end

function try_encode!(e::ZstdEncodeOptions, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}; kwargs...)::MaybeSize
    check_contiguous(dst)
    check_contiguous(src)
    src_size::Int64 = length(src)
    dst_size::Int64 = length(dst)
    check_in_range(decoded_size_range(e); src_size)
    # create ZSTD_CCtx
    cctx = ccall((:ZSTD_createCCtx, libzstd), Ptr{ZSTD_CCtx}, ())
    if cctx == C_NULL
        throw(OutOfMemoryError())
    end
    try
        # set parameters
        ZSTD_c_compressionLevel = Cint(100)
        ZSTD_c_checksumFlag = Cint(201)
        _set_parameter(cctx, ZSTD_c_compressionLevel, Cint(e.compressionLevel))
        if e.checksum
            _set_parameter(cctx, ZSTD_c_checksumFlag, Cint(1))
        end
        for (param, value) in e.advanced_parameters
            _set_parameter(cctx, Cint(param), Cint(value))
        end
        # do compression
        ret = ccall((:ZSTD_compress2, libzstd), Csize_t,
            (Ptr{ZSTD_CCtx}, Ref{UInt8}, Csize_t, Ref{UInt8}, Csize_t,),
            cctx, dst, dst_size, src, src_size,
        )
        if ZSTD_isError(ret)
            err_code = ZSTD_getErrorCode(ret)
            if err_code == Integer(ZSTD_error_dstSize_tooSmall)
                return NOT_SIZE
            elseif err_code == Integer(ZSTD_error_memory_allocation)
                throw(OutOfMemoryError())
            else
                error("unexpected libzstd error code $(err_code) from ZSTD_compress2.")
            end
        else
            return Int64(ret)
        end
    finally
        # free ZSTD_CCtx
        ccall((:ZSTD_freeCCtx, libzstd), Csize_t, (Ptr{ZSTD_CCtx},), cctx)
    end
end
