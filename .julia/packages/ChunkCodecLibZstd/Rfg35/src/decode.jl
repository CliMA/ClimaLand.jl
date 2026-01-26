"""
    ZstdDecodingError(x)

Error for data that cannot be decoded.
"""
struct ZstdDecodingError <: DecodingError
    x::Union{Symbol, Csize_t}
end

function Base.showerror(io::IO, err::ZstdDecodingError)
    print(io, "ZstdDecodingError: ")
    if err.x isa Symbol
        print(io, err.x)
    else
        print(io, ZSTD_getErrorName(err.x))
    end
    nothing
end

"""
    struct ZstdDecodeOptions <: DecodeOptions
    ZstdDecodeOptions(; kwargs...)

Zstandard decompression using libzstd: www.zstd.net

Zstandard's format is documented in [RFC8878](https://datatracker.ietf.org/doc/html/rfc8878)

# Keyword Arguments

- `codec::ZstdCodec=ZstdCodec()`
- `advanced_parameters::Vector{Pair{Int32, Int32}}=[]`: Pairs of `param => value`.

  Warning, some parameters are experimental and may change in new versions of libzstd,
  so you may need to check [`ZSTD_versionNumber`](@ref) and [`ZSTD_dParam_getBounds`](@ref).
  See comments in zstd.h.
  Additional parameters are set with `ZSTD_DCtx_setParameter`.
  The vector must not be mutated after calling the constructor.
"""
struct ZstdDecodeOptions <: DecodeOptions
    codec::ZstdCodec
    advanced_parameters::Vector{Pair{Int32, Int32}}
end
function ZstdDecodeOptions(;
        codec::ZstdCodec=ZstdCodec(),
        advanced_parameters::Vector{Pair{Int32, Int32}}=Pair{Int32, Int32}[],
        kwargs...
    )
    ZstdDecodeOptions(codec, advanced_parameters)
end

is_thread_safe(::ZstdDecodeOptions) = true

can_concatenate(::ZstdDecodeOptions) = true

# find_decompressed_size is modified from CodecZstd.jl
# https://github.com/JuliaIO/CodecZstd.jl/blob/2f7d084b8b157d83ed85e9d15105f0a708038e45/src/libzstd.jl#L157C1-L215C4
# From mkitti's PR https://github.com/JuliaIO/CodecZstd.jl/pull/63
function try_find_decoded_size(::ZstdDecodeOptions, src::AbstractVector{UInt8})::Union{Nothing, Int64}
    check_contiguous(src)
    srcSize::Int64 = length(src)
    frameOffset::Int64 = 0
    decompressedSize::Int64 = 0
    while frameOffset < srcSize
        remainingSize = srcSize - frameOffset
        # Obtain the decompressed frame content size of the next frame, accumulate
        frameContentSize = ccall(
            (:ZSTD_getFrameContentSize, libzstd), Culonglong,
            (Ref{UInt8}, Csize_t),
            Ref(src, firstindex(src) + frameOffset), remainingSize,
        )
        if frameContentSize == ZSTD_CONTENTSIZE_UNKNOWN
            return nothing
        end
        if frameContentSize > typemax(Int64) # also handles ZSTD_CONTENTSIZE_ERROR
            throw(ZstdDecodingError(:decoded_size_error))
        end
        decompressedSize, overflow = Base.Checked.add_with_overflow(decompressedSize, frameContentSize%Int64)
        if overflow
            throw(ZstdDecodingError(:decoded_size_overflow))
        end
        # Advance the offset forward by the size of the compressed frame
        # this is required if there are more than on frame
        ret = ccall(
            (:ZSTD_findFrameCompressedSize, libzstd), Csize_t,
            (Ref{UInt8}, Csize_t),
            Ref(src, firstindex(src) + frameOffset), remainingSize,
        )
        if ZSTD_isError(ret)
            err_code = ZSTD_getErrorCode(ret)
            if err_code == Integer(ZSTD_error_memory_allocation)
                throw(OutOfMemoryError())
            else
                throw(ZstdDecodingError(ret))
            end
        end
        @assert ret ∈ 1:remainingSize
        frameOffset += Int64(ret)
    end
    @assert frameOffset == srcSize
    return decompressedSize
end

function try_decode!(d::ZstdDecodeOptions, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}; kwargs...)::MaybeSize
    check_contiguous(dst)
    check_contiguous(src)
    if isempty(src)
        throw(ZstdDecodingError(:src_empty))
    end
    # create ZSTD_DCtx
    dctx = ccall((:ZSTD_createDCtx, libzstd), Ptr{ZSTD_DCtx}, ())
    if dctx == C_NULL
        throw(OutOfMemoryError())
    end
    try
        # set parameters
        for (param, value) in d.advanced_parameters
            _set_parameter(dctx, Cint(param), Cint(value))
        end
        # do decompression
        ret = ccall((:ZSTD_decompressDCtx, libzstd), Csize_t,
            (Ptr{ZSTD_DCtx}, Ref{UInt8}, Csize_t, Ref{UInt8}, Csize_t,),
            dctx, dst, length(dst), src, length(src),
        )
        if ZSTD_isError(ret)
            err_code = ZSTD_getErrorCode(ret)
            if err_code == Integer(ZSTD_error_dstSize_tooSmall)
                return NOT_SIZE
            elseif err_code == Integer(ZSTD_error_memory_allocation)
                throw(OutOfMemoryError())
            else
                throw(ZstdDecodingError(ret))
            end
        else
            @assert ret ∈ 0:length(dst)
            return Int64(ret)
        end
    finally
        # free ZSTD_DCtx
        ccall((:ZSTD_freeDCtx, libzstd), Csize_t, (Ptr{ZSTD_DCtx},), dctx)
    end
end

# For now rely on fallback `try_resize_decode!`
# Incase `try_find_decoded_size` returns `nothing`, the fallback repeatedly 
# calls `try_decode!` with larger and larger `dst`.
# This isn't ideal, but in a chunk decoding context
# the decoded size is typically found by `try_find_decoded_size`.
