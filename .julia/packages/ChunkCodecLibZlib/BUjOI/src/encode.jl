"""
    struct ZlibEncodeOptions <: EncodeOptions
    ZlibEncodeOptions(; kwargs...)

zlib compression using libzlib: https://www.zlib.net/

This is the zlib format described in RFC 1950

# Keyword Arguments

- `codec::ZlibCodec=ZlibCodec()`
- `level::Integer=-1`: The compression level must be -1, or between 0 and 9.

  1 gives best speed, 9 gives best compression, 0 gives no compression at all
  (the input data is simply copied a block at a time). -1
  requests a default compromise between speed and compression (currently
  equivalent to level 6).
"""
struct ZlibEncodeOptions <: EncodeOptions
    codec::ZlibCodec
    level::Int32
end
function ZlibEncodeOptions(;
        codec::ZlibCodec=ZlibCodec(),
        level::Integer=-1,
        kwargs...
    )
    ZlibEncodeOptions(
        codec,
        Int32(clamp(level, -1, 9)),
    )
end

"""
    struct DeflateEncodeOptions <: EncodeOptions
    DeflateEncodeOptions(; kwargs...)

deflate compression using libzlib: https://www.zlib.net/

This is the deflate format described in RFC 1951

# Keyword Arguments

- `codec::DeflateCodec=DeflateCodec()`
- `level::Integer=-1`: The compression level must be -1, or between 0 and 9.

  1 gives best speed, 9 gives best compression, 0 gives no compression at all
  (the input data is simply copied a block at a time). -1
  requests a default compromise between speed and compression (currently
  equivalent to level 6).
"""
struct DeflateEncodeOptions <: EncodeOptions
    codec::DeflateCodec
    level::Int32
end
function DeflateEncodeOptions(;
        codec::DeflateCodec=DeflateCodec(),
        level::Integer=-1,
        kwargs...
    )
    DeflateEncodeOptions(
        codec,
        Int32(clamp(level, -1, 9)),
    )
end

"""
    struct GzipEncodeOptions <: EncodeOptions
    GzipEncodeOptions(; kwargs...)

gzip compression using libzlib: https://www.zlib.net/

This is the gzip (.gz) format described in RFC 1952

# Keyword Arguments

- `codec::GzipCodec=GzipCodec()`
- `level::Integer=-1`: The compression level must be -1, or between 0 and 9.

  1 gives best speed, 9 gives best compression, 0 gives no compression at all
  (the input data is simply copied a block at a time). -1
  requests a default compromise between speed and compression (currently
  equivalent to level 6).
"""
struct GzipEncodeOptions <: EncodeOptions
    codec::GzipCodec
    level::Int32
end
function GzipEncodeOptions(;
        codec::GzipCodec=GzipCodec(),
        level::Integer=-1,
        kwargs...
    )
    GzipEncodeOptions(
        codec,
        Int32(clamp(level, -1, 9)),
    )
end

const _AllEncodeOptions = Union{ZlibEncodeOptions, DeflateEncodeOptions, GzipEncodeOptions}

is_thread_safe(::_AllEncodeOptions) = true

# Modified from the deflateBound function in zlib/deflate.c
# Since the default gzip header, windowBits, and memLevel are always
# used in this package, the code can be simplified.
# If alternate settings are used this must be modified.
function encode_bound(e::_AllEncodeOptions, src_size::Int64)::Int64
    wraplen = _wraplen(e)
    clamp(
        widen(src_size) +
        widen(src_size>>12 + src_size>>14 + src_size>>25 + Int64(7) + wraplen),
        Int64,
    )
end
_wraplen(::ZlibEncodeOptions)    = Int64(6)
_wraplen(::DeflateEncodeOptions) = Int64(0)
_wraplen(::GzipEncodeOptions)    = Int64(18)

# max to prevent overflows in encode_bound
# From ChunkCodecTests.find_max_decoded_size(::EncodeOptions)
_max_decoded_size(::ZlibEncodeOptions)::Int64 = 0x7ff60087fa602c72
_max_decoded_size(::DeflateEncodeOptions)::Int64 = 0x7ff60087fa602c78
_max_decoded_size(::GzipEncodeOptions)::Int64 = 0x7ff60087fa602c66

decoded_size_range(e::_AllEncodeOptions) = Int64(0):Int64(1):_max_decoded_size(e)

function try_encode!(e::_AllEncodeOptions, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}; kwargs...)::MaybeSize
    # -15: deflate, 15: zlib, 15+16: gzip
    # smaller windowBits might break encode bound
    windowBits = _windowBits(e.codec)
    @assert windowBits ∈ (-15, 15, 15+16)
    check_contiguous(dst)
    check_contiguous(src)
    src_size::Int64 = length(src)
    dst_size::Int64 = length(dst)
    check_in_range(decoded_size_range(e); src_size)
    if iszero(dst_size)
        return NOT_SIZE
    end
    stream = ZStream()
    deflateInit2(stream, e.level, windowBits)
    try
        # deflate loop
        cconv_src = Base.cconvert(Ptr{UInt8}, src)
        cconv_dst = Base.cconvert(Ptr{UInt8}, dst)
        GC.@preserve cconv_src cconv_dst begin
            src_p = Base.unsafe_convert(Ptr{UInt8}, cconv_src)
            dst_p = Base.unsafe_convert(Ptr{UInt8}, cconv_dst)
            stream.next_in = src_p
            stream.next_out = dst_p
            src_left::Int64 = src_size
            dst_left::Int64 = dst_size
            while true
                start_avail_in = clamp(src_left, Cuint)
                start_avail_out = clamp(dst_left, Cuint)
                stream.avail_in = start_avail_in
                stream.avail_out = start_avail_out
                @assert !iszero(stream.avail_out)
                action = if stream.avail_in == src_left
                    Z_FINISH
                else
                    Z_NO_FLUSH
                end
                ret = ccall(
                    (:deflate, libz),
                    Cint,
                    (Ref{ZStream}, Cint),
                    stream, action,
                )
                @assert ret != Z_STREAM_ERROR # state not clobbered
                @assert stream.avail_in ≤ start_avail_in
                @assert stream.avail_out ≤ start_avail_out
                # there must be progress
                @assert stream.avail_in < start_avail_in || stream.avail_out < start_avail_out
                src_left -= start_avail_in - stream.avail_in
                dst_left -= start_avail_out - stream.avail_out
                @assert src_left ∈ 0:src_size
                @assert dst_left ∈ 0:dst_size
                @assert stream.next_in == src_p + (src_size - src_left)
                @assert stream.next_out == dst_p + (dst_size - dst_left)
                if ret == Z_STREAM_END
                    # yay done just return compressed size
                    @assert dst_left ∈ 0:dst_size
                    @assert iszero(src_left)
                    return dst_size - dst_left
                end
                if iszero(dst_left)
                    # no more space, but not Z_STREAM_END
                    return NOT_SIZE
                end
                @assert ret == Z_OK
            end
        end
    finally
        deflateEnd(stream)
    end
end
