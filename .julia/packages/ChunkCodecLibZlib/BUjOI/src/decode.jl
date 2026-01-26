"""
    LibzDecodingError(msg)

Error for data that cannot be decoded.
"""
struct LibzDecodingError <: DecodingError
    msg::String
end

function Base.showerror(io::IO, err::LibzDecodingError)
    print(io, "LibzDecodingError: ")
    print(io, err.msg)
    nothing
end

"""
    struct ZlibDecodeOptions <: DecodeOptions
    ZlibDecodeOptions(; kwargs...)

zlib decompression using libzlib: https://www.zlib.net/

This is the zlib format described in RFC 1950

# Keyword Arguments

- `codec::ZlibCodec=ZlibCodec()`
"""
struct ZlibDecodeOptions <: DecodeOptions
    codec::ZlibCodec
end
function ZlibDecodeOptions(;
        codec::ZlibCodec=ZlibCodec(),
        kwargs...
    )
    ZlibDecodeOptions(codec)
end

"""
    struct DeflateDecodeOptions <: DecodeOptions
    DeflateDecodeOptions(; kwargs...)

deflate decompression using libzlib: https://www.zlib.net/

This is the deflate format described in RFC 1951

# Keyword Arguments

- `codec::DeflateCodec=DeflateCodec()`
"""
struct DeflateDecodeOptions <: DecodeOptions
    codec::DeflateCodec
end
function DeflateDecodeOptions(;
        codec::DeflateCodec=DeflateCodec(),
        kwargs...
    )
    DeflateDecodeOptions(codec)
end


"""
    struct GzipDecodeOptions <: DecodeOptions
    GzipDecodeOptions(; kwargs...)

gzip decompression using libzlib: https://www.zlib.net/

This is the gzip (.gz) format described in RFC 1952

# Keyword Arguments

- `codec::GzipCodec=GzipCodec()`
"""
struct GzipDecodeOptions <: DecodeOptions
    codec::GzipCodec
end
function GzipDecodeOptions(;
        codec::GzipCodec=GzipCodec(),
        kwargs...
    )
    GzipDecodeOptions(codec)
end

can_concatenate(::GzipDecodeOptions) = true

const _AllDecodeOptions = Union{ZlibDecodeOptions, DeflateDecodeOptions, GzipDecodeOptions}

is_thread_safe(::_AllDecodeOptions) = true

function try_find_decoded_size(::_AllDecodeOptions, src::AbstractVector{UInt8})::Nothing
    nothing
end
function try_decode!(d::_AllDecodeOptions, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}; kwargs...)::MaybeSize
    try_resize_decode!(d, dst, src, Int64(length(dst)))
end

function try_resize_decode!(d::_AllDecodeOptions, dst::AbstractVector{UInt8}, src::AbstractVector{UInt8}, max_size::Int64; kwargs...)::MaybeSize
    dst_size::Int64 = length(dst)
    src_size::Int64 = length(src)
    src_left::Int64 = src_size
    dst_left::Int64 = dst_size
    check_contiguous(dst)
    check_contiguous(src)
    if isempty(src)
        throw(LibzDecodingError("unexpected end of stream"))
    end
    cconv_src = Base.cconvert(Ptr{UInt8}, src)
    # This outer loop is to decode a concatenation of multiple compressed streams.
    # If `can_concatenate(d)` is false, this outer loop doesn't rerun.
    while true
        stream = ZStream()
        inflateInit2(stream, _windowBits(d.codec))
        try
            # This inner loop is needed because libz can work on at most 
            # 2^32 - 1 bytes at a time.
            while true
                # dst may get resized, so cconvert needs to be redone on each iteration.
                cconv_dst = Base.cconvert(Ptr{UInt8}, dst)
                GC.@preserve cconv_src cconv_dst begin
                    src_p = Base.unsafe_convert(Ptr{UInt8}, cconv_src)
                    dst_p = Base.unsafe_convert(Ptr{UInt8}, cconv_dst)
                    start_avail_in = clamp(src_left, Cuint)
                    start_avail_out = clamp(dst_left, Cuint)
                    stream.avail_in = start_avail_in
                    stream.avail_out = start_avail_out
                    stream.next_in = src_p + (src_size - src_left)
                    stream.next_out = dst_p + (dst_size - dst_left)
                    # @info (;stream.avail_in, stream.avail_out)
                    ret = ccall(
                        (:inflate, libz),
                        Cint,
                        (Ref{ZStream}, Cint),
                        stream, Z_NO_FLUSH,
                    )
                    # @info (;ret, stream.avail_in, stream.avail_out)
                    if ret == Z_OK || ret == Z_STREAM_END || ret == Z_BUF_ERROR
                        @assert stream.avail_in ≤ start_avail_in
                        @assert stream.avail_out ≤ start_avail_out
                        dst_left -= start_avail_out - stream.avail_out
                        src_left -= start_avail_in - stream.avail_in
                        @assert src_left ∈ 0:src_size
                        @assert dst_left ∈ 0:dst_size
                    end
                    if ret == Z_OK || ret == Z_BUF_ERROR
                        if (
                                iszero(stream.avail_out) && !iszero(dst_left) ||
                                iszero(stream.avail_in) && !iszero(src_left) ||
                                !iszero(src_left) && !iszero(dst_left)
                            )
                            # there must be progress
                            @assert ret != Z_BUF_ERROR
                            @assert stream.avail_in < start_avail_in || stream.avail_out < start_avail_out
                        elseif iszero(dst_left) # needs more output
                            local next_size = grow_dst!(dst, max_size)
                            if isnothing(next_size)
                                return NOT_SIZE
                            end
                            dst_left += next_size - dst_size
                            dst_size = next_size
                            @assert dst_left > 0
                        else
                            @assert iszero(src_left) && !iszero(stream.avail_out)
                            throw(LibzDecodingError("unexpected end of stream"))
                        end
                    elseif ret == Z_STREAM_END
                        if iszero(src_left)
                            # yay done return decompressed size
                            real_dst_size = dst_size - dst_left
                            @assert real_dst_size ∈ 0:length(dst)
                            return real_dst_size
                        else
                            if can_concatenate(d)
                                # try and decompress next stream if the decoder can_concatenate
                                # there must be progress
                                @assert stream.avail_in < start_avail_in || stream.avail_out < start_avail_out
                                break
                            else
                                # Otherwise, throw an error
                                throw(LibzDecodingError("unexpected $(src_left) bytes after stream"))
                            end
                        end
                    elseif ret == Z_STREAM_ERROR
                        error("Z_STREAM_ERROR this should be unreachable")
                    elseif ret == Z_NEED_DICT
                        throw(LibzDecodingError("Z_NEED_DICT: a preset dictionary is needed at this point"))
                    elseif ret == Z_DATA_ERROR
                        if stream.msg == C_NULL
                            throw(LibzDecodingError("unknown Z_DATA_ERROR"))
                        else
                            throw(LibzDecodingError(unsafe_string(stream.msg)))
                        end
                    elseif ret == Z_MEM_ERROR
                        throw(OutOfMemoryError())
                    else
                        error("unknown libz error code: $(ret)")
                    end
                end
            end
        finally
            inflateEnd(stream)
        end
    end
end
