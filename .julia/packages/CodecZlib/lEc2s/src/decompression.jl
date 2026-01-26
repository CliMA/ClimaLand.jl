# Decompressor Codecs
# ====================

abstract type DecompressorCodec <: TranscodingStreams.Codec end

function Base.show(io::IO, codec::DecompressorCodec)
    print(io, summary(codec), "(windowbits=$(codec.windowbits))")
end


# Gzip
# ----

struct GzipDecompressor <: DecompressorCodec
    zstream::ZStream
    windowbits::Int
end

"""
    GzipDecompressor(;windowbits=$(Z_DEFAULT_WINDOWBITS), gziponly=false)

Create a gzip decompressor codec.

If `gziponly` is `false`, this codec can decompress the zlib format as well.

Arguments
---------
- `windowbits` (8..15): Changing `windowbits` from its default of 15 will prevent decoding data using a history buffer larger than `2^windowbits`.
- `gziponly`: flag to inactivate data format detection

!!! warning
    `serialize` and `deepcopy` will not work with this codec due to stored raw pointers.
"""
function GzipDecompressor(;windowbits::Integer=Z_DEFAULT_WINDOWBITS, gziponly::Bool=false)
    if !(8 ≤ windowbits ≤ 15)
        throw(ArgumentError("windowbits must be within 8..15"))
    end
    zstream = ZStream()
    finalizer(decompress_finalizer!, zstream)
    return GzipDecompressor(zstream, windowbits+(gziponly ? 16 : 32))
end

const GzipDecompressorStream{S} = TranscodingStream{GzipDecompressor,S} where S<:IO

"""
    GzipDecompressorStream(stream::IO; kwargs...)

Create a gzip decompression stream (see `GzipDecompressor` for `kwargs`).

!!! warning
    `serialize` and `deepcopy` will not work with this stream due to stored raw pointers.
"""
function GzipDecompressorStream(stream::IO; kwargs...)
    x, y = splitkwargs(kwargs, (:windowbits, :gziponly))
    return TranscodingStream(GzipDecompressor(;x...), stream; y...)
end


# Zlib
# ----

struct ZlibDecompressor <: DecompressorCodec
    zstream::ZStream
    windowbits::Int
end

"""
    ZlibDecompressor(;windowbits=$(Z_DEFAULT_WINDOWBITS))

Create a zlib decompression codec.

Arguments
---------
- `windowbits` (8..15): Changing `windowbits` from its default of 15 will prevent decoding data using a history buffer larger than `2^windowbits`.

!!! warning
    `serialize` and `deepcopy` will not work with this codec due to stored raw pointers.
"""
function ZlibDecompressor(;windowbits::Integer=Z_DEFAULT_WINDOWBITS)
    if !(8 ≤ windowbits ≤ 15)
        throw(ArgumentError("windowbits must be within 8..15"))
    end
    zstream = ZStream()
    finalizer(decompress_finalizer!, zstream)
    return ZlibDecompressor(zstream, windowbits)
end

const ZlibDecompressorStream{S} = TranscodingStream{ZlibDecompressor,S} where S<:IO

"""
    ZlibDecompressorStream(stream::IO; kwargs...)

Create a deflate decompression stream (see `ZlibDecompressor` for `kwargs`).

!!! warning
    `serialize` and `deepcopy` will not work with this stream due to stored raw pointers.
"""
function ZlibDecompressorStream(stream::IO; kwargs...)
    x, y = splitkwargs(kwargs, (:windowbits,))
    return TranscodingStream(ZlibDecompressor(;x...), stream; y...)
end


# Deflate
# -------

struct DeflateDecompressor <: DecompressorCodec
    zstream::ZStream
    windowbits::Int
end

"""
    DeflateDecompressor(;windowbits=$(Z_DEFAULT_WINDOWBITS))

Create a deflate decompression codec.

Arguments
---------
- `windowbits` (8..15): Changing `windowbits` from its default of 15 will prevent decoding data using a history buffer larger than `2^windowbits`.

!!! warning
    `serialize` and `deepcopy` will not work with this codec due to stored raw pointers.
"""
function DeflateDecompressor(;windowbits::Integer=Z_DEFAULT_WINDOWBITS)
    if !(8 ≤ windowbits ≤ 15)
        throw(ArgumentError("windowbits must be within 8..15"))
    end
    zstream = ZStream()
    finalizer(decompress_finalizer!, zstream)
    return DeflateDecompressor(zstream, -Int(windowbits))
end

const DeflateDecompressorStream{S} = TranscodingStream{DeflateDecompressor,S} where S<:IO

"""
    DeflateDecompressorStream(stream::IO; kwargs...)

Create a deflate decompression stream (see `DeflateDecompressor` for `kwargs`).

!!! warning
    `serialize` and `deepcopy` will not work with this stream due to stored raw pointers.
"""
function DeflateDecompressorStream(stream::IO; kwargs...)
    x, y = splitkwargs(kwargs, (:windowbits,))
    return TranscodingStream(DeflateDecompressor(;x...), stream; y...)
end


# Methods
# -------

function TranscodingStreams.startproc(codec::DecompressorCodec, ::Symbol, error_ref::Error)
    # indicate that no input data is being provided for future zlib compat
    codec.zstream.next_in = C_NULL
    codec.zstream.avail_in = 0
    if codec.zstream.state == C_NULL
        code = inflate_init!(codec.zstream, codec.windowbits)
        # errors in inflate_init! do not require clean up, so just throw
        if code == Z_OK
            return :ok
        elseif code == Z_MEM_ERROR
            throw(OutOfMemoryError())
        elseif code == Z_STREAM_ERROR
            error("Z_STREAM_ERROR: invalid parameter, this should be caught in the codec constructor")
        elseif code == Z_VERSION_ERROR
            error("Z_VERSION_ERROR: zlib library version is incompatible")
        else
            error("unexpected libz error code: $(code)")
        end
    else
        code = inflate_reset!(codec.zstream)
        # errors in inflate_reset! do not require clean up, so just throw
        if code == Z_OK
            return :ok
        elseif code == Z_STREAM_ERROR
            error("Z_STREAM_ERROR: the source stream state was inconsistent")
        else
            error("unexpected libz error code: $(code)")
        end
    end
end

function TranscodingStreams.process(codec::DecompressorCodec, input::Memory, output::Memory, error_ref::Error)
    zstream = codec.zstream
    if zstream.state == C_NULL
        error("startproc must be called before process")
    end
    zstream.next_in = input.ptr

    avail_in = min(input.size, typemax(UInt32))
    zstream.avail_in = avail_in
    zstream.next_out = output.ptr
    avail_out = min(output.size, typemax(UInt32))
    zstream.avail_out = avail_out
    code = inflate!(zstream, Z_NO_FLUSH)
    @assert code != Z_STREAM_ERROR # state not clobbered
    Δin = Int(avail_in - zstream.avail_in)
    Δout = Int(avail_out - zstream.avail_out)
    if code == Z_OK
        return Δin, Δout, :ok
    elseif code == Z_STREAM_END
        return Δin, Δout, :end
    elseif code == Z_MEM_ERROR
        throw(OutOfMemoryError())
    elseif code == Z_BUF_ERROR && iszero(input.size)
        error_ref[] = ZlibError("the compressed stream may be truncated")
        return Δin, Δout, :error
    else
        error_ref[] = ZlibError(zlib_error_message(zstream, code))
        return Δin, Δout, :error
    end
end
