# Compressor Codecs
# ==================

abstract type CompressorCodec <: TranscodingStreams.Codec end

function Base.show(io::IO, codec::CompressorCodec)
    print(io, summary(codec), "(level=$(codec.level), windowbits=$(codec.windowbits))")
end


# Gzip
# ----

struct GzipCompressor <: CompressorCodec
    zstream::ZStream
    level::Int
    windowbits::Int
end

"""
    GzipCompressor(;level=$(Z_DEFAULT_COMPRESSION), windowbits=$(Z_DEFAULT_WINDOWBITS))

Create a gzip compression codec.

Arguments
---------
- `level` (-1..9): compression level. 1 gives best speed, 9 gives best compression, 0 gives no compression at all (the input data is simply copied a block at a time). -1 requests a default compromise between speed and compression (currently equivalent to level 6).
- `windowbits` (9..15): size of history buffer is `2^windowbits`.

!!! warning
    `serialize` and `deepcopy` will not work with this codec due to stored raw pointers.
"""
function GzipCompressor(;level::Integer=Z_DEFAULT_COMPRESSION,
                         windowbits::Integer=Z_DEFAULT_WINDOWBITS)
    if !(-1 ≤ level ≤ 9)
        throw(ArgumentError("compression level must be within -1..9"))
    elseif !(9 ≤ windowbits ≤ 15)
        throw(ArgumentError("windowbits must be within 9..15"))
    end
    zstream = ZStream()
    finalizer(compress_finalizer!, zstream)
    return GzipCompressor(zstream, level, windowbits+16)
end

const GzipCompressorStream{S} = TranscodingStream{GzipCompressor,S} where S<:IO

"""
    GzipCompressorStream(stream::IO; kwargs...)

Create a gzip compression stream (see `GzipCompressor` for `kwargs`).

!!! warning
    `serialize` and `deepcopy` will not work with this stream due to stored raw pointers.
"""
function GzipCompressorStream(stream::IO; kwargs...)
    x, y = splitkwargs(kwargs, (:level, :windowbits))
    return TranscodingStream(GzipCompressor(;x...), stream; y...)
end


# Zlib
# ----

struct ZlibCompressor <: CompressorCodec
    zstream::ZStream
    level::Int
    windowbits::Int
end

"""
    ZlibCompressor(;level=$(Z_DEFAULT_COMPRESSION), windowbits=$(Z_DEFAULT_WINDOWBITS))

Create a zlib compression codec.

Arguments
---------
- `level` (-1..9): compression level. 1 gives best speed, 9 gives best compression, 0 gives no compression at all (the input data is simply copied a block at a time). -1 requests a default compromise between speed and compression (currently equivalent to level 6).
- `windowbits` (9..15): size of history buffer is `2^windowbits`.

!!! warning
    `serialize` and `deepcopy` will not work with this codec due to stored raw pointers.
"""
function ZlibCompressor(;level::Integer=Z_DEFAULT_COMPRESSION,
                         windowbits::Integer=Z_DEFAULT_WINDOWBITS)
    if !(-1 ≤ level ≤ 9)
        throw(ArgumentError("compression level must be within -1..9"))
    elseif !(9 ≤ windowbits ≤ 15)
        throw(ArgumentError("windowbits must be within 9..15"))
    end
    zstream = ZStream()
    finalizer(compress_finalizer!, zstream)
    return ZlibCompressor(zstream, level, windowbits)
end

const ZlibCompressorStream{S} = TranscodingStream{ZlibCompressor,S} where S<:IO

"""
    ZlibCompressorStream(stream::IO)

Create a zlib compression stream (see `ZlibCompressor` for `kwargs`).

!!! warning
    `serialize` and `deepcopy` will not work with this stream due to stored raw pointers.
"""
function ZlibCompressorStream(stream::IO; kwargs...)
    x, y = splitkwargs(kwargs, (:level, :windowbits))
    return TranscodingStream(ZlibCompressor(;x...), stream; y...)
end


# Deflate
# -------

struct DeflateCompressor <: CompressorCodec
    zstream::ZStream
    level::Int
    windowbits::Int
end

"""
    DeflateCompressor(;level=$(Z_DEFAULT_COMPRESSION), windowbits=$(Z_DEFAULT_COMPRESSION))

Create a deflate compression codec.

Arguments
---------
- `level` (-1..9): compression level. 1 gives best speed, 9 gives best compression, 0 gives no compression at all (the input data is simply copied a block at a time). -1 requests a default compromise between speed and compression (currently equivalent to level 6).
- `windowbits` (9..15): size of history buffer is `2^windowbits`.

!!! warning
    `serialize` and `deepcopy` will not work with this codec due to stored raw pointers.
"""
function DeflateCompressor(;level::Integer=Z_DEFAULT_COMPRESSION,
                        windowbits::Integer=Z_DEFAULT_WINDOWBITS)
    if !(-1 ≤ level ≤ 9)
        throw(ArgumentError("compression level must be within -1..9"))
    elseif !(9 ≤ windowbits ≤ 15)
        throw(ArgumentError("windowbits must be within 9..15"))
    end
    zstream = ZStream()
    finalizer(compress_finalizer!, zstream)
    return DeflateCompressor(zstream, level, -Int(windowbits))
end

const DeflateCompressorStream{S} = TranscodingStream{DeflateCompressor,S} where S<:IO

"""
    DeflateCompressorStream(stream::IO; kwargs...)

Create a deflate compression stream (see `DeflateCompressor` for `kwargs`).

!!! warning
    `serialize` and `deepcopy` will not work with this stream due to stored raw pointers.
"""
function DeflateCompressorStream(stream::IO; kwargs...)
    x, y = splitkwargs(kwargs, (:level, :windowbits))
    return TranscodingStream(DeflateCompressor(;x...), stream; y...)
end


# Methods
# -------

function TranscodingStreams.startproc(codec::CompressorCodec, state::Symbol, error_ref::Error)
    if codec.zstream.state == C_NULL
        code = deflate_init!(codec.zstream, codec.level, codec.windowbits)
        # errors in deflate_init! do not require clean up, so just throw
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
        code = deflate_reset!(codec.zstream)
        # errors in deflate_reset! do not require clean up, so just throw
        if code == Z_OK
            return :ok
        elseif code == Z_STREAM_ERROR
            error("Z_STREAM_ERROR: the source stream state was inconsistent")
        else
            error("unexpected libz error code: $(code)")
        end
    end
end

function TranscodingStreams.process(codec::CompressorCodec, input::Memory, output::Memory, error_ref::Error)
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
    code = deflate!(zstream, zstream.avail_in > 0 ? Z_NO_FLUSH : Z_FINISH)
    @assert code != Z_STREAM_ERROR # state not clobbered
    Δin = Int(avail_in - zstream.avail_in)
    Δout = Int(avail_out - zstream.avail_out)
    if code == Z_OK
        return Δin, Δout, :ok
    elseif code == Z_STREAM_END
        return Δin, Δout, :end
    else
        error_ref[] = ErrorException(zlib_error_message(zstream, code))
        return Δin, Δout, :error
    end
end
