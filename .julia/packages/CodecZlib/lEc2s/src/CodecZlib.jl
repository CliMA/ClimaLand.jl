module CodecZlib

export
    # gzip
    GzipCompressor,
    GzipCompressorStream,
    GzipDecompressor,
    GzipDecompressorStream,

    # zlib
    ZlibCompressor,
    ZlibCompressorStream,
    ZlibDecompressor,
    ZlibDecompressorStream,

    # deflate
    DeflateCompressor,
    DeflateCompressorStream,
    DeflateDecompressor,
    DeflateDecompressorStream

import TranscodingStreams:
    TranscodingStreams,
    TranscodingStream,
    Memory,
    Error,
    initialize,
    finalize,
    splitkwargs
using Zlib_jll: libz

include("libz.jl")
include("compression.jl")
include("decompression.jl")

end # module
