CodecZlib.jl
============

[![CI](https://github.com/JuliaIO/CodecZlib.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaIO/CodecZlib.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/JuliaIO/CodecZlib.jl/graph/badge.svg?token=6V3Z847Ywr)](https://codecov.io/gh/JuliaIO/CodecZlib.jl)

## Installation

```julia
Pkg.add("CodecZlib")
```

## Usage

```julia
using CodecZlib

# Some text.
text = """
Lorem ipsum dolor sit amet, consectetur adipiscing elit. Aenean sollicitudin
mauris non nisi consectetur, a dapibus urna pretium. Vestibulum non posuere
erat. Donec luctus a turpis eget aliquet. Cras tristique iaculis ex, eu
malesuada sem interdum sed. Vestibulum ante ipsum primis in faucibus orci luctus
et ultrices posuere cubilia Curae; Etiam volutpat, risus nec gravida ultricies,
erat ex bibendum ipsum, sed varius ipsum ipsum vitae dui.
"""

# Streaming API.
stream = GzipCompressorStream(IOBuffer(text))
for line in eachline(GzipDecompressorStream(stream))
    println(line)
end
close(stream)

# Array API.
compressed = transcode(GzipCompressor, text)
@assert sizeof(compressed) < sizeof(text)
@assert transcode(GzipDecompressor, compressed) == Vector{UInt8}(text)
```

This package exports following codecs and streams:

| Codec                  | Stream                       |
| ---------------------- | ---------------------------- |
| `GzipCompressor`       | `GzipCompressorStream`       |
| `GzipDecompressor`     | `GzipDecompressorStream`     |
| `ZlibCompressor`       | `ZlibCompressorStream`       |
| `ZlibDecompressor`     | `ZlibDecompressorStream`     |
| `DeflateCompressor`    | `DeflateCompressorStream`    |
| `DeflateDecompressor`  | `DeflateDecompressorStream`  |

See docstrings and [TranscodingStreams.jl](https://github.com/JuliaIO/TranscodingStreams.jl) for details.
