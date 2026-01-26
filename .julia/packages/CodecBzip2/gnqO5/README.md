# CodecBzip2.jl

[![CI](https://github.com/JuliaIO/CodecBzip2.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaIO/CodecBzip2.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/JuliaIO/CodecBzip2.jl/graph/badge.svg?token=eToD28jdA7)](https://codecov.io/gh/JuliaIO/CodecBzip2.jl)

## Installation

```julia
Pkg.add("CodecBzip2")
```

## Usage

```julia
using CodecBzip2

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
stream = Bzip2CompressorStream(IOBuffer(text))
for line in eachline(Bzip2DecompressorStream(stream))
    println(line)
end
close(stream)

# Array API.
compressed = transcode(Bzip2Compressor, text)
@assert sizeof(compressed) < sizeof(text)
@assert transcode(Bzip2Decompressor, compressed) == Vector{UInt8}(text)
```

This package exports following codecs and streams:

| Codec                | Stream                     |
| -------------------- | -------------------------- |
| `Bzip2Compressor`    | `Bzip2CompressorStream`    |
| `Bzip2Decompressor`  | `Bzip2DecompressorStream`  |

See docstrings and [TranscodingStreams.jl](https://github.com/JuliaIO/TranscodingStreams.jl) for details.
