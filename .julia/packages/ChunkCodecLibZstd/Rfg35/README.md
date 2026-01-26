# ChunkCodecLibZstd

This package implements the ChunkCodec interface for the following encoders and decoders
using the zstd C library <www.zstd.net>

1. `ZstdCodec`, `ZstdEncodeOptions`, `ZstdDecodeOptions`

## Example

```julia-repl
julia> using ChunkCodecLibZstd

julia> data = [0x00, 0x01, 0x02, 0x03];

julia> compressed_data = encode(ZstdEncodeOptions(;compressionLevel=3), data);

julia> decompressed_data = decode(ZstdCodec(), compressed_data; max_size=length(data), size_hint=length(data));

julia> data == decompressed_data
true
```

The low level interface is defined in the `ChunkCodecCore` package.

