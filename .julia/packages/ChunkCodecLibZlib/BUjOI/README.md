# ChunkCodecLibZlib

This package implements the ChunkCodec interface for the following encoders and decoders
using the zlib C library <https://www.zlib.net/>

1. `ZlibCodec`, `ZlibEncodeOptions`, `ZlibDecodeOptions`
2. `DeflateCodec`, `DeflateEncodeOptions`, `DeflateDecodeOptions`
3. `GzipCodec`, `GzipEncodeOptions`, `GzipDecodeOptions`

## Example

```julia-repl
julia> using ChunkCodecLibZlib

julia> data = [0x00, 0x01, 0x02, 0x03];

julia> compressed_data = encode(GzipEncodeOptions(;level=6), data);

julia> decompressed_data = decode(GzipCodec(), compressed_data; max_size=length(data), size_hint=length(data));

julia> data == decompressed_data
true
```

The low level interface is defined in the `ChunkCodecCore` package.

