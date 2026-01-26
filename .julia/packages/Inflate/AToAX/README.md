# Inflate.jl

Inflate provides a pure Julia implementation of
[zlib](https://zlib.net) *de*compression functionality, with both in-
memory and streaming interfaces. This covers decompression of the
Deflate algorithm and the Zlib and Gzip wrapper formats, as specified
in [RFC 1950](https://www.ietf.org/rfc/rfc1950.txt),
[RFC 1951](https://www.ietf.org/rfc/rfc1951.txt), and
[RFC 1952](https://www.ietf.org/rfc/rfc1952.txt).

The main reasons to choose Inflate over
[CodecZlib](https://github.com/bicycle1885/CodecZlib.jl) are:
* 100% Julia code - great for Julia purists.
* No binary dependencies.
* Actually no dependencies at all.
* Can read gzip headers.

You should choose CodecZlib over Inflate if the points above are not
compelling or one or more of the following applies to you:
* Need to compress, not only decompress.
* Want higher speed.
* Want a full-featured streaming interface.
* Want a battle-proven library.

## In-Memory Decompression

In-memory decompression is done by the following functions:

| function | decompresses |
| -------- | ------------ |
| `inflate(data::Vector{UInt8})` | raw Deflate data |
| `inflate_zlib(data::Vector{UInt8})` | Zlib data |
| `inflate_gzip(data::Vector{UInt8})` | Gzip data |

They all take a `Vector{UInt8}` with compressed data as input and
return a `Vector{UInt8}` of decompressed data. Additionally
```
gzip_headers = Dict{String, Any}()
out = inflate_gzip(data, headers = gzip_headers)
```
fills in `gzip_headers` with the Gzip headers present in `data`.

Both `inflate_zlib` and `inflate_gzip` accept the keyword argument
`ignore_checksum`, which if set to true skips consistency checking by
means of Adler and CRC checksums respectively. This disables the
computation of the checksums, saving time.

Finally,
there is also a convenience function to read a compressed text file in
gzip format
```
out = inflate_gzip(filename::String)
```
This returns the decompressed file as a string.


## Streaming Decompression

Streaming decompression is done using the following types:

| stream | decompresses |
| ------ | ------------ |
| `InflateStream(stream::IO)` | raw Deflate stream |
| `InflateZlibStream(stream::IO)` | Zlib stream |
| `InflateGzipStream(stream::IO)` | Gzip stream |

The stream types are subtypes of `IO` and decompression is done by
reading from instances of the types.

Example:
```
f = open("compressed_file.gz", "r")
gz = InflateGzipStream(f)
for line in readlines(gz)
    println(line)
end
close(f)
```
The streaming interface is minimalistic. If you need a full-featured
interface, the CodecZlib package is likely to be a better fit.

Reading of Gzip headers can be done from the streaming interface too.
```
gzip_headers = Dict{String, Any}()
gz = InflateGzipStream(stream, headers = gzip_headers)
```
The retrieved headers will be available immediately upon construction
of the `InflateGzipStream`. It is not necessary to read any data
first.

Likewise both `InflateZlibStream` and `InflateGzipStream` accept the
keyword argument `ignore_checksum` in the same way as the
non-streaming functions `inflate_zlib` and `inflate_gzip`.
