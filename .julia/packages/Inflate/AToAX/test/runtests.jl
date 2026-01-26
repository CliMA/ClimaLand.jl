using Test
using Random
using Inflate
using CodecZlib: GzipCompressorStream, ZlibCompressorStream,
                 DeflateCompressorStream, CodecZlib

include("utils.jl")

empty_string = ""
short_string = "This is a short string."
medium_string = read(pathof(Inflate), String)
long_string = join(fill(medium_string, 1000), short_string)

@testset "Text strings" begin
    for s in [empty_string, short_string, medium_string, long_string]
        @test String(inflate(read(DeflateCompressorStream(IOBuffer(s))))) == s
        @test String(inflate_zlib(read(ZlibCompressorStream(IOBuffer(s))))) == s
        @test String(inflate_gzip(read(GzipCompressorStream(IOBuffer(s))))) == s
        @test read(InflateStream(DeflateCompressorStream(IOBuffer(s))), String) == s
        @test read(InflateZlibStream(ZlibCompressorStream(IOBuffer(s))), String) == s
        @test read(InflateGzipStream(GzipCompressorStream(IOBuffer(s))), String) == s
    end
end

@testset "Incompressible data" begin
    Random.seed!(1)
    for n in [0, 1, 10, 100, 1000, 10000, 100000, 1000000]
        data = rand(UInt8, n)
        @test inflate(read(DeflateCompressorStream(IOBuffer(data)))) == data
        @test inflate_zlib(read(ZlibCompressorStream(IOBuffer(data)))) == data
        @test inflate_gzip(read(GzipCompressorStream(IOBuffer(data)))) == data
        @test read(InflateStream(DeflateCompressorStream(IOBuffer(data)))) == data
        @test read(InflateZlibStream(ZlibCompressorStream(IOBuffer(data)))) == data
        @test read(InflateGzipStream(GzipCompressorStream(IOBuffer(data)))) == data
    end
end

@testset "Huffman compressible data" begin
    Random.seed!(1)
    for n in [0, 1, 10, 100, 1000, 10000, 100000, 1000000]
        data = rand(UInt8, n) .& 0x0f
        @test inflate(read(DeflateCompressorStream(IOBuffer(data)))) == data
        @test inflate_zlib(read(ZlibCompressorStream(IOBuffer(data)))) == data
        @test inflate_gzip(read(GzipCompressorStream(IOBuffer(data)))) == data
        @test read(InflateStream(DeflateCompressorStream(IOBuffer(data)))) == data
        @test read(InflateZlibStream(ZlibCompressorStream(IOBuffer(data)))) == data
        @test read(InflateGzipStream(GzipCompressorStream(IOBuffer(data)))) == data
    end
end

# Test gzip headers, including header CRC.
include("gzip_with_header.jl")
os = 255   # Unknown
mtime = 1537006040
fextra = UInt8[0x00, 0xbe, 0xef, 0x00]
fname = "foo.txt"
fcomment = "example string"
gz = gzip_with_header("foo", mtime, os, fextra, fname, fcomment, true)

@testset "Gzip headers" begin
    headers = Dict{String, Any}()
    @test String(inflate_gzip(gz, headers = headers)) == "foo"
    @test headers["mtime"] == mtime
    @test headers["os"] == os
    @test headers["fextra"] == fextra
    @test headers["fname"] == fname
    @test headers["fcomment"] == fcomment

    headers = Dict{String, Any}()
    @test read(InflateGzipStream(IOBuffer(gz), headers = headers), String) == "foo"
    @test headers["mtime"] == mtime
    @test headers["os"] == os
    @test headers["fextra"] == fextra
    @test headers["fname"] == fname
    @test headers["fcomment"] == fcomment
end

# Test readline interface and get coverage of read(stream, UInt8).
@testset "readline" begin
    s = "one line of text\n"
    @test readline(InflateGzipStream(GzipCompressorStream(IOBuffer(s))),
                   keep = true) == s
end

# Test gzip decompression of very large data (>2^32 zeros),
# specifically to stress test the length check.
#
# Note: These tests are disabled unless you manually change to `if
# true`. Since they inflate to 4 GB each they take some time to run
# but more importantly they would put undue memory strain on CI.
@testset "large gzip" begin
    if false
        @test all(==(0), inflate_gzip(all_zeros_gzip(2^32 + 1)))
        @test all(==(0), read(InflateGzipStream(IOBuffer(all_zeros_gzip(2^32 + 1)))))
    end
end

# Test failure cases, mostly corrupt data.
include("provoke_errors.jl")
