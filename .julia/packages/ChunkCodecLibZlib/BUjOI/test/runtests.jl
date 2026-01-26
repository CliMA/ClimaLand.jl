using Random: Random
using ChunkCodecCore: encode_bound, decoded_size_range, encode, decode
using ChunkCodecLibZlib:
    ChunkCodecLibZlib,
    ZlibCodec,
    ZlibEncodeOptions,
    ZlibDecodeOptions,
    DeflateCodec,
    DeflateEncodeOptions,
    DeflateDecodeOptions,
    GzipCodec,
    GzipEncodeOptions,
    GzipDecodeOptions,
    LibzDecodingError
using ChunkCodecTests: test_codec
using Test: @testset, @test_throws, @test
using Aqua: Aqua

Aqua.test_all(ChunkCodecLibZlib; persistent_tasks = false)

Random.seed!(1234)

tests = [
    (
        ZlibCodec,
        ZlibEncodeOptions,
        ZlibDecodeOptions,
    ),
    (
        DeflateCodec,
        DeflateEncodeOptions,
        DeflateDecodeOptions,
    ),
    (
        GzipCodec,
        GzipEncodeOptions,
        GzipDecodeOptions,
    ),
]

@testset "tests for $(codec)" for (codec, encode_opt, decode_opt) in tests
    @testset "encode_bound" begin
        local a = last(decoded_size_range(encode_opt()))
        @test encode_bound(encode_opt(), a) == typemax(Int64) - 1
    end
    @testset "default" begin
        test_codec(codec(), encode_opt(), decode_opt(); trials=50)
    end
    @testset "level options" begin
        # level should get clamped to -1 to 9
        @test encode_opt(; level=-10).level == -1
        @test encode_opt(; level=10).level == 9
        @test encode_opt(; level=-2).level == -1
        for i in -1:9
            test_codec(codec(), encode_opt(; level=i), decode_opt(); trials=5)
        end
    end
    @testset "unexpected eof" begin
        local d = decode_opt()
        local u = [0x00, 0x01, 0x02]
        local c = encode(encode_opt(), u)
        @test decode(d, c) == u
        for i in 1:length(c)
            @test_throws LibzDecodingError("unexpected end of stream") decode(d, c[1:i-1])
        end
    end
end
@testset "unexpected eof with concatenation" begin
    e = GzipEncodeOptions()
    d = GzipDecodeOptions()
    u = [0x00, 0x01, 0x02]
    c = encode(e, u)
    @test decode(d, c) == u
    @test_throws LibzDecodingError decode(d, u)
    c[end] ⊻= 0xFF
    @test_throws LibzDecodingError decode(d, c)
    @test_throws LibzDecodingError decode(d, [encode(e, u); c])
    @test_throws LibzDecodingError decode(d, [encode(e, u); 0x00])
end
@testset "$(codec) doesn't support concatenation" for (codec, encode_opt, decode_opt) in tests[1:2]
    e = encode_opt()
    d = decode_opt()
    u = [0x00, 0x01, 0x02]
    c = encode(e, u)
    @test decode(d, c) == u
    @test_throws LibzDecodingError decode(d, u)
    c[begin] ⊻= 0xFF
    @test_throws LibzDecodingError decode(d, c)
    @test_throws LibzDecodingError("unexpected $(length(c)) bytes after stream") decode(d, [encode(e, u); c])
    @test_throws LibzDecodingError("unexpected $(length(c)) bytes after stream") decode(d, [encode(e, u); encode(e, u)])
    @test_throws LibzDecodingError("unexpected 1 bytes after stream") decode(d, [encode(e, u); 0x00])
end
@testset "Z_NEED_DICT error" begin
    @test_throws(
        LibzDecodingError("Z_NEED_DICT: a preset dictionary is needed at this point"),
        decode(
            ZlibDecodeOptions(),
            UInt8[0x78, 0xbb, 0x00, 0x00, 0x00, 0x01, 0x03, 0x00, 0x00, 0x00, 0x00, 0x01],
        )
    )
end
@testset "errors" begin
    @test sprint(Base.showerror, LibzDecodingError("test error message")) ==
        "LibzDecodingError: test error message"
end
