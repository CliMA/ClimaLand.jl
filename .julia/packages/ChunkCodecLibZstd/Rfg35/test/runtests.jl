using Random: Random
using ChunkCodecCore: encode_bound, decoded_size_range, encode, decode
using ChunkCodecLibZstd:
    ChunkCodecLibZstd,
    ZstdCodec,
    ZstdEncodeOptions,
    ZstdDecodeOptions,
    ZstdDecodingError,
    ZSTD_minCLevel,
    ZSTD_maxCLevel,
    ZSTD_defaultCLevel,
    ZSTD_versionNumber,
    ZSTD_isError,
    ZSTD_bounds,
    ZSTD_cParam_getBounds,
    ZSTD_dParam_getBounds
using ChunkCodecTests: test_codec
using Test: @testset, @test_throws, @test
using Aqua: Aqua

Aqua.test_all(ChunkCodecLibZstd; persistent_tasks = false)

Random.seed!(1234)

@testset "library info functions" begin
    ZSTD_c_compressionLevel = Int32(100)
    bounds = ZSTD_cParam_getBounds(ZSTD_c_compressionLevel)
    @test !ZSTD_isError(bounds.error)
    @test bounds.lowerBound == ZSTD_minCLevel()
    @test bounds.upperBound == ZSTD_maxCLevel()
    @test ZSTD_defaultCLevel() ∈ ZSTD_minCLevel():ZSTD_maxCLevel()
    @test ZSTD_versionNumber().major == 1
    ZSTD_d_windowLogMax = Int32(100)
    bounds = ZSTD_dParam_getBounds(ZSTD_d_windowLogMax)
    @test bounds isa ZSTD_bounds
end
@testset "encode_bound" begin
    local a = last(decoded_size_range(ZstdEncodeOptions()))
    @test encode_bound(ZstdEncodeOptions(), a) == typemax(Int64) - 1
    # zstd has adds a margin to encode bound for sizes less than 128 KB
    # Ensure this doesn't break monotonicity
    for i in 1:(Int64(128)<<10 + 100)
        @test encode_bound(ZstdEncodeOptions(), i) ≥ encode_bound(ZstdEncodeOptions(), i-1)
    end
end
@testset "default" begin
    test_codec(ZstdCodec(), ZstdEncodeOptions(), ZstdDecodeOptions(); trials=100)
end
@testset "compressionLevel options" begin
    # Compression level is clamped
    for i in [typemin(Int128); -30:25;]
        test_codec(ZstdCodec(), ZstdEncodeOptions(; compressionLevel=i), ZstdDecodeOptions(); trials=3)
    end
end
@testset "other options" begin
    # checksum
    test_codec(ZstdCodec(), ZstdEncodeOptions(; checksum=true), ZstdDecodeOptions(); trials=50)
    test_codec(ZstdCodec(), ZstdEncodeOptions(; checksum=false), ZstdDecodeOptions(); trials=50)
    # advanced parameters
    # As an example turn off the content size flag
    # From zstd.h:
    ZSTD_c_contentSizeFlag = Int32(200)
    # /* Content size will be written into frame header _whenever known_ (default:1)
    # * Content size must be known at the beginning of compression.
    # * This is automatically the case when using ZSTD_compress2(),
    # * For streaming scenarios, content size must be provided with ZSTD_CCtx_setPledgedSrcSize() */
    e = ZstdEncodeOptions(;advanced_parameters=[
        ZSTD_c_contentSizeFlag=>Int32(0),
    ])
    test_codec(ZstdCodec(), e, ZstdDecodeOptions(); trials=50)
end
@testset "unexpected eof" begin
    e = ZstdEncodeOptions()
    d = ZstdDecodeOptions()
    u = [0x00, 0x01, 0x02]
    c = encode(e, u)
    @test decode(d, c) == u
    for i in 1:length(c)
        @test_throws ZstdDecodingError decode(d, c[1:i-1])
    end
    @test_throws ZstdDecodingError decode(d, u)
    c[2] = 0x00
    @test_throws ZstdDecodingError decode(d, c)
    @test_throws ZstdDecodingError decode(d, [encode(e, u); c])
    @test_throws ZstdDecodingError decode(d, [encode(e, u); 0x00])
    e = ZstdEncodeOptions(;checksum=true)
    c = encode(e, u)
    c[end] = 0x00
    # This fails checksum
    @test_throws ZstdDecodingError decode(d, c)
end
@testset "errors" begin
    @test sprint(Base.showerror, ZstdDecodingError(:foo)) ==
        "ZstdDecodingError: foo"
    @test startswith(
        sprint(Base.showerror, ZstdDecodingError(-1%Csize_t)),
        "ZstdDecodingError: ",
    )
end
@testset "public" begin
    if VERSION >= v"1.11.0-DEV.469"
        for sym in (
            :ZSTD_minCLevel,
            :ZSTD_maxCLevel,
            :ZSTD_defaultCLevel,
            :ZSTD_versionNumber,
            :ZSTD_isError,
            :ZSTD_bounds,
            :ZSTD_cParam_getBounds,
            :ZSTD_dParam_getBounds
        )
            @test Base.ispublic(ChunkCodecLibZstd, sym)
        end
    end
end
