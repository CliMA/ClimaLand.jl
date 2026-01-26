# This file contains tests that require a large amount of memory (at least 15 GB)
# and take a long time to run. The tests are designed to check the 
# compression and decompression functionality of the ChunkCodecLibZlib package 
# with very large inputs. These tests are not run with CI

using ChunkCodecLibZlib:
    ChunkCodecLibZlib,
    GzipCodec,
    GzipEncodeOptions,
    ZlibEncodeOptions,
    ZlibCodec,
    encode,
    decode
using Test: @testset, @test

@testset "Big Memory Tests" begin
    Sys.WORD_SIZE == 64 || error("tests require 64 bit word size")
    @info "compressing zeros"
    let
        local n = 2^32-327671
        local u = decode(ZlibCodec(), encode(ZlibEncodeOptions(;level=0), zeros(UInt8, n)); size_hint=n, max_size=n)
        all_zero = all(iszero, u)
        len_n = length(u) == n
        @test all_zero && len_n
    end
    for n in (2^32 - 1, 2^32, 2^32 +1, 2^33)
        @info "compressing"
        local c = encode(GzipEncodeOptions(), zeros(UInt8, n))
        @info "decompressing"
        local u = decode(GzipCodec(), c; size_hint=n)
        c = nothing
        all_zero = all(iszero, u)
        len_n = length(u) == n
        @test all_zero && len_n
    end

    @info "compressing random"
    for n in (2^32 - 1, 2^32, 2^32 +1)
        local u = rand(UInt8, n)
        @info "compressing"
        local c = encode(GzipEncodeOptions(), u)
        @info "decompressing"
        local u2 = decode(GzipCodec(), c)
        c = nothing
        are_equal = u == u2
        @test are_equal
    end

    @info "decompressing huge concatenation"
    uncompressed = rand(UInt8, 2^20)
    @info "compressing"
    compressed = encode(GzipEncodeOptions(), uncompressed)
    total_compressed = UInt8[]
    sizehint!(total_compressed, length(compressed)*2^12)
    total_uncompressed = UInt8[]
    sizehint!(total_uncompressed, length(uncompressed)*2^12)
    for i in 1:2^12
        append!(total_uncompressed, uncompressed)
        append!(total_compressed, compressed)
    end
    @test length(total_compressed) > 2^32
    @info "decompressing"
    @test total_uncompressed == decode(GzipCodec(), total_compressed; size_hint=length(total_uncompressed))
end
