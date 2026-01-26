# This file contains tests that require a large amount of memory (at least 25 GB)
# and take a long time to run. The tests are designed to check the 
# compression and decompression functionality of the CodecBzip2 package 
# with very large inputs. These tests are not run with CI

using Test
using CodecBzip2

@testset "memory leak" begin
    function foo()
        for i in 1:1000000
            c = transcode(Bzip2Compressor(), zeros(UInt8,16))
            u = transcode(Bzip2Decompressor(), c)
        end
    end
    foo()
end

@testset "Big Memory Tests" begin
    Sys.WORD_SIZE == 64 || error("tests require 64 bit word size")
    @info "compressing zeros"
    for n in (2^32 - 1, 2^32, 2^32 +1)
        @info "compressing"
        local c = transcode(Bzip2Compressor, zeros(UInt8, n))
        @info "decompressing"
        local u = transcode(Bzip2Decompressor, c)
        c = nothing
        all_zero = all(iszero, u)
        len_n = length(u) == n
        @test all_zero && len_n
    end

    @info "compressing random"
    for n in (2^32 - 1, 2^32, 2^32 +1)
        local u = rand(UInt8, n)
        @info "compressing"
        local c = transcode(Bzip2Compressor, u)
        @info "decompressing"
        local u2 = transcode(Bzip2Decompressor, c)
        c = nothing
        are_equal = u == u2
        @test are_equal
    end

    @info "decompressing huge concatenation"
    uncompressed = rand(UInt8, 2^20)
    @info "compressing"
    compressed = transcode(Bzip2Compressor, uncompressed)
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
    @test total_uncompressed == transcode(Bzip2Decompressor, total_compressed)
end
