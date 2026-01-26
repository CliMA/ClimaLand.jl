# Deflate compression of empty string.
empty_deflate = [0x03, 0x00]

# Zlib compression of empty string.
empty_zlib = [0x78, 0x9c, 0x03, 0x00, 0x00, 0x00, 0x00, 0x01]

# Gzip compression of empty string with no headers other than header crc.
empty_gzip = [0x1f, 0x8b, 0x08, 0x03, 0x00, 0x00, 0x00, 0x00,
              0x02, 0x03, 0x91, 0x1e, 0x03, 0x00, 0x00, 0x00,
              0x00, 0x00, 0x00, 0x00, 0x00, 0x00]

@testset "Empty messages" begin
    @test inflate(empty_deflate) == UInt8[]
    @test inflate_zlib(empty_zlib) == UInt8[]
    @test inflate_gzip(empty_gzip) == UInt8[]
end

@testset "Deflate corruption" begin
    d1 = UInt8[0x01, 0x00, 0x00, 0x00, 0x00] # corrupted compression mode 0 data
    d2 = UInt8[0x07]                         # invalid compression mode 3
    d3 = UInt8[0xed, 0x1c, 0xed, 0x72,
               0xdb, 0x48, 0xf2, 0x3f]       # incomplete code table
    for d in [d1, d2, d3]
        @test_throws ErrorException inflate(d)
        @test_throws ErrorException read(InflateStream(IOBuffer(d)))
    end
end

@testset "Zlib header corruption" begin
    z1 = copy(empty_zlib); z1[1] = 0x77   # unsupported compression method
    z2 = copy(empty_zlib); z2[1] = 0x88   # invalid LZ77 window size
    z3 = copy(empty_zlib); z3[2] = 0xbc   # preset dictionary not supported
    z4 = copy(empty_zlib); z4[2] = 0x9d   # header checksum error
    for z in [z1, z2, z3, z4]
        @test_throws ErrorException inflate_zlib(z)
        @test_throws ErrorException InflateZlibStream(IOBuffer(z))
    end
end

@testset "Zlib trailer corruption" begin
    z5 = copy(empty_zlib); z5[5] = 0x01   # adler checksum error
    z6 = copy(empty_zlib); z6[6] = 0x01   # adler checksum error
    z7 = copy(empty_zlib); z7[7] = 0x01   # adler checksum error
    z8 = copy(empty_zlib); z8[8] = 0x02   # adler checksum error
    for z in [z5, z6, z7, z8]
        @test_throws ErrorException inflate_zlib(z)
        @test_throws ErrorException read(InflateZlibStream(IOBuffer(z)))
        @test isempty(inflate_zlib(z, ignore_checksum = true))
        @test eof(InflateZlibStream(IOBuffer(z), ignore_checksum = true))
    end
end

@testset "Gzip header corruption" begin
    g1 = copy(empty_gzip); g1[1] = 0x1e   # not gzipped data
    g2 = copy(empty_gzip); g2[2] = 0x8a   # not gzipped data
    g3 = copy(empty_gzip); g3[3] = 0x07   # unsupported compression method
    g4 = copy(empty_gzip); g4[4] = 0x23   # reserved FLG bit set
    g5 = copy(empty_gzip); g5[11] = 0x92  # header crc error
    g6 = copy(empty_gzip); g6[12] = 0x1f  # header crc error
    for g in [g1, g2, g3, g4, g5, g6]
        @test_throws ErrorException inflate_gzip(g)
        @test_throws ErrorException InflateGzipStream(IOBuffer(g))
    end
    for g in [g5, g6]
        @test isempty(inflate_gzip(g, ignore_checksum = true))
        @test eof(InflateGzipStream(IOBuffer(g), ignore_checksum = true))
    end
end

@testset "Gzip trailer corruption" begin
    g7 = copy(empty_gzip); g7[15] = 0x01   # crc check failed
    g8 = copy(empty_gzip); g8[16] = 0x01   # crc check failed
    g9 = copy(empty_gzip); g9[17] = 0x01   # crc check failed
    g10 = copy(empty_gzip); g10[18] = 0x01   # crc check failed
    g11 = copy(empty_gzip); g11[19] = 0x01   # length check failed
    g12 = copy(empty_gzip); g12[20] = 0x01   # length check failed
    g13 = copy(empty_gzip); g13[21] = 0x01   # length check failed
    g14 = copy(empty_gzip); g14[22] = 0x01   # length check failed
    for g in [g7, g8, g9, g10, g11, g12, g13, g14]
        @test_throws ErrorException inflate_gzip(g)
        @test_throws ErrorException read(InflateGzipStream(IOBuffer(g)))
    end
    for g in [g7, g8, g9, g10]
        @test isempty(inflate_gzip(g, ignore_checksum = true))
        @test eof(InflateGzipStream(IOBuffer(g), ignore_checksum = true))
    end
end

@testset "Reading past end of file" begin
    s1 = InflateStream(IOBuffer(empty_deflate))
    s2 = InflateZlibStream(IOBuffer(empty_zlib))
    s3 = InflateGzipStream(IOBuffer(empty_gzip))
    for s in [s1, s2, s3]
        @test eof(s)
        @test_throws EOFError read(s, UInt8)
    end
end
