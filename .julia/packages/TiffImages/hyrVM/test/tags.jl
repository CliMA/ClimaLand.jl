@testset "Unspecified type" begin
    tf = TiffImages.TiffFile{UInt32}()

    write(tf, UInt16(TiffImages.IMAGEWIDTH))
    write(tf, 0x0007)
    write(tf, UInt32(2))
    write(tf, UInt16[8, 1])

    seekstart(tf.io)
    tag = read(tf, TiffImages.Tag)
    @test tag == TiffImages.Tag(TiffImages.IMAGEWIDTH, Any[8, 0, 1, 0])
    @test sizeof(tag) == 4
end

@testset "Data array only part of data field" begin
    tf = TiffImages.TiffFile{UInt64}()

    write(tf, UInt16(TiffImages.BITSPERSAMPLE))
    write(tf, 0x0003)
    write(tf, UInt64(2))
    write(tf, UInt16[8, 8])
    write(tf, UInt32(0))

    seekstart(tf.io)
    tag = read(tf, TiffImages.Tag)
    @test length(tag.data) == 2 # not 4
    @test sizeof(tag) == 4 # number of bytes is still 2*2
end

@testset "Rational, full space" begin
    tf = TiffImages.TiffFile{UInt64}()
    write(tf, UInt16(TiffImages.XRESOLUTION))
    write(tf, 0x0005)
    write(tf, UInt64(1))
    ratio = Rational{UInt32}(1, 20)
    write(tf, ratio)

    seekstart(tf.io)
    tag = read(tf, TiffImages.Tag)
    @test tag.data == ratio
    @test sizeof(tag) == 8 # two UInt64s 
end

@testset "Additional sizeof tests" begin
    tag = TiffImages.Tag(TiffImages.IMAGEDESCRIPTION, "test")
    @test sizeof(tag) == 5 # 4 bytes + NUL terminator
end