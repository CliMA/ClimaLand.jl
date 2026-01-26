@testset "tags" begin
    @testset "exact fit tags" begin
        tf = TiffImages.TiffFile{UInt32}()
        t1 = TiffImages.Tag(UInt16(TiffImages.IMAGEWIDTH), 0x00000200)
        @test write(tf, t1) # should fit so true

        seekstart(tf.io)
        t2 = read(tf, TiffImages.Tag)
        @test t1 == t2
    end

    @testset "small tags" begin
        tf = TiffImages.TiffFile{UInt32}()
        t1 = TiffImages.Tag(UInt16(TiffImages.COMPRESSION), 0x0001)
        @test write(tf, t1) # should fit but needs padding
        @test position(tf.io) == 12 # should be full length

        seekstart(tf.io)
        t2 = read(tf, TiffImages.Tag)
        @test t2.data == 0x0001
    end

    @testset "large tags" begin
        tf = TiffImages.TiffFile{UInt32}()
        offsets = UInt32[8, 129848, 259688, 389528]
        t1 = TiffImages.Tag(UInt16(TiffImages.STRIPOFFSETS), offsets)

        @test !write(tf, t1) # should fail to write into tf since it's too large
        @test position(tf.io) == 12

        data = Array{UInt8}(undef, 12)
        seekstart(tf.io)
        read!(tf.io, data)

        @test all(data .== 0x00)

        seekstart(tf.io)
        @test write(tf, t1, 12)

        seekstart(tf.io)
        t2 = read(tf, TiffImages.Tag)

        # make sure the data field contains our single offset from above
        @test typeof(t2) <: TiffImages.Tag{<: TiffImages.RemoteData}
        @test getfield(t2, :data).position == 12
    end

    @testset "String length equal to offset size" begin
        tf = TiffImages.TiffFile{UInt32}()
        t1 = TiffImages.Tag(TiffImages.SOFTWARE, "test")

        # should fail because it's too large to fit
        @test !write(tf, t1)
        @test position(tf.io) == 12

        t3 = TiffImages.Tag(TiffImages.SOFTWARE, "tes\0")
        seekstart(tf.io)
        @test write(tf, t3)
        @test position(tf.io) == 12

        seekstart(tf.io)
        t4 = read(tf, TiffImages.Tag)

        # Since this is a NUL-terminated string of length 4 it should fit
        @test getfield(t4, :data) == "tes\0"
        @test t4.data == "tes"
    end

    @testset "UTF-8 strings" begin
        ifd = TiffImages.IFD(UInt32)
        ifd[TiffImages.IMAGEDESCRIPTION] = "αβγ" # 2 bytes each, 6 total
        ifd[TiffImages.SOFTWARE] = "∫" # 3 byte character

        tf = TiffImages.TiffFile{UInt32}()
        write(tf, ifd)

        ifd2 = first(tf)
        TiffImages.load!(tf, ifd2)

        @test ifd[TiffImages.SOFTWARE] == ifd2[TiffImages.SOFTWARE]
        @test ifd[TiffImages.IMAGEDESCRIPTION] == ifd2[TiffImages.IMAGEDESCRIPTION] # test remote data utf-8
    end
end

@testset "ifds" begin
    tf = TiffImages.TiffFile{UInt32}()
    ifd = TiffImages.IFD(UInt32)
    ifd[TiffImages.IMAGEDESCRIPTION] = "Testing IFD read/write"
    ifd[TiffImages.IMAGEWIDTH] = UInt32(512)
    ifd[TiffImages.IMAGELENGTH] = UInt32(512)
    ifd[TiffImages.XRESOLUTION] = Rational{UInt32}(72, 1)
    ifd[TiffImages.COMPRESSION] = TiffImages.COMPRESSION_NONE
    ifd[TiffImages.SAMPLESPERPIXEL] = 1
    ifd[TiffImages.PHOTOMETRIC] = TiffImages.PHOTOMETRIC_MINISBLACK
    ifd[TiffImages.STRIPBYTECOUNTS] = UInt32[512*512]
    ifd[TiffImages.BITSPERSAMPLE] = 8
    ifd[TiffImages.ICCPROFILE] = Any[0x00, 0x00, 0x0b, 0xe8, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x6d]

    write(tf, ifd)

    seekstart(tf.io)
    read_ifd, next_ifd = read(tf, TiffImages.IFD)
    
    seekend(tf.io)
    @test sizeof(ifd) == position(tf.io) # predicted size should match actual size
    @test sizeof(read_ifd) == position(tf.io) # unloaded IFDs should also compute properly

    TiffImages.load!(tf, read_ifd)

    @test all(sort(collect(ifd), by = x -> x[1]) .== sort(collect(read_ifd), by = x -> x[1]))
end

@testset "Simple 2D image" begin
    filepath = get_example("house.tif")
    img = TiffImages.load(filepath)

    img2 = TiffImages.DenseTaggedImage(Gray{Float64}.(img.data));
    path, io = mktemp()
    write(io, img2)
    img3 = TiffImages.load(path)
    @test eltype(img3) == Gray{Float64}
end

@testset "3D image" begin
    filepath = get_example("mri.tif")
    img = TiffImages.load(filepath)

    img2 = TiffImages.DenseTaggedImage(RGB{Float64}.(img.data));
    path, io = mktemp()
    write(io, img2)
    img3 = TiffImages.load(path)
    @test eltype(img3) == RGB{Float64}
end

@testset "BigTIFF saving" begin
    filepath = get_example("house.tif")
    img = TiffImages.load(filepath)

    ifd = TiffImages._constructifd(img.data, UInt64);
    img2 = TiffImages.DenseTaggedImage(img.data, ifd)
    path, io = mktemp()
    write(io, img2)
    img3 = TiffImages.load(path)
    @test eltype(img3) == GrayA{N0f8}
end

@testset "DenseTaggedImage with OffsetArray (#46)" begin
    for T in (Gray, RGB)
        data = OffsetArray(rand(T, 10, 10), -1, -1);
        filepath = tempname()
        TiffImages.save(filepath, data)

        data_saveload = TiffImages.load(filepath)
        # the offset is stripped when saving
        @test axes(data_saveload) != axes(data)
        @test data_saveload == collect(data)
    end
end

@testset "DenseTaggedImage with AxisArray" begin
    for T in (Gray, RGB)
        data = AxisArray(rand(Gray, 100, 100, 3), :h, :w)
        filepath = tempname()
        TiffImages.save(filepath, data)

        # the axis information is stripped when saving
        data_saveload = TiffImages.load(filepath)
        @test axes(data_saveload) == axes(data)
        @test data_saveload == collect(data)
    end
end

@testset "save function" begin
    data = rand(RGB{N0f8}, 128, 128)
    filepath = tempname()

    TiffImages.save(filepath, data)
    @test all(data .== TiffImages.load(filepath))
end


@testset "Software tag, issue #60" begin
    filepath = get_example("house.tif")
    img = TiffImages.load(filepath)
    ifd = TiffImages.IFD(UInt32)
    ifd[TiffImages.SOFTWARE] = "test"

    img2 = TiffImages.DenseTaggedImage(Gray{Float64}.(img.data), ifd);
    path, io = mktemp()
    write(io, img2)
    img3 = TiffImages.load(path)
    @test occursin("test;", ifds(img3)[TiffImages.SOFTWARE].data)
end

@testset "DenseTaggedImage with unusual bit depth (#181)" begin
    filepath = tempname()
    for C in (Gray, RGB)
        for T in (N0f8, N3f5, N4f12, N13f19, N12f52)
            for size in (1,10,100,567)
                data = C{T}.(rand(size, size))
                TiffImages.save(filepath, data)

                data_saveload = TiffImages.load(filepath)

                @test data_saveload.data == data
            end
        end
    end
end
