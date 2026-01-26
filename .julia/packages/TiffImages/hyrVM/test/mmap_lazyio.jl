@testset "Extant File lazyio" begin
    filepath = get_example("julia.tif")
    img = TiffImages.load(filepath, lazyio=true)
    @test size(img) == (300, 500, 1)
    @test all(img[3, 1:50] .== RGB{N0f8}(1, 1, 1))
    c = img[70,65]   # in the blue dot
    @test blue(c) > max(red(c), green(c))
    # force close the stream behind the file to see if it's properly reopened
    close(img.data.file.io)
    @test all(img[3, 1:50] .== RGB{N0f8}(1, 1, 1))

    # TODO: inefficient convert method
    img_cvt = convert(Array{eltype(img), ndims(img)}, img)
    @test img_cvt == img.data
    @test img_cvt !== img.data
end

@testset "Extant File mmap" begin
    filepath = get_example("flagler.tif")
    let img = TiffImages.load(filepath, mmap=true)
        @test eltype(img) === TiffImages.WidePixel{RGB{N0f8}, Tuple{N0f8}}
        @test size(img) === (200, 541)
        p = img[34, 105]
        @test green(p.color) > blue(p.color) > red(p.color)
        @test p.extra == (1.0,)
        p = img[61, 218]
        @test red(p.color) > blue(p.color) > green(p.color)
        @test_throws BoundsError img[201, 541]
        @test_throws BoundsError img[200, 542]
        @test img[200,541] isa TiffImages.WidePixel{RGB{N0f8}, Tuple{N0f8}}
        imgeager = TiffImages.load(filepath)
        @test img == imgeager

        # Test that attempting to write the mmapped version throws an error,
        # unless we open with write permissions
        imgeager[61, 218] = zero(p)
        @test imgeager[61, 218] === zero(p)
        @test_throws ReadOnlyMemoryError img[61, 218] = zero(p)

        img = TiffImages.load(filepath, mmap=true, mode="r+")
        img[61, 218] = zero(p)
        @test img[61, 218] === zero(p)
        # Put the file back
        img[61, 218] = p
    end

    # 3d
    img0 = rand(Gray{N0f16}, 1000, 1000, 40)
    filepath = tempname() * ".tif"
    TiffImages.save(filepath, img0)
    img1 = TiffImages.load(filepath);
    let img2 = TiffImages.load(filepath; mmap=true), img3 = TiffImages.load(filepath; lazyio=true)
        @test size(img1) == size(img2) == size(img3) == size(img0)
        @test img1[1,2,3] == img2[1,2,3] == img3[1,2,3] == img0[1,2,3]
        c = img0[1,2,3]
        img1[1,2,3] = complement(c)
        @test img1[1,2,3] == complement(c)
        @test_throws ReadOnlyMemoryError img2[1,2,3] = complement(c)
        @test_throws ErrorException img3[1,2,3] = complement(c)
    end

    # with N0f8 (a special case for sizing slice buffers)
    img0 = Gray{N0f8}[0.2 0.4;
                      0   1]
    filepath2 = tempname() * ".tif"
    TiffImages.save(filepath2, img0)
    let img = TiffImages.load(filepath2; mmap=true)
        @test img == img0
    end
end

@testset "De novo construction" begin
    rm("test.tif", force = true)
    img = empty(LazyBufferedTIFF, Gray{N0f8}, "test.tif")

    # a newly initialized file should have every dimension equal to zero and
    # error if accessed
    @test size(img) == (0, 0, 0)
    @test_throws ErrorException img[1, 1, 1]

    push!(img, rand(Gray{N0f8}, 100, 100))
    @test size(img) == (100, 100, 1)
    @test TiffImages.offset(img) == UInt32

    @test_throws AssertionError push!(img, rand(Gray{N0f8}, 99, 99)) # wrong size

    @testset "BigTIFF" begin
        rm("test.btif", force = true)
        img = empty(LazyBufferedTIFF, Gray{N0f8}, "test.btif"; bigtiff = true)

        push!(img, rand(Gray{N0f8}, 100, 100))
        @test size(img) == (100, 100, 1)
        @test TiffImages.offset(img) == UInt64
    end

    img.readonly = true
    @test_throws ErrorException push!(img, rand(Gray{N0f8}, 100, 100))
end
