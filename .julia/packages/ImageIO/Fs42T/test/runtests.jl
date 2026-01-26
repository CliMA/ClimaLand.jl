using Test
using ImageIO
using FileIO: File, DataFormat, Stream, @format_str
using ImageCore: N0f8, RGB, Gray, RGBA, GrayA, n0f8
using IndirectArrays
using ImageQualityIndexes

tmpdir = mktempdir()
Threads.nthreads() <= 1 && @info "Threads.nthreads() = $(Threads.nthreads()), multithread tests will be disabled"
@testset "ImageIO" begin

    @testset "PNGs" begin
        if Threads.nthreads() > 1 ## NOTE: This test must go first, to test that PNGFiles loading behaves correctly
            @testset "Threaded save" begin
                # test that loading of PNGFiles happens sequentially and doesn't segfault
                img = rand(UInt8, 10, 10)
                Threads.@threads for i in 1:Threads.nthreads()
                    f = File{DataFormat{:PNG}}(joinpath(tmpdir, "test_fpath_$i.png"))
                    ImageIO.save(f, img)
                end
            end
        end
        for typ in [UInt8, N0f8, Gray{N0f8}, RGB{N0f8}] # TODO: Fix 0:1, Gray{Float64}, RGB{Float64} in PNGFiles
            @testset "$typ PNG" begin
                img = rand(typ, 10, 10)
                f = File{DataFormat{:PNG}}(joinpath(tmpdir, "test_fpath.png"))
                ImageIO.save(f, img)
                img_saveload = ImageIO.load(f)
                if typ == UInt8
                    @test all(img .== reinterpret(UInt8, img_saveload))
                else
                    @test img == img_saveload
                end
                @test img_saveload isa Array
                img_saveload = ImageIO.load(f; gamma=1.0)
                if typ == UInt8
                    @test all(img .== reinterpret(UInt8, img_saveload))
                else
                    @test img == img_saveload
                end

                open(io->ImageIO.save(Stream{format"PNG"}(io), img, permute_horizontal=false), joinpath(tmpdir, "test_io.png"), "w")
                img_saveload = open(io->ImageIO.load(Stream{format"PNG"}(io)), joinpath(tmpdir, "test_io.png"))
                if typ == UInt8
                    @test all(img .== reinterpret(UInt8, img_saveload))
                else
                    @test img == img_saveload
                end
                @test img_saveload isa Array

                ImageIO.save(f, img; compression_level=1)
                img_saveload = ImageIO.load(f)
                if typ == UInt8
                    @test all(img .== reinterpret(UInt8, img_saveload))
                else
                    @test img == img_saveload
                end
            end
        end

        @testset "indexed image" begin
            f = File{DataFormat{:PNG}}(joinpath(tmpdir, "test_fpath.png"))
            img = IndirectArray(rand(1:5, 4, 4), rand(RGB{N0f8}, 5))
            ImageIO.save(f, img)
            img_saveload = ImageIO.load(f)

            @test img == img_saveload
            @test img isa IndirectArray
        end
    end

    @testset "Portable bitmap" begin
        @testset "Bicolor pbm" begin
            img = rand(0:1, 10, 10)
            for fmt in (format"PBMBinary", format"PBMText")
                f = File{fmt}(joinpath(tmpdir, "test_fpath.pbm"))
                ImageIO.save(f, img)
                img_saveload = ImageIO.load(f)
                @test img == img_saveload
                @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)

                open(io->ImageIO.save(Stream{fmt}(io), img), joinpath(tmpdir, "test_io.pbm"), "w")
                img_saveload = open(io->ImageIO.load(Stream{fmt}(io)), joinpath(tmpdir, "test_io.pbm"))
                @test img == img_saveload
                @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)
            end
        end

        @testset "Gray pgm" begin
            img = rand(N0f8, 10, 10)
            for fmt in (format"PGMBinary", format"PGMText")
                f = File{fmt}(joinpath(tmpdir, "test_fpath.pgm"))
                ImageIO.save(f, img)
                img_saveload = ImageIO.load(f)
                @test img == img_saveload
                @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)

                open(io->ImageIO.save(Stream{fmt}(io), img), joinpath(tmpdir, "test_io.pgm"), "w")
                img_saveload = open(io->ImageIO.load(Stream{fmt}(io)), joinpath(tmpdir, "test_io.pgm"))
                @test img == img_saveload
                @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)
            end
        end

        @testset "Color ppm" begin
            img = rand(RGB{N0f8}, 10, 10)
            for fmt in (format"PPMBinary", format"PPMText")
                f = File{fmt}(joinpath(tmpdir, "test_fpath.ppm"))
                ImageIO.save(f, img)
                img_saveload = ImageIO.load(f)
                @test img == img_saveload
                @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)

                open(io->ImageIO.save(Stream{fmt}(io), img), joinpath(tmpdir, "test_io.ppm"), "w")
                img_saveload = open(io->ImageIO.load(Stream{fmt}(io)), joinpath(tmpdir, "test_io.ppm"))
                @test img == img_saveload
                @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)
            end
        end
    end

    @testset "TIFF" begin
        for typ in [Gray{N0f8}, Gray{Float64}, RGB{N0f8}, RGB{Float64}] # TODO: Add UInt8, N0f8 support in TiffImages
            @testset "$typ TIFF" begin
                img = rand(typ, 10, 10)
                f = File{format"TIFF"}(joinpath(tmpdir, "test_fpath_$(typ).tiff"))
                ImageIO.save(f, img)
                img_saveload = ImageIO.load(f)
                @test img == img_saveload
                @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)
                img_saveload = ImageIO.load(f; mmap=true)
                @test img == reshape(img_saveload, size(img))

                open(io->ImageIO.save(Stream{format"TIFF"}(io), img), joinpath(tmpdir, "test_io_$(typ).tiff"), "w")
                img_saveload = open(io->ImageIO.load(Stream{format"TIFF"}(io)), joinpath(tmpdir, "test_io_$(typ).tiff"))
                @test img == img_saveload
                @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)

                # mmapped images should not canonicalize by default, and can be controlled manually
                img_saveload = ImageIO.load(f; mmap=true)
                @test typeof(img_saveload) != ImageIO.canonical_type(f, img_saveload)
                img_saveload = ImageIO.load(f; canonicalize=false)
                @test typeof(img_saveload) != ImageIO.canonical_type(f, img_saveload)
                img_saveload = open(io->ImageIO.load(Stream{format"TIFF"}(io); mmap=true), joinpath(tmpdir, "test_io_$(typ).tiff"))
                @test typeof(img_saveload) != ImageIO.canonical_type(f, img_saveload)
                img_saveload = open(io->ImageIO.load(Stream{format"TIFF"}(io); canonicalize=false), joinpath(tmpdir, "test_io_$(typ).tiff"))
                @test typeof(img_saveload) != ImageIO.canonical_type(f, img_saveload)
            end
        end
    end

    @testset "EXR" begin
        for typ in [RGBA{Float16}, RGB{Float16}, GrayA{Float16}, Gray{Float16}]
            img = rand(typ, 10, 10)
            f = File{format"EXR"}(joinpath(tmpdir, "test_fpath.exr"))
            ImageIO.save(f, img)
            img_saveload = ImageIO.load(f)
            @test img == img_saveload
            @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)
        end
    end

    @testset "QOI" begin
        for typ in [RGBA{N0f8}, RGB{N0f8}]
            img = rand(typ, 10, 10)
            f = File{format"QOI"}(joinpath(tmpdir, "test_fpath.qoi"))
            ImageIO.save(f, img)
            img_saveload = ImageIO.load(f)
            @test img == img_saveload
            @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)
        end
    end

    @testset "sixel" begin
        for typ in [Gray{N0f8}, Gray{Float64}, RGB{N0f8}, RGB{Float64}]
            @testset "$typ sixel" begin
                img = repeat(typ.(0:0.1:0.9), inner=(10, 50))
                f = File{format"SIXEL"}(joinpath(tmpdir, "test_fpath.sixel"))
                ImageIO.save(f, img)
                img_saveload = ImageIO.load(f)
                @test eltype(img_saveload) == RGB{N0f8} # currently Sixel forces eltype to be RGB{N0f8}
                @test assess_psnr(img, eltype(img).(img_saveload)) > 30 # Sixel encode involves lossy quantization and dither operations
                @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)

                open(io->ImageIO.save(Stream{format"SIXEL"}(io), img), joinpath(tmpdir, "test_io.sixel"), "w")
                img_saveload = open(io->ImageIO.load(Stream{format"SIXEL"}(io)), joinpath(tmpdir, "test_io.sixel"))
                @test eltype(img_saveload) == RGB{N0f8} # currently Sixel forces eltype to be RGB{N0f8}
                @test assess_psnr(img, eltype(img).(img_saveload)) > 30 # Sixel encode involves lossy quantization and dither operations
                @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)
            end
        end
    end

    @testset "JPEG" begin
        for typ in [Gray{N0f8}, Gray{Float64}, RGB{N0f8}, RGB{Float64}]
            @testset "$typ JPEG" begin
                img = repeat(typ.(0:0.1:0.9), inner=(10, 50))
                f = File{format"JPEG"}(joinpath(tmpdir, "test_fpath.jpg"))
                ImageIO.save(f, img)
                img_saveload = ImageIO.load(f)
                @test eltype(img_saveload) == n0f8(typ) # JpegTurbo uses 8bit
                @test assess_psnr(img, eltype(img).(img_saveload)) > 51
                @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)

                open(io->ImageIO.save(Stream{format"JPEG"}(io), img), joinpath(tmpdir, "test_io.jpg"), "w")
                img_saveload = open(io->ImageIO.load(Stream{format"JPEG"}(io)), joinpath(tmpdir, "test_io.jpg"))
                @test eltype(img_saveload) == n0f8(typ) # JpegTurbo uses 8bit
                @test assess_psnr(img, eltype(img).(img_saveload)) > 51
                @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)
            end
        end
    end

    @testset "WebP" begin
        for typ in [RGBA{N0f8}, RGB{N0f8}]
            @testset "$typ WebP" begin
                img = rand(typ, 10, 10)
                f = File{format"WebP"}(joinpath(tmpdir, "test_fpath.webp"))
                ImageIO.save(f, img)
                img_saveload = ImageIO.load(f)
                @test eltype(img_saveload) == n0f8(typ) # WebP uses 8bit
                @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)

                open(io->ImageIO.save(Stream{format"WebP"}(io), img), joinpath(tmpdir, "test_io.webp"), "w")
                img_saveload = open(io->ImageIO.load(Stream{format"WebP"}(io)), joinpath(tmpdir, "test_io.webp"))
                @test eltype(img_saveload) == n0f8(typ) # WebP uses 8bit
                @test typeof(img_saveload) == ImageIO.canonical_type(f, img_saveload)
            end
        end
    end
end
