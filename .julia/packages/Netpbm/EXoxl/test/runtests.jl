using Netpbm, IndirectArrays, ImageCore, FileIO, OffsetArrays, ImageMetadata
using Test

@testset "IO" begin
    workdir = joinpath(tempdir(), "Images")
    isdir(workdir) && rm(workdir, recursive=true)
    mkdir(workdir)

    function compare_metadata(fn, meta)
        info = readlines(filename(fn))
        expect = sort!(["# $k: $v" for (k, v) in properties(meta)])
        actual = sort(view(info, 2:(1+length(properties(meta)))))
        @test expect == actual
    end

    @testset "Bicolor pbm" begin
        # 20 columns = 2.5 bytes
        af = rand(0:1, 3, 20)
        for fmt in (format"PBMBinary", format"PBMText")
            for T in (Bool, Int)
                # test save/load without metadata
                ac = convert(Array{T}, af)
                fn = File{fmt}(joinpath(workdir, "20by3.pbm"))
                Netpbm.save(fn, ac)
                b = Netpbm.load(fn)
                @test b == ac

                # test save with metadata, test metadata content, and correct loading
                meta = ImageMeta(ac, author = "anonymous", time = "now", comment = "something interesting")
                Netpbm.save(fn, meta)
                compare_metadata(fn, meta)
                b = Netpbm.load(fn)
                @test b == ac
            end
        end
    end

    @testset "Gray pgm" begin
        af = rand(2, 3)
        for fmt in (format"PGMBinary", format"PGMText")
            for T in (N0f8, N4f12, N0f16,
                      Gray{N0f8}, Gray{N4f12}, Gray{N0f16})
                # test save/load without metadata
                ac = convert(Array{T}, af)
                fn = File{fmt}(joinpath(workdir, "3by2.pgm"))
                Netpbm.save(fn, ac)
                b = Netpbm.load(fn)
                @test b == ac

                # test save with metadata, test metadata content, and correct loading
                meta = ImageMeta(ac, author = "anonymous", time = "now", comment = "something interesting")
                Netpbm.save(fn, meta)
                compare_metadata(fn, meta)
                b = Netpbm.load(fn)
                @test b == ac
            end
        end
        a8 = convert(Array{N0f8}, af)
        for fmt in (format"PGMBinary", format"PGMText")
            for T in (Float32, Float64, Gray{Float32}, Gray{Float64})
                # test save/load without metadata
                ac = convert(Array{T}, af)
                fn = File{fmt}(joinpath(workdir, "3by2.pgm"))
                Netpbm.save(fn, ac)
                b = Netpbm.load(fn)
                @test b == a8

                # test save with metadata, test metadata content, and correct loading
                meta = ImageMeta(ac, author = "anonymous", time = "now", comment = "something interesting")
                Netpbm.save(fn, meta)
                compare_metadata(fn, meta)
                b = Netpbm.load(fn)
                @test b == a8
            end
        end

        # Issue #11
        @test_nowarn load("test_header_space.pgm")
    end

    @testset "Color ppm" begin
        af = rand(RGB{Float64}, 2, 3)
        for fmt in (format"PPMBinary", format"PPMText")
            for T in (RGB{N0f8}, RGB{N4f12}, RGB{N0f16})
                # test save/load without metadata
                ac = convert(Array{T}, af)
                fn = File{fmt}(joinpath(workdir, "3by2.ppm"))
                Netpbm.save(fn, ac)
                b = Netpbm.load(fn)
                @test b == ac

                # test save with metadata, test metadata content, and correct loading
                meta = ImageMeta(ac, author = "anonymous", time = "now", comment = "something interesting")
                Netpbm.save(fn, meta)
                compare_metadata(fn, meta)
                b = Netpbm.load(fn)
                @test b == ac
            end
        end
        a8 = convert(Array{RGB{N0f8}}, af)
        for fmt in (format"PPMBinary", format"PPMText")
            for T in (RGB{Float32}, RGB{Float64}, HSV{Float64})
                # test save/load without metadata
                ac = convert(Array{T}, af)
                fn = File{fmt}(joinpath(workdir, "3by2.ppm"))
                Netpbm.save(fn, ac)
                b = Netpbm.load(fn)
                @test b == a8

                # test save with metadata, test metadata content, and correct loading
                meta = ImageMeta(ac, author = "anonymous", time = "now", comment = "something interesting")
                Netpbm.save(fn, meta)
                compare_metadata(fn, meta)
                b = Netpbm.load(fn)
                @test b == a8
           end
        end
    end

    @testset "Colormap" begin
        datafloat = reshape(range(0.5, stop=1.5, length=6), 2, 3)
        dataint = map(x->round(UInt8, x), 254*(datafloat .- 0.5) .+ 1) # ranges 1 to 255
        # build our colormap
        b = RGB(0,0,1)
        w = RGB(1,1,1)
        r = RGB(1,0,0)
        cmaprgb = Array{RGB{N0f8}}(undef, 255)
        f = range(0, stop=1, length=128)
        cmaprgb[1:128] = [(1-x)*b + x*w for x in f]
        cmaprgb[129:end] = [(1-x)*w + x*r for x in f[2:end]]
        img = IndirectArray(dataint, cmaprgb)
        fn = File{format"PPMBinary"}(joinpath(workdir,"cmap.ppm"))
        Netpbm.save(fn, img)
        imgr = Netpbm.load(fn)
        @test imgr == img
        cmaprgb = Array{RGB}(undef, 255) # poorly-typed cmap, Images issue #336
        @test !isconcretetype(eltype(cmaprgb))
        cmaprgb[1:128] = RGB{N0f16}[(1-x)*b + x*w for x in f]
        cmaprgb[129:end] = RGB{N4f12}[(1-x)*w + x*r for x in f[2:end]]
        img = IndirectArray(dataint, cmaprgb)
        @test_throws ErrorException Netpbm.save(fn, img) # widens to unsupported type
        cmaprgb[129:end] = RGB{N0f8}[(1-x)*w + x*r for x in f[2:end]]
        img = IndirectArray(dataint, cmaprgb)
        Netpbm.save(fn, img)
        imgr = Netpbm.load(fn)
        @test imgr == img
    end

    # Images issue #256
    @testset "Clamping" begin
        A = rand(2,3)
        A[1,1] = -0.4
        fn = File{format"PGMBinary"}(joinpath(workdir, "2by3.pgm"))
        @test_throws InexactError Netpbm.save(fn, A)
        Netpbm.save(fn, A, mapf=clamp01nan)
        B = Netpbm.load(fn)
        A[1,1] = 0
        @test B == Gray{N0f8}.(A)
    end

    @testset "OffsetArrays" begin
        ac = OffsetArray(rand(Bool, 3, 20), -1:1, 0:19)
        fn = File{format"PBMText"}(joinpath(workdir, "20by3.pbm"))
        Netpbm.save(fn, ac)
        b = Netpbm.load(fn)
        @test b[1:end] == ac[1:end]

        ac = OffsetArray(rand(Gray{N0f8}, 2, 3), 0:1, -1:1)
        fn = File{format"PGMText"}(joinpath(workdir, "3by2.pgm"))
        Netpbm.save(fn, ac)
        b = Netpbm.load(fn)
        @test b[1:end] == ac[1:end]

        ac = OffsetArray(rand(RGB{N0f8}, 2, 3), 0:1, -1:1)
        fn = File{format"PPMText"}(joinpath(workdir, "3by2.ppm"))
        Netpbm.save(fn, ac)
        b = Netpbm.load(fn)
        @test b[1:end] == ac[1:end]
    end
end

nothing
