@testset "encoder" begin
    @testset "Real image encode test" begin
        @info "Check if all real images are visually good"
        img_gray = testimage("mandril_gray")
        img_color = testimage("mandril_color")
        img3d = testimage("mri-stack")
        test_sixel_display() do
            sixel_encode(img_gray)
            sixel_encode(img_color)
            sixel_encode(img3d)
            # higher dimensional array is stacked as if it is 3d
            sixel_encode(reshape(img3d, size(img3d)[1:2]..., 3, :))
        end
    end

    @testset "1d vector" begin
        for img in (
            repeat(Gray.(0:0.1:0.9), inner=10),
            repeat(distinguishable_colors(10), inner=(10, ))
        )
            @info "1d vector: $(eltype(img))"
            @test size(img) == (100, )

            enc = Sixel.default_encoder(img)
            @test enc isa Sixel.LibSixel.LibSixelEncoder
            test_sixel_display() do
                sixel_encode(img, enc)
                println()
            end

            io = IOBuffer()
            sixel_encode(io, img, enc)
            bufferdata = String(take!(io))
            h, w = size(img, 1), 1
            @test startswith(bufferdata, "\ePq\"1;1;$w;$h") && endswith(bufferdata, "\e\\")
            test_sixel_display() do
                println(bufferdata)
            end

            tmp_file = tempname()
            open(tmp_file, "w") do io
                sixel_encode(io, img, enc)
            end
            filedata = read(tmp_file, String)
            @test startswith(filedata, "\ePq\"1;1;$w;$h") && endswith(filedata, "\e\\")
            test_sixel_display() do
                println(filedata)
            end

            @info "transpose check: the first two should look the same, while the third is transposed."
            test_sixel_display() do
                sixel_encode(img, enc)
                sixel_encode(img, enc; transpose=false)
                sixel_encode(img, enc; transpose=true)
            end
        end
    end

    @testset "2d matrix" begin
        for img in (
            repeat(Gray.(0:0.1:0.9), inner=(10, 50)),
            repeat(distinguishable_colors(10), inner=(10, 50))
        )
            @info "2d matrix: $(eltype(img))"
            h, w = size(img)
            @test size(img) == (100, 50)

            enc = Sixel.default_encoder(img)
            @test enc isa Sixel.LibSixel.LibSixelEncoder
            test_sixel_display() do
                sixel_encode(img, enc)
                println()
            end
            io = IOBuffer()
            sixel_encode(io, img, enc)
            bufferdata = String(take!(io))
            @test startswith(bufferdata, "\ePq\"1;1;$w;$h") && endswith(bufferdata, "\e\\")
            test_sixel_display() do
                println(bufferdata)
            end

            tmp_file = tempname()
            open(tmp_file, "w") do io
                sixel_encode(io, img, enc)
            end
            filedata = read(tmp_file, String)
            @test startswith(filedata, "\ePq\"1;1;$w;$h") && endswith(filedata, "\e\\")
            test_sixel_display() do
                println(filedata)
            end

            @info "transpose check: the first two should look the same, while the third is transposed."
            test_sixel_display() do
                sixel_encode(img, enc)
                sixel_encode(img, enc; transpose=false)
                sixel_encode(img, enc; transpose=true)
            end
        end
    end

    @testset "3d array" begin
        for img in (
            repeat(Gray.(0:0.1:0.9), inner=(10, 50, 3)),
            repeat(distinguishable_colors(5), inner=(20, 50, 3))
        )
            @info "3d array: $(eltype(img))"
            @test size(img) == (100, 50, 3)

            enc = Sixel.default_encoder(img)
            @test enc isa Sixel.LibSixel.LibSixelEncoder
            test_sixel_display() do
                sixel_encode(img, enc)
                println()
            end
            io = IOBuffer()
            sixel_encode(io, img, enc)
            bufferdata = String(take!(io))
            h, w, c = size(img)
            w = w * c
            @test startswith(bufferdata, "\ePq\"1;1;$w;$h") && endswith(bufferdata, "\e\\")
            test_sixel_display() do
                println(bufferdata)
            end

            tmp_file = tempname()
            open(tmp_file, "w") do io
                sixel_encode(io, img, enc)
            end
            filedata = read(tmp_file, String)
            h, w, c = size(img)
            w = w * c
            @test startswith(filedata, "\ePq\"1;1;$w;$h") && endswith(filedata, "\e\\")
            test_sixel_display() do
                println(filedata)
            end

            @info "transpose check: the first two should look the same, while the third is transposed."
            test_sixel_display() do
                sixel_encode(img, enc)
                sixel_encode(img, enc; transpose=false)
                sixel_encode(img, enc; transpose=true)
            end

            io = IOBuffer()
            sixel_encode(io, img, enc; transpose=true)
            bufferdata = String(take!(io))
            h, w, c = size(img)
            h = h * c
            h, w = w, h # transpose
            @test startswith(bufferdata, "\ePq\"1;1;$w;$h") && endswith(bufferdata, "\e\\")
        end
    end

    @testset "various eltype" begin
        enc = Sixel.LibSixel.LibSixelEncoder()
        img = repeat(distinguishable_colors(5), inner=(20, 50))

        # test different color spaces
        @info "The images below should visually look the same"
        for T in (
            RGB{N0f8}, RGB{Float32}, RGB24,
            HSV, Lab, ARGB, RGBA, BGR, BGRA
        )
            test_sixel_display() do
                sixel_encode(T.(img), enc)
            end
        end

        # test different storage format
        @info "The images below should visually look the same"
        img = gray.(Gray.(img))
        for T in (Gray, N0f8, Float32, Float64, x->round(Int, 255x))
            test_sixel_display() do
                sixel_encode(T.(img), enc)
            end
        end
    end

    @testset "untypical array types" begin
        # lazy array that does not occupy full memory
        img = Diagonal(repeat(distinguishable_colors(5), inner=20))
        w, h = size(img)
        io = IOBuffer()
        sixel_encode(io, img)
        bufferdata = String(take!(io))
        @test startswith(bufferdata, "\ePq\"1;1;$w;$h") && endswith(bufferdata, "\e\\")
        @info "Test Diagonal array"
        test_sixel_display() do
            println(bufferdata)
        end

        # array that has un-contiguous memory layout
        ori_img = repeat(distinguishable_colors(5), inner=(20, 50))
        img = PermutedDimsArray(permutedims(ori_img, (2, 1)), (2, 1))
        @test size(img) == (100, 50) == size(ori_img)
        @test stride(img, 1) == 50
        @test stride(ori_img, 1) == 1

        io = IOBuffer()
        sixel_encode(io, ori_img)
        bufferdata = String(take!(io))
        h, w = size(ori_img)
        @test startswith(bufferdata, "\ePq\"1;1;$w;$h") && endswith(bufferdata, "\e\\")
        io = IOBuffer()
        sixel_encode(io, img)
        @test bufferdata == String(take!(io))

        io = IOBuffer()
        sixel_encode(io, ori_img, transpose=false)
        bufferdata = String(take!(io))
        h, w = size(ori_img)
        @test startswith(bufferdata, "\ePq\"1;1;$w;$h") && endswith(bufferdata, "\e\\")
        io = IOBuffer()
        sixel_encode(io, img, transpose=false)
        @test bufferdata == String(take!(io))
    end

    @testset "OffsetArray" begin
        enc = Sixel.LibSixelEncoder()
        for img in [
            repeat(distinguishable_colors(10), inner=(10, )),
            repeat(distinguishable_colors(10), inner=(10, 50)),
            repeat(distinguishable_colors(5), inner=(20, 50, 3))
        ]
            io = IOBuffer()
            sixel_encode(io, img, enc)
            ref = String(take!(io))

            imgo = OffsetArray(img, OffsetArrays.Origin(0))
            io = IOBuffer()
            sixel_encode(io, imgo, enc)
            actual = String(take!(io))
            @test ref == actual
        end
    end
end

@testset "decoder" begin
    @testset "API" begin
        dec = Sixel.default_decoder()
        @test dec isa Sixel.LibSixelDecoder

        tmp_file = tempname()        
        img = repeat(Gray.(0:0.1:0.9), inner=(10, 50))
        sz = size(img)
        open(tmp_file, "w") do io
            sixel_encode(io, img)
        end

        # IO
        img_readback = open(tmp_file, "r") do io
            sixel_decode(io)
        end
        @test sz == size(img_readback)
        @test img_readback isa IndirectArray
        @test eltype(img_readback) == RGB{N0f8}

        img_readback = open(tmp_file, "r") do io
            sixel_decode(io, dec)
        end
        @test sz == size(img_readback)
        @test img_readback isa IndirectArray
        @test eltype(img_readback) == RGB{N0f8}

        img_readback = open(tmp_file, "r") do io
            sixel_decode(Gray{N0f8}, io)
        end
        @test sz == size(img_readback)
        @test img_readback isa IndirectArray
        @test eltype(img_readback) == Gray{N0f8}

        # String and Vector{UInt8}
        bytes_data = open(tmp_file, "r") do io
            read(io)
        end
        string_data = read(tmp_file, String)
        for src in (bytes_data, string_data)
            img_readback = sixel_decode(src, dec)
            @test sz == size(img_readback)
            @test img_readback isa IndirectArray
            @test eltype(img_readback) == RGB{N0f8}

            img_readback = sixel_decode(Gray{N0f8}, src, dec)
            @test sz == size(img_readback)
            @test img_readback isa IndirectArray
            @test eltype(img_readback) == Gray{N0f8}
        end
    end

    @testset "Quality test" begin
        enc = Sixel.LibSixelEncoder()
        dec = Sixel.LibSixelDecoder()
        tmp_file = tempname()

        for img in (
            repeat(Gray.(0:0.1:0.9), inner=10),
            repeat(distinguishable_colors(10), inner=(10, )),
            repeat(Gray.(0:0.1:0.9), inner=(10, 50)),
            repeat(distinguishable_colors(10), inner=(10, 50)),
            repeat(Gray.(0:0.1:0.9), inner=(10, 50, 3)),
            repeat(distinguishable_colors(5), inner=(20, 50, 3))
        )
            open(tmp_file, "w") do io
                sixel_encode(io, img, enc)
            end
            
            img_readback = open(tmp_file, "r") do io
                sixel_decode(eltype(img), io, dec)
            end

            img_readback = reshape(img_readback, size(img))
            # 30 is actually pretty good given that sixel encode always do quantization
            @test assess_psnr(img, img_readback) > 30
        end
    end
end
