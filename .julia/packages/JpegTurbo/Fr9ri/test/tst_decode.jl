@testset "jpeg_decode" begin
    img_rgb = testimage("lighthouse")
    img_rgb_bytes = jpeg_encode(img_rgb)

    # ensure default keyword values are not changed by accident
    @test jpeg_decode(img_rgb_bytes) ≈
        jpeg_decode(RGB, img_rgb_bytes; transpose=false, scale_ratio=1) ≈
        jpeg_decode(img_rgb_bytes; transpose=false, scale_ratio=1)

    @testset "filename and IOStream" begin
        tmpfile = joinpath(tmpdir, "tmp.jpg")
        jpeg_encode(tmpfile, img_rgb)
        @test read(tmpfile) == img_rgb_bytes

        # IOStream
        img = open(tmpfile, "r") do io
            jpeg_decode(io)
        end
        @test img == jpeg_decode(img_rgb_bytes)

        img = open(tmpfile, "r") do io
            jpeg_decode(Gray, io; scale_ratio=0.5)
        end
        @test img == jpeg_decode(Gray, img_rgb_bytes; scale_ratio=0.5)

        # filename
        @test jpeg_decode(tmpfile) == jpeg_decode(img_rgb_bytes)
        @test jpeg_decode(Gray, tmpfile; scale_ratio=0.5) == jpeg_decode(Gray, img_rgb_bytes; scale_ratio=0.5)
    end

    @testset "colorspace" begin
        native_color_spaces = [Gray, RGB, BGR, RGBA, BGRA, ABGR, ARGB]
        ext_color_spaces = [YCbCr, RGBX, XRGB, Lab, YIQ] # supported by Colors.jl
        for CT in [native_color_spaces..., ext_color_spaces...]
            data = jpeg_decode(CT, img_rgb_bytes)
            @test eltype(data) <: CT
            if CT == Gray
                @test assess_psnr(data, Gray.(img_rgb)) > 34.92
            else
                @test assess_psnr(RGB.(data), img_rgb) > 33.87
            end
        end
    end

    @testset "scale_ratio" begin
        data = jpeg_decode(img_rgb_bytes; scale_ratio=0.25)
        @test size(data) == (128, 192) == 0.25 .* size(img_rgb)

        # `jpeg_decode` will map input `scale_ratio` to allowed values.
        data = jpeg_decode(img_rgb_bytes; scale_ratio=0.3)
        @test size(data) == (128, 192) != 0.3 .* size(img_rgb)
    end

    @testset "preferred_size" begin
        # lighthouse: (512, 768)
        actual_size = size(img_rgb)
        data = jpeg_decode(img_rgb_bytes; scale_ratio=0.25)
        data_new = jpeg_decode(img_rgb_bytes; preferred_size=size(data))
        @test size(data_new) == size(data)

        data_new = jpeg_decode(img_rgb_bytes; preferred_size=(>=, size(data)))
        @test size(data_new) == size(data)

        data_new = jpeg_decode(img_rgb_bytes; preferred_size=(>, size(data)))
        @test size(data_new) == 3/8 .* actual_size

        data_new = jpeg_decode(img_rgb_bytes; preferred_size=(<=, size(data)))
        @test size(data_new) == size(data)

        data_new = jpeg_decode(img_rgb_bytes; preferred_size=(<, size(data)))
        @test size(data_new) == 1/8 .* actual_size

        data_new = @suppress_err jpeg_decode(img_rgb_bytes; preferred_size=3 .* actual_size)
        @test size(data_new) == 2 .* actual_size
        msg = @capture_err jpeg_decode(img_rgb_bytes; preferred_size=3 .* actual_size)
        @test occursin("Warning: Failed to infer appropriate scale ratio, use `scale_ratio=2` instead.", msg)

        data_new = @suppress_err jpeg_decode(img_rgb_bytes; preferred_size=(<=, 1/16 .* actual_size))
        @test size(data_new) == 1/8 .* actual_size
        msg = @capture_err jpeg_decode(img_rgb_bytes; preferred_size=(<=, 1/16 .* actual_size))
        @test occursin("Warning: Failed to infer appropriate scale ratio, use `scale_ratio=1/8` instead.", msg)

        @test_throws ArgumentError jpeg_decode(img_rgb_bytes; preferred_size=size(data), scale_ratio=1)
    end

    @testset "transpose" begin
        data = jpeg_decode(jpeg_encode(img_rgb; transpose=true); transpose=true)
        @test assess_psnr(data, img_rgb) > 33.95
    end

    @testset "integrity check" begin
        @test_throws ArgumentError jpeg_decode(UInt8[])
        @test_throws ArgumentError jpeg_decode(img_rgb_bytes[1:600])
    end
end
