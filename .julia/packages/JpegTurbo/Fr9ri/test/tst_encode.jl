@testset "jpeg_encode" begin

img_rgb = testimage("lighthouse")

@testset "basic" begin
    for CT in [Gray, RGB, #=YCbCr,=# #=RGBX,=# BGR, #=XRGB,=# RGBA, BGRA, ABGR, ARGB]
        img = CT.(img_rgb)
        data = jpeg_decode(jpeg_encode(img))
        @test eltype(data) <: Union{Gray, RGB}
        @test size(data) == size(img)
        @test data ≈ jpeg_decode(jpeg_encode(float32.(img)))

        # ensure default keyword values are not changed by accident
        @test data == jpeg_decode(jpeg_encode(img, transpose=false))
        @test jpeg_decode(jpeg_encode(img, transpose=true)) ==
            jpeg_decode(jpeg_encode(img', transpose=false))
    end

    # numerical array is treated as Gray image
    img = Gray.(img_rgb)
    @test jpeg_encode(Float32.(img)) == jpeg_encode(img)

    # out-of-range values are mapped into [0, 1]
    img_or = 1.5 .* img
    img_or[1] = Gray(NaN)
    @test jpeg_encode(img_or) == jpeg_encode(Float64.(img_or)) == jpeg_encode(clamp01nan!(img_or))

    @test_throws ArgumentError("empty image is not allowed") jpeg_encode(Array{Float64,2}(undef, 0, 0))
end

# keyword checks
@testset "quality" begin
    img = testimage("cameraman")
    psnr_refs = [
        1   => 24.63,
        10  => 31.34,
        50  => 38.87,
        100 => 56.74,
    ]
    for (q, r) in psnr_refs
        v = assess_psnr(img, jpeg_decode(jpeg_encode(img, quality=q)))
        @test v >= r
    end
end

@testset "progressive" begin
    tmpfile = joinpath(tmpdir, "tmp.jpg")
    img = testimage("cameraman")
    progressive_modes = [nothing, false, true]
    for progressive_mode ∈ progressive_modes
        bytes = if isnothing(progressive_mode)
            JpegTurbo.jpeg_encode(tmpfile, img)
            jpeg_encode(img)
        else
            JpegTurbo.jpeg_encode(tmpfile, img; progressive_mode=progressive_mode)
            jpeg_encode(img; progressive_mode=progressive_mode)
        end
        if progressive_mode ∈ [nothing, false]
            @test is_progressive_jpeg(tmpfile) == false
            @test is_progressive_jpeg(bytes) == false
        else
            @test is_progressive_jpeg(tmpfile) == true
            @test is_progressive_jpeg(bytes) == true
        end
    end
end

end # @testset "jpeg_encode"
