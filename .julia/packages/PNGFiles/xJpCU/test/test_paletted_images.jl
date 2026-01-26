synth_paletted_imgs = [
    ("RGB_paletted", IndirectArray(rand(1:100, 127, 257), rand(RGB{N0f8}, 100))),
    ("BRG_paletted", IndirectArray(rand(1:100, 127, 257), rand(BGR{N0f8}, 100))),
    ("ARGB_paletted", IndirectArray(rand(1:100, 127, 257), rand(ARGB{N0f8}, 100))),
    ("ABGR_paletted", IndirectArray(rand(1:100, 127, 257), rand(ABGR{N0f8}, 100))),
    ("RGBA_paletted", IndirectArray(rand(1:100, 127, 257), rand(RGBA{N0f8}, 100))),
    ("RGB_paletted_float", IndirectArray(rand(1:100, 127, 257), rand(RGB{Float64}, 100))),
    ("BRG_paletted_float", IndirectArray(rand(1:100, 127, 257), rand(BGR{Float64}, 100))),
    ("ARGB_paletted_float", IndirectArray(rand(1:100, 127, 257), rand(ARGB{Float64}, 100))),
    ("ABGR_paletted_float", IndirectArray(rand(1:100, 127, 257), rand(ABGR{Float64}, 100))),
    ("RGBA_paletted_float", IndirectArray(rand(1:100, 127, 257), rand(RGBA{Float64}, 100))),
]

expected_img(x::Matrix{<:TransparentColor}) = RGBA{N0f8}.(x)
expected_img(x::Matrix{<:AbstractRGB}) = RGB{N0f8}.(x)

@testset "synthetic paletted images" begin
    for expand_paletted in (true, false)
        @testset "expand_paletted = $(expand_paletted)" begin
            for (case, image) in synth_paletted_imgs
                @testset "$(case)" begin
                    expected = expected_img(collect(image))
                    fpath = joinpath(PNG_TEST_PATH, "test_img_$(case).png")
                    @testset "write" begin
                        open(io -> PNGFiles.save(io, image), fpath, "w")
                        @test PNGFiles.save(fpath, image) === nothing
                    end
                    @testset "read" begin
                        global read_in_pngf = PNGFiles.load(fpath, expand_paletted=expand_paletted)
                        @test typeof(read_in_pngf) <: AbstractMatrix
                    end
                    @testset "compare" begin
                        @test all(expected .≈ read_in_pngf)
                    end
                    global read_in_immag = _standardize_grayness(ImageMagick.load(fpath))
                    @testset "$(case): ImageMagick read type equality" begin
                        # The lena image is Grayscale saved as RGB...
                        @test eltype(_standardize_grayness(read_in_pngf)) == eltype(read_in_immag)
                    end
                    @testset "$(case): ImageMagick read values equality" begin
                        imdiff_val = imdiff(read_in_pngf, read_in_immag)
                        onfail(@test imdiff_val < 0.01) do
                            PNGFiles._inspect_png_read(fpath)
                            _add_debugging_entry(fpath, case, imdiff_val)
                        end
                    end
                    path, ext = splitext(fpath)
                    newpath = path * "_new" * ext
                    PNGFiles.save(newpath, read_in_pngf)
                    @testset "$(case): IO is idempotent" begin
                        @test all(read_in_pngf .≈ PNGFiles.load(newpath, expand_paletted=expand_paletted))
                        @test all(read_in_pngf .≈ open(io -> PNGFiles.load(io, expand_paletted=expand_paletted), newpath))
                    end
                end
            end
        end
    end
end

@testset "palleted image with fewer transparency values than color values" begin
    expected_palette = RGBA{N0f8}[
        RGBA{N0f8}(0.0, 0.0, 0.0, 0.0),
        RGBA{N0f8}(0.282, 0.133, 0.169, 0.024),
        RGBA{N0f8}(0.502, 0.431, 0.459, 0.337),
        RGBA{N0f8}(0.49, 0.427, 0.478, 0.51),
        RGBA{N0f8}(0.282, 0.133, 0.169, 0.208),
        RGBA{N0f8}(0.0, 0.0, 0.0, 0.004),
        RGBA{N0f8}(0.714, 0.518, 0.153, 0.722),
        RGBA{N0f8}(0.6, 0.49, 0.322, 0.953),
        RGBA{N0f8}(0.298, 0.18, 0.208, 0.718),
        RGBA{N0f8}(0.196, 0.11, 0.11, 0.408),
        RGBA{N0f8}(0.0, 0.0, 0.0, 0.039),
        RGBA{N0f8}(0.655, 0.765, 0.18, 0.227),
        RGBA{N0f8}(0.592, 0.686, 0.18, 0.212),
        RGBA{N0f8}(0.961, 0.839, 0.153, 0.957),
        RGBA{N0f8}(0.58, 0.412, 0.11, 0.831),
        RGBA{N0f8}(0.0, 0.0, 0.0, 0.345),
        RGBA{N0f8}(0.0, 0.0, 0.0, 0.188),
        RGBA{N0f8}(0.0, 0.0, 0.0, 0.02),
        RGBA{N0f8}(0.667, 0.757, 0.208, 0.749),
        RGBA{N0f8}(0.616, 0.69, 0.196, 0.988),
        RGBA{N0f8}(0.584, 0.365, 0.169, 0.996),
        RGBA{N0f8}(0.427, 0.439, 0.208, 0.573),
        RGBA{N0f8}(0.537, 0.533, 0.18, 0.341),
        RGBA{N0f8}(0.51, 0.482, 0.196, 0.122),
        RGBA{N0f8}(0.545, 0.678, 0.196, 0.804),
        RGBA{N0f8}(0.208, 0.133, 0.082, 0.694),
        RGBA{N0f8}(0.298, 0.239, 0.133, 0.569),
        RGBA{N0f8}(0.082, 0.0, 0.0, 0.275),
        RGBA{N0f8}(0.0, 0.0, 0.0, 0.008),
        RGBA{N0f8}(0.518, 0.627, 0.196, 0.659),
        RGBA{N0f8}(0.427, 0.522, 0.153, 0.8),
        RGBA{N0f8}(0.412, 0.494, 0.133, 0.765),
        RGBA{N0f8}(0.239, 0.18, 0.082, 0.769),
        RGBA{N0f8}(0.196, 0.0, 0.082, 0.635),
        RGBA{N0f8}(0.0, 0.0, 0.0, 0.027),
        RGBA{N0f8}(0.082, 0.0, 0.0, 0.035),
        RGBA{N0f8}(0.196, 0.22, 0.0, 0.082),
        RGBA{N0f8}(0.082, 0.11, 0.0, 0.239),
        RGBA{N0f8}(0.208, 0.231, 0.0, 0.204),
        RGBA{N0f8}(0.133, 0.0, 0.0, 0.098),
        RGBA{N0f8}(0.133, 0.0, 0.0, 0.4),
        RGBA{N0f8}(0.0, 0.0, 0.0, 0.035),
        RGBA{N0f8}(0.796, 0.659, 0.22, 1.0), # <- Implicitly opaque
        RGBA{N0f8}(0.553, 0.698, 0.208, 1.0),
        RGBA{N0f8}(0.478, 0.561, 0.251, 1.0),
        RGBA{N0f8}(0.376, 0.314, 0.169, 1.0),
        RGBA{N0f8}(1.0, 1.0, 1.0, 1.0),
    ]
    img = PNGFiles.load("test_images/paletted_implicit_transparency.png")
    @testset "palette equality" begin
        @test img.values.parent == expected_palette
    end

    @testset begin "idempotency"
        new_path = "test_images/paletted_implicit_transparency_new.png"
        PNGFiles.save(new_path, img)
        new_img = PNGFiles.load(new_path)
        @test img == new_img
        Base.rm(new_path)
    end
end
