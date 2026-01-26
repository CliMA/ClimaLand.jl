invalid_imgs = [
    ("too_few_dimensions", MethodError, rand(127)),
    ("too_many_channels", AssertionError, rand(127, 257, 5)),
    ("too_many_dimensions", MethodError, rand(127, 257, 3, 1)),
]

invalid_palette_arguments = [
    ("palette_too_large", ArgumentError, IndirectArray(rand(1:257, 127, 256), rand(RGB{N0f8}, 257))),
    ("palette_index_too_few_dimensions", MethodError, IndirectArray(rand(1:4, 127), rand(RGB{N0f8}, 4))),
    ("palette_index_too_many_dimensions", ArgumentError, IndirectArray(rand(1:4, 127, 256, 1), rand(RGB{N0f8}, 4))),
    ("palette_eltype_too_wide", ArgumentError, IndirectArray(rand(1:4, 127, 256), rand(RGB{N0f16}, 4))),
    ("palette_non_RGB_colorant", ArgumentError, IndirectArray(rand(1:4, 127, 256), rand(Gray{N0f8}, 4))),
]

incompatible_background_arguments = [
    ("graya_img_rgb_background", ArgumentError, rand(127, 257, 2), "temp/PngSuite/bgbn4a08.png", RGB{N0f8}(1.0, 0.5, 0.0)),
    ("graya_img_index_background", ArgumentError, rand(127, 257, 2), "temp/PngSuite/bgbn4a08.png", 0x00),
    ("rgba_img_index_background", ArgumentError, rand(127, 257, 4), "temp/PngSuite/bgwn6a08.png", 0x00),
]

warn_background_arguments = [
    ("rgb_no_alpha_background", rand(127, 257, 3), "temp/PngSuite/basn2c08.png", RGB{N0f8}(0.0)),
    ("gray_no_alpha_background", rand(127, 257, 1), "temp/PngSuite/basn0g08.png", Gray{N0f8}(0.0)),
]

@testset "Invalid inputs" begin
    @testset "Dimensionality" begin
        for (case, exception, image) in invalid_imgs
            @testset "$(case) throws" begin
                @test_throws exception PNGFiles.save(joinpath(PNG_TEST_PATH, "test_img_err_$(case).png"),image)
                @test_throws exception open(io->PNGFiles.save(io,image), joinpath(PNG_TEST_PATH, "test_img_err_$(case).png"), "w")
            end
        end
    end
    
    @testset "Paletted" begin
        for (case, exception, image) in invalid_palette_arguments
            @testset "$(case) throws" begin
                @test_throws exception PNGFiles.save(joinpath(PNG_TEST_PATH, "test_img_err_$(case).png"),image)
                @test_throws exception open(io->PNGFiles.save(io,image), joinpath(PNG_TEST_PATH, "test_img_err_$(case).png"), "w")
            end
        end
    end

    @testset "Background" begin
        for (case, exception, image, load_path, bg) in incompatible_background_arguments
            @testset "$(case) throws save" begin
                @test_throws exception PNGFiles.save(joinpath(PNG_TEST_PATH, "test_img_err_$(case).png"), image, background=bg)
                @test_throws exception open(io->PNGFiles.save(io, image, background=bg), joinpath(PNG_TEST_PATH, "test_img_err_$(case).png"), "w")
            end
            @testset "$(case) throws load" begin
                @test_throws exception PNGFiles.load(load_path, background=bg)
            end
        end
        @testset "palette_img_index_background_not_expandded" begin
            @test_throws ArgumentError PNGFiles.load("test_images/paletted_implicit_transparency.png", background=0x00)
        end
    end

    warn_msg = "Background color for non-transparent image: PNGFiles._png_color_16_struct(0x00, 0x0000, 0x0000, 0x0000, 0x0000)"
    for (case, image, load_path, bg) in warn_background_arguments
        @testset "$(case) throws save" begin
            @test PNGFiles.save(joinpath(PNG_TEST_PATH, "test_img_err_$(case).png"), image, background=bg) === nothing
        end
        @testset "$(case) throws load" begin
            @test_logs (:warn, warn_msg) PNGFiles.load(load_path, background=bg)
        end
    end
    @testset "non_transparent_image_with_background_warns" begin
        @test_throws ArgumentError PNGFiles.load("test_images/nontransparent_image_with_background.png", background=0x00)
    end
end