# See also test_pngsuite.jl file

@testset "Images with background" begin
    img_no_bg = PNGFiles.load("test_images/paletted_implicit_transparency.png", background=false)
    img_file_bg = PNGFiles.load("test_images/paletted_implicit_transparency.png", background=true)
    img_index_bg = PNGFiles.load("test_images/paletted_implicit_transparency.png", background=UInt8(46), expand_paletted=true)
    img_gray_bg = PNGFiles.load("test_images/paletted_implicit_transparency.png", background=Gray{N0f8}(1.0))
    img_rgb_bg = PNGFiles.load("test_images/paletted_implicit_transparency.png", background=RGB{N0f8}(1.0))
    img_gray_bg16 = PNGFiles.load("test_images/paletted_implicit_transparency.png", background=Gray{N0f16}(1.0))
    img_rgb_bg16 = PNGFiles.load("test_images/paletted_implicit_transparency.png", background=RGB{N0f16}(1.0))

    @testset "load" begin
        @testset "Background is not ignored" begin 
            @test img_file_bg[1, 1] == RGBA{N0f8}(1.0, 1.0, 1.0, 1.0)
            @test img_no_bg[1, 1] == RGBA{N0f8}(0.0, 0.0, 0.0, 0.0)
        end

        @testset "Background as index into palette" begin
            @test img_file_bg ≈ img_index_bg
        end
        @testset "Background as RGB 8 bit" begin
            @test img_file_bg ≈ img_rgb_bg
        end
        @testset "Background as Gray 8 bit" begin
            @test img_file_bg ≈ img_gray_bg
        end
        @testset "Background as RGB 16 bit" begin
            @test img_file_bg ≈ img_rgb_bg16
        end
        @testset "Background as Gray 16 bit" begin
            @test img_file_bg ≈ img_gray_bg16
        end
    end

    img_no_bg_no_gamma = PNGFiles.load("test_images/paletted_implicit_transparency.png", background=false, gamma=1.0)
    img_file_bg_no_gamma = PNGFiles.load("test_images/paletted_implicit_transparency.png", background=true, gamma=1.0)
    new_path = "test_images/paletted_implicit_transparency_bg_new.png"
    @testset "Paletted background rountripping" begin
        PNGFiles.save(new_path, img_no_bg_no_gamma, background=UInt8(46), file_gamma=1.0)
        new_img = PNGFiles.load(new_path, background=true, gamma=1.0)
        @test img_file_bg_no_gamma ≈ new_img
    end
end