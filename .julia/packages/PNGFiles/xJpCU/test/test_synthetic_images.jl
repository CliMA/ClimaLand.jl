include("test_images/synth_images.jl")

@testset "synthetic images" begin
    for (case, image) in synth_imgs
        @testset "$(case)" begin
            expected = collect(PNGFiles._prepare_buffer(image))
            fpath = joinpath(PNG_TEST_PATH, "test_img_$(case).png")
            @testset "write" begin
                open(io->PNGFiles.save(io, image), fpath, "w") #test IO method
                @test PNGFiles.save(fpath, image) === nothing
            end
            @testset "read" begin
                global read_in_pngf = PNGFiles.load(fpath)
                @test typeof(read_in_pngf) <: AbstractMatrix
            end
            @testset "compare" begin
                @test eltype(eltype(expected)) <: N0f8 ? expected == read_in_pngf : all(expected .≈ read_in_pngf)
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
            @testset "$(case): Fpath save IO is idempotent" begin
                @test all(read_in_pngf .≈ PNGFiles.load(newpath))
                @test all(read_in_pngf .≈ open(io->PNGFiles.load(io), newpath))
            end
            open(io->PNGFiles.save(io, read_in_pngf), newpath, "w") #test IO method
            @testset "$(case): Stream save IO is idempotent" begin
                @test all(read_in_pngf .≈ PNGFiles.load(newpath))
                @test all(read_in_pngf .≈ open(io->PNGFiles.load(io), newpath))
            end
        end
    end

    for (case, func_in, image) in edge_case_imgs
        @testset "$(case)" begin
            fpath = joinpath(PNG_TEST_PATH, "test_img_$(case).png")
            @testset "write" begin
                @test PNGFiles.save(fpath, image) === nothing
            end
            @testset "read" begin
                global read_in_pngf = PNGFiles.load(fpath)
                @test read_in_pngf isa Matrix
            end
            @testset "compare" begin
                @test all(read_in_pngf .== func_in(image))
            end
            global read_in_immag = _standardize_grayness(ImageMagick.load(fpath))
            @testset "$(case): ImageMagick read type equality" begin
                @test eltype(read_in_pngf) == eltype(read_in_immag)
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
                @test imdiff(read_in_pngf, PNGFiles.load(newpath)) < 0.01
                @test imdiff(read_in_pngf, open(io->PNGFiles.load(io), newpath)) < 0.01
            end
        end
    end
end
