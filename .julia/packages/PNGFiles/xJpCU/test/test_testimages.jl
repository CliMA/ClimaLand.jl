using TestImages

real_imgs = [
    first(splitext(img_name)) => testimage(img_name)
    for img_name
    in TestImages.remotefiles
    if endswith(img_name, ".png")
]

@testset "TestImages.jl" begin
    for (case, image) in real_imgs
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
                @test all(expected .≈ read_in_pngf)
            end
            global read_in_immag = _standardize_grayness(ImageMagick.load(fpath))
            @testset "$(case): ImageMagick read type equality" begin
                # The lena image is Grayscale saved as RGB...
                @test eltype(_standardize_grayness(read_in_pngf)) == eltype(read_in_immag)
            end
            @testset "$(case): ImageMagick read values equality" begin
                # The lena image is Grayscale saved as RGB...
                imdiff_val = imdiff(_standardize_grayness(read_in_pngf), read_in_immag)
                onfail(@test imdiff_val < 0.01) do
                    PNGFiles._inspect_png_read(fpath)
                    _add_debugging_entry(fpath, case, imdiff_val)
                end
            end
            path, ext = splitext(fpath)
            newpath = path * "_new" * ext
            PNGFiles.save(newpath, read_in_pngf)
            @testset "$(case): IO is idempotent" begin
                @test all(read_in_pngf .≈ PNGFiles.load(newpath))
                @test all(read_in_pngf .≈ open(io->PNGFiles.load(io), newpath))
            end
        end
    end
end
