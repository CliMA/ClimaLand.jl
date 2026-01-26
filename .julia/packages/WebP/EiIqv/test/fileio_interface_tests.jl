using FileIO
using Test
using TestImages
using WebP

@testset "FileIO interface" begin
    expected_image = testimage("lighthouse")

    mktempdir() do tmp_dir_path
        file_path = joinpath(tmp_dir_path, "lighthouse.webp")

        @testset "File{format\"WebP\"}" begin
            f = File{format"WebP"}(file_path)

            @testset "fileio_load" begin
                WebP.write_webp(file_path, expected_image)
                image = WebP.fileio_load(f)
                @test size(image) == size(expected_image)
            end

            @testset "fileio_save" begin
                WebP.fileio_save(f, expected_image)
                image = WebP.read_webp(file_path)
                @test size(image) == size(expected_image)
            end
        end

        @testset "Stream{format\"WebP\"}" begin
            s = Stream{format"WebP"}(IOBuffer())

            @testset "fileio_load" begin
                WebP.write_webp(file_path, expected_image)
                open(file_path, "r") do io
                    s = Stream{format"WebP"}(io)
                    image = WebP.fileio_load(s)
                    @test size(image) == size(expected_image)
                end
            end

            @testset "fileio_save" begin
                open(file_path, "w") do io
                    s = Stream{format"WebP"}(io)
                    WebP.fileio_save(s, expected_image)
                end
                image = WebP.read_webp(file_path)
                @test size(image) == size(expected_image)
            end
        end
    end
end
