using ColorTypes
using FixedPointNumbers
using Downloads
using Test
using WebP

@testset "Decoding" begin
    webp_galleries = (
        lossy = (
            url = "https://www.gstatic.com/webp/gallery",
            data = Dict(
                "1.webp" => (368, 550),
                "2.webp" => (404, 550),
                "3.webp" => (720, 1280),
                "4.webp" => (772, 1024),
                "5.webp" => (752, 1024),
            ),
        ),
        lossless = (
            url = "https://www.gstatic.com/webp/gallery3",
            data = Dict(
                "1_webp_ll.webp" => (301, 400),
                "2_webp_ll.webp" => (395, 386),
                "3_webp_ll.webp" => (600, 800),
                "4_webp_ll.webp" => (163, 421),
                "5_webp_ll.webp" => (300, 300),
            ),
        ),
    )
    for gallery in webp_galleries
        for (filename, image_size) in gallery.data
            mktempdir() do tmp_dir_path
                file_path = joinpath(tmp_dir_path, filename)
                Downloads.download(joinpath(gallery.url, filename), file_path)

                for kwargs in (NamedTuple(), (transpose = true,), (transpose = false,))
                    if hasproperty(kwargs, :transpose) && kwargs.transpose
                        expected_image_size = reverse(image_size)
                    else
                        expected_image_size = image_size
                    end

                    @testset "WebP.read_webp($(joinpath(gallery.url, filename)); $kwargs)" begin
                        image = WebP.read_webp(file_path; kwargs...)
                        @test size(image) == expected_image_size
                    end

                    for TColor in [
                        ARGB{N0f8}, BGR{N0f8}, BGRA{N0f8}, RGB{N0f8}, RGBA{N0f8}, Gray{N0f8}
                    ]
                        @testset "WebP.read_webp($TColor, $(joinpath(gallery.url, filename)); $kwargs)" begin
                            image = WebP.read_webp(TColor, file_path; kwargs...)
                            @test size(image) == expected_image_size
                        end
                    end
                end
            end
        end
    end
end
