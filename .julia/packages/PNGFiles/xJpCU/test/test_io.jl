using TestImages, Base64

test_name, test_image = first(
    first(splitext(img_name)) => testimage(img_name)
    for img_name
    in TestImages.remotefiles
    if endswith(img_name, ".png")
)


@testset "IOBuffer" begin
    @test test_image == let
        b = IOBuffer()
        PNGFiles.save(b, test_image)
        seek(b, 0)
        PNGFiles.load(b)
    end
end

@testset "Base64EncodePipe" begin
    @test test_image == let
        io = Base64EncodePipe(IOBuffer())
        PNGFiles.save(io, test_image)
        close(io)
        decoded = IOBuffer(base64decode(take!(io.io)))
        PNGFiles.load(decoded)
    end
end

