@testset "FileIO" begin
    img = testimage("cameraman")
    ref = jpeg_decode(jpeg_encode(img))
    tmpfile = File{format"JPEG"}(joinpath(tmpdir, "tmp.jpg"))

    JpegTurbo.fileio_save(tmpfile, img)
    data = JpegTurbo.fileio_load(tmpfile)
    @test data == ref

    open(tmpfile, "w") do s
        JpegTurbo.fileio_save(s, img)
    end
    data = open(tmpfile) do s
        JpegTurbo.fileio_load(s)
    end
    @test data == ref
end
