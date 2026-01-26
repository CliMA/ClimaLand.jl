@testset "dpi" begin
    io_none = IOBuffer()
    io_300 = IOBuffer()
    io_300_300 = IOBuffer()
    io_300_500 = IOBuffer()

    img = rand(RGB{N0f8}, 2, 2)

    PNGFiles.save(io_none, img; dpi = nothing)
    PNGFiles.save(io_300, img; dpi = 300)
    PNGFiles.save(io_300_300, img; dpi = (300f0, 300.0))
    PNGFiles.save(io_300_500, img; dpi = (300, 500))

    for io in [io_none, io_300, io_300_300, io_300_500]
        seekstart(io)
        @test PNGFiles.load(io) == img
    end

    s_none = String(take!(io_none))
    s_300 = String(take!(io_300))
    s_300_300 = String(take!(io_300_300))
    s_300_500 = String(take!(io_300_500))

    function physblock(xdpi, ydpi)
        io = IOBuffer()
        write(io, "pHYs")
        write(io, hton(round(UInt32, xdpi / 0.0254)))
        write(io, hton(round(UInt32, ydpi / 0.0254)))
        write(io, hton(UInt8(1)))
        return String(take!(io))
    end
    
    @test !occursin("pHYs", s_none)
    @test occursin(physblock(300, 300), s_300)
    @test occursin(physblock(300, 300), s_300_300)
    @test occursin(physblock(300, 500), s_300_500)
end
