@testset "multithreads" begin
    @test Threads.nthreads() > 1

    img = testimage("lighthouse")
    @testset "jpeg_encode" begin
        ref = jpeg_encode(img)

        out = [fill(zero(UInt8), size(ref)) for _ in 1:Threads.nthreads()]
        Threads.@threads for i in 1:Threads.nthreads()
            out[i] = jpeg_encode(img)
        end
        @test all(out .== Ref(ref))
    end

    @testset "jpeg_decode" begin
        out = [similar(img) for _ in 1:Threads.nthreads()]
        tmpdir = mktempdir()

        for i in 1:Threads.nthreads()
            jpeg_encode(joinpath(tmpdir, "$i.jpg"), img)
        end
        ref = jpeg_decode(joinpath(tmpdir, "1.jpg"))
        Threads.@threads for i in 1:Threads.nthreads()
            out[i] = jpeg_decode(joinpath(tmpdir, "$i.jpg"))
        end
        @test all(out .== Ref(ref))
    end
end
