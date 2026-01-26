using ImageCore
using FFTW
using Logging
using Test

@testset "functions" begin
    ag = rand(Gray{Float32}, 4, 5)
    ac = rand(RGB{Float32}, 4, 5)
    for (f, args) in ((fft, (ag,)), (fft, (ag, 1:2)), (plan_fft, (ag,)),
                      (rfft, (ag,)), (rfft, (ag, 1:2)), (plan_rfft, (ag,)),
                      (fft, (ac,)), (fft, (ac, 1:2)), (plan_fft, (ac,)),
                      (rfft, (ac,)), (rfft, (ac, 1:2)), (plan_rfft, (ac,)))
        dims_str = eltype(args[1])<:Gray ? "1:2" : "2:3"
        @test_throws "channelview, and likely $dims_str" f(args...)
    end
    for (a, dims) in ((ag, 1:2), (ac, 2:3))
        @test ifft(fft(channelview(a), dims), dims) ≈ channelview(a)
        ret = @test_throws "channelview, and likely $dims" rfft(a)
        @test irfft(rfft(channelview(a), dims), 4, dims) ≈ channelview(a)
    end
    # Ensure that the hint doesn't error
    @test_logs min_level=Logging.Warn try
        [3](1, 2)
    catch err
        showerror(devnull, err)
    end

    a = [RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)]
    @test a' == [RGB(1,0,0) RGB(0,1,0) RGB(0,0,1)]

    a = [RGB(1,0,0) RGB(0,1,0); RGB(0,0,1) RGB(0, 0, 0)]
    @test a' == [RGB(1,0,0) RGB(0,0,1); RGB(0,1,0) RGB(0, 0, 0)]
end

nothing
