using ImageCore, FixedPointNumbers, Colors, ColorVectorSpace
using Test

@testset "map" begin
    @testset "clamp01" begin
        @test clamp01(0.1) === 0.1
        @test clamp01(-0.1) === 0.0
        @test clamp01(0.9) === 0.9
        @test clamp01(1.1) === 1.0
        @test clamp01(Inf) === 1.0
        @test clamp01(-Inf) === 0.0
        @test isnan(clamp01(NaN))
        @test clamp01(N0f8(0.1)) === N0f8(0.1)
        @test clamp01(N4f12(1.2)) === N4f12(1.0)
        @test clamp01(Gray(-0.2)) === Gray(0.0)
        @test clamp01(Gray(0.7)) === Gray(0.7)
        @test clamp01(Gray(N4f12(1.2))) === Gray(N4f12(1.0))
        @test clamp01(RGB(0.2,-0.2,1.8)) === RGB(0.2,0.0,1.0)
        A = [-1.2,0.4,800.3]
        f = takemap(clamp01, A)
        fA = f.(A)
        @test eltype(fA) == Float64
        @test fA == [0.0, 0.4, 1.0]
        f = takemap(clamp01, N0f8, A)
        fA = f.(A)
        @test eltype(fA) == N0f8
        @test fA == [N0f8(0), N0f8(0.4), N0f8(1)]

        A = [-1.2,0.4,800.3]
        clamp01!(A)
        @test eltype(A) == Float64
        @test A == [0, 0.4, 1]
    end

    @testset "clamp01nan" begin
        @test clamp01nan(0.1) === 0.1
        @test clamp01nan(-0.1) === 0.0
        @test clamp01nan(0.9) === 0.9
        @test clamp01nan(1.1) === 1.0
        @test clamp01nan(Inf) === 1.0
        @test clamp01nan(-Inf) === 0.0
        @test clamp01nan(NaN) === 0.0
        @test clamp01nan(N0f8(0.1)) === N0f8(0.1)
        @test clamp01nan(N4f12(1.2)) === N4f12(1.0)
        @test clamp01nan(Gray(-0.2)) === Gray(0.0)
        @test clamp01nan(Gray(0.7)) === Gray(0.7)
        @test clamp01nan(Gray(NaN32)) === Gray(0.0f0)
        @test clamp01nan(Gray(N4f12(1.2))) === Gray(N4f12(1.0))
        @test clamp01nan(RGB(0.2,-0.2,1.8)) === RGB(0.2,0.0,1.0)
        @test clamp01nan(RGB(0.2,NaN,1.8)) === RGB(0.2,0.0,1.0)
        A = [-1.2,0.4,-Inf,NaN,Inf,800.3]
        f = takemap(clamp01nan, A)
        fA = f.(A)
        @test eltype(fA) == Float64
        @test fA == Float64[0, 0.4, 0, 0, 1, 1]
        f = takemap(clamp01nan, N0f8, A)
        fA = f.(A)
        @test eltype(fA) == N0f8
        @test fA == [N0f8(0), N0f8(0.4), N0f8(0), N0f8(0), N0f8(1), N0f8(1)]

        A = [-1.2,0.4,-Inf,NaN,Inf,800.3]
        clamp01nan!(A)
        @test eltype(A) == Float64
        @test A == [0, 0.4, 0, 0, 1, 1]
    end

    @testset "scaleminmax" begin
        A = [0, 1, 100, 1000, 2000, -7]
        target = map(x->clamp(x, 0, 1), A/1000)
        for (f, tgt) in ((scaleminmax(0, 1000), target),
                          (scaleminmax(0, 1000.0), target),
                          (scaleminmax(N0f8, 0, 1000), N0f8.(target)),
                          (scaleminmax(N0f8, 0, 1000.0), N0f8.(target)),
                          (scaleminmax(Gray, 0, 1000), Gray{Float64}.(target)),
                          (scaleminmax(Gray{N0f8}, 0, 1000.0), Gray{N0f8}.(target)))
            fA = @inferred(map(f, A))
            @test fA == tgt
            @test eltype(fA) == eltype(tgt)
        end
        B = A.+10
        f = scaleminmax(10, 1010)
        @test f.(B) == target
        A = [0, 1, 100, 1000]
        target = A/1000
        for (f, tgt) in ((takemap(scaleminmax, A), target),
                          (takemap(scaleminmax, N0f8, A), N0f8.(target)),
                          (takemap(scaleminmax, Gray{N0f8}, A), Gray{N0f8}.(target)))
            fA = @inferred(map(f, A))
            @test fA == tgt
            @test eltype(fA) == eltype(tgt)
        end
        A = [Gray(-0.1),Gray(0.1)]
        f = scaleminmax(Gray, -0.1, 0.1)
        @test f.(A) == [Gray(0.0),Gray(1.0)]
        A = reinterpretc(RGB, [0.0 128.0; 255.0 0.0; 0.0 0.0])
        f = scaleminmax(RGB, 0, 255)
        @test f.(A) == [RGB(0,1.0,0), RGB(128/255,0,0)]
        f = scaleminmax(RGB{N0f8}, 0, 255)
        @test f.(A) == [RGB(0,1,0), RGB{N0f8}(128/255,0,0)]
        f = takemap(scaleminmax, A)
        @test f.(A) == [RGB(0,1.0,0), RGB(128/255,0,0)]
        f = takemap(scaleminmax, RGB{N0f8}, A)
        @test f.(A) == [RGB(0,1,0), RGB{N0f8}(128/255,0,0)]
    end

    @testset "scalesigned" begin
        A = [-100,1000]
        target = A/1000
        for (f, tgt) in ((scalesigned(1000), target),
                         (scalesigned(-1000, 0, 1000), target),
                         (scalesigned(-1000, 0, 1000.0), target),
                         (takemap(scalesigned, A), target))
            fA = f.(A)
            @test fA == tgt
            @test eltype(fA) == eltype(tgt)
        end
    end

    @testset "colorsigned" begin
        g, w, m = colorant"green1", colorant"white", colorant"magenta"
        for f in (colorsigned(),
                  colorsigned(g, m),
                  colorsigned(g, w, m))
            @test f(-1) == g
            @test f( 0) == w
            @test f( 1) == m
            @test f(-0.5) ≈ mapc(N0f8, 0.5g+0.5w)
            @test f( 0.5) ≈ mapc(N0f8, 0.5w+0.5m)
            @test f(-0.25) ≈ mapc(N0f8, 0.25g+0.75w)
            @test f( 0.75) ≈ mapc(N0f8, 0.75m+0.25w)
            @test f(Gray(0.75)) ≈ mapc(N0f8, 0.75m+0.25w)
        end
        g, w, m = RGBA(0.,1.,0.,1.), RGBA(1.,1.,1.,0.), RGBA(1.,0.,1.,1.)
        f = colorsigned(g, w, m)
        @test f(-1) == g
        @test f( 0) == w
        @test f( 1) == m
        @test f(-0.5) ≈ RGBA(0.5,1.0,0.5,0.5)
        @test f(-0.5) ≈ 0.5g+0.5w
        @test f( 0.5) ≈ 0.5w+0.5m
        @test f(-0.25) ≈ 0.25g+0.75w
        @test f( 0.75) ≈ 0.75m+0.25w
        w = RGBA{Float64}(colorant"white")
        f = colorsigned(g, m)
        @test f(-1) == g
        @test f( 0) == w
        @test f( 1) == m
        @test f(-0.5) ≈ 0.5g+0.5w
        @test f( 0.5) ≈ 0.5w+0.5m
        @test f(-0.25) ≈ 0.25g+0.75w
        @test f( 0.75) ≈ 0.75m+0.25w
    end
end

nothing
