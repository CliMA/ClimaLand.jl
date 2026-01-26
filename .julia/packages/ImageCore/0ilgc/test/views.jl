# some views are in colorchannels.jl
using Colors, FixedPointNumbers, ImageCore, OffsetArrays, Test
using OffsetArrays: IdentityUnitRange
using PaddedViews: filltype

@testset "rawview" begin
    a = map(N0f8, rand(3,5))
    a[2,2] = N0f8(0.5)
    v = rawview(a)
    @test v[2,2] === a[2,2].i
    v[1,3] = 0xff
    @test a[1,3] === N0f8(1)
    v[1,3] = 0x01
    @test a[1,3] === N0f8(1/255)
    s = view(a, 1:2, 1:2)
    v = rawview(s)
    @test v[2,2] === a[2,2].i
    v[2,2] = 0x0f
    @test a[2,2].i == 0x0f
    @test rawview(v) === v
end

@testset "normedview" begin
    a = rand(UInt8, 3, 5)
    a[2,2] = 0x80
    v = normedview(a)
    @test v[2,2] === N0f8(0.5)
    v[1,3] = 1
    @test a[1,3] === 0xff
    v[1,3] = 1/255
    @test a[1,3] === 0x01
    s = view(a, 1:2, 1:2)
    v = normedview(s)
    @test v[2,2] === N0f8(0.5)
    v[2,2] = 15/255
    @test a[2,2] == 0x0f
    @test normedview(v) === v
    @test normedview(N0f8, v) === v
end

@testset "StackedView" begin
    for (A, B, T) = (([1 3;2 4], [-1 -5; -2 -3], Int),
                     ([1 3;2 4], [-1.0 -5.0; -2.0 -3.0], Float64))
        V = @inferred(StackedView(A, B))
        @test eltype(V) == T
        @test size(V) == (2, 2, 2)
        @test axes(V) === (Base.OneTo(2), Base.OneTo(2), Base.OneTo(2))
        @test @inferred(V[1,1,1]) === T(1)
        @test @inferred(V[2,1,1]) === T(-1)
        @test V[1,:,:] == A
        @test V[2,:,:] == B
        @test_throws BoundsError V[0,1,1]
        @test_throws BoundsError V[3,1,1]
        @test_throws BoundsError V[2,0,1]
        @test_throws BoundsError V[2,3,1]
        @test_throws BoundsError V[1,1,0]
        @test_throws BoundsError V[1,1,3]
        V32 = @inferred(StackedView{Float32}(A, B))
        @test eltype(V32) == Float32
        @test V32[1,1,2] == Float32(3)
        V[1,2,2] = 0
        @test A[2,2] == 0
        V[2,1,2] = 11
        @test B[1,2] == 11

        V = @inferred(StackedView(A, zeroarray, B))
        @test eltype(V) == T
        @test size(V) == (3, 2, 2)
        @test axes(V) === (Base.OneTo(3), Base.OneTo(2), Base.OneTo(2))
        @test V[1,:,:] == A
        @test all(iszero, V[2,:,:])
        @test V[3,:,:] == B
        err = VERSION >= v"1.8.0-rc1" ? CanonicalIndexError : ErrorException
        @test_throws err V[2,1,1] = 7
        V32 = @inferred(StackedView{Float32}(A, zeroarray, B))
        @test eltype(V32) == Float32
        @test V32[1,1,2] == Float32(3)
    end

    # With mixed grayscale/real arrays
    a, b = Gray{N0f8}[0.1 0.2; 0.3 0.4], [0.5 0.6; 0.7 0.8]
    V = @inferred(StackedView(a, b))
    @test eltype(V) == Float64
    @test V[1,1,1] === Float64(N0f8(0.1))
    @test V[2,1,1] === 0.5
    V = @inferred(StackedView{N0f8}(a, b))
    @test eltype(V) == N0f8
    @test V[1,1,1] === N0f8(0.1)
    @test V[2,1,1] === N0f8(0.5)
    @test b[1,1] === 0.5
    V = @inferred(StackedView(a, zeroarray, b))
    @test eltype(V) == Float64
    V = @inferred(StackedView{N0f8}(a, zeroarray, b))
    @test eltype(V) == N0f8

    # With colorview
    a = [0.1 0.2; 0.3 0.4]
    b = [0.5 0.6; 0.7 0.8]
    z = zeros(2,2)  # because setindex! won't work with zeroarray
    # GrayA
    v = @inferred(colorview(GrayA, a, b))
    @test eltype(v) == GrayA{Float64}
    @test @inferred(v[2,1]) === GrayA(0.3,0.7)
    v = @inferred(colorview(GrayA{N0f8}, a, b))
    @test @inferred(v[2,1]) === GrayA{N0f8}(0.3,0.7)
    v[1,2] = GrayA(0.25, 0.25)
    @test @inferred(v[1,2]) === GrayA{N0f8}(0.25, 0.25)
    v = @inferred(colorview(GrayA{N0f8}, a, zeroarray))
    @test @inferred(v[2,1]) === GrayA{N0f8}(0.3,0)
    err = VERSION >= v"1.8.0-rc1" ? CanonicalIndexError : ErrorException
    @test_throws err (v[1,2] = GrayA(0.25, 0.25))
    # RGB
    v = @inferred(colorview(RGB{N0f8}, a, zeroarray, b))
    @test @inferred(v[2,1]) === RGB{N0f8}(0.3,0,0.7)
    v = colorview(RGB{N0f8}, a, z, b)
    @test @inferred(v[2,1]) === RGB{N0f8}(0.3,0,0.7)
    v[2,1] = RGB(0,0.9,0)
    @test @inferred(v[2,1]) === RGB{N0f8}(0,0.9,0)
    @test a[2,1] == b[2,1] == 0
    @test z[2,1] == N0f8(0.9)
    # RGBA
    v = @inferred(colorview(RGBA{N0f8}, a, zeroarray, b, b))
    @test @inferred(v[2,2]) === RGBA{N0f8}(0.4,0,0.8,0.8)
    v = colorview(RGBA{N0f8}, a, copy(z), b, copy(z))
    v[2,2] = RGBA(0.75,0.8,0.75,0.8)
    @test @inferred(v[2,2]) === RGBA{N0f8}(0.75,0.8,0.75,0.8)
    @test a[2,2] == b[2,2] == N0f8(0.75)

    @test eltype(zeroarray) == Union{}
    @test_throws ErrorException colorview(RGB{N0f8}, zeroarray, zeroarray, zeroarray)
    @test_throws DimensionMismatch StackedView(rand(2,3), rand(2,5))

    # With padding
    r = reshape([0.1,0.2], 2, 1)
    g = [false true]
    cv = colorview(RGB, paddedviews(0, r, g, zeroarray)...)
    @test cv[1,1] === RGB(0.1,0,0)
    @test cv[2,1] === RGB(0.2,0,0)
    @test cv[1,2] === RGB(0,1.0,0)
    @test cv[2,2] === RGB(0,0.0,0)
    @test axes(cv) === (Base.OneTo(2), Base.OneTo(2))

    a = [0.1 0.2; 0.3 0.4]
    b = OffsetArray([0.1 0.2; 0.3 0.4], 0:1, 2:3)
    @test_throws DimensionMismatch colorview(RGB, a, b, zeroarray)
    cv = colorview(RGB, paddedviews(0, a, b, zeroarray)...)
    @test axes(cv) == (0:2, 1:3)
    @test red.(cv[axes(a)...]) == a
    @test green.(cv[axes(b)...]) == b
    @test parent(copy(cv)) == [RGB(0,0,0)   RGB(0,0.1,0)   RGB(0,0.2,0);
                               RGB(0.1,0,0) RGB(0.2,0.3,0) RGB(0,0.4,0);
                               RGB(0.3,0,0) RGB(0.4,0,0)   RGB(0,0,0)]
    chanv = channelview(cv)
    @test @inferred(axes(chanv)) == (IdentityUnitRange(1:3), IdentityUnitRange(0:2), IdentityUnitRange(1:3))
    @test chanv[1,1,1] == 0.1
    @test chanv[2,1,2] == 0.3

    @test_throws ErrorException paddedviews(0, zeroarray, zeroarray)

    # Just to boost coverage. These methods are necessary to
    # make inference happy.
    @test_throws ErrorException ImageCore._unsafe_getindex(1, (1,1))
    @test_throws ErrorException ImageCore._unsafe_setindex!(1, (1,1), 0)
end

@testset "Multi-component colorview" begin
    r, g, b = rand(Gray{N0f16},5,5), rand(Gray{N0f16},5,5), rand(Gray{N0f16},5,5)
    A = colorview(RGB, r, g, b)
    @test A[1,1] == RGB(r[1,1], g[1,1], b[1,1])
    o = ones(2, 2)
    A = colorview(RGB, o, o, zeroarray)
    @test A[1,1] == RGB(1,1,0)
    @test size(A) == (2,2)
    @test axes(A) == (Base.OneTo(2),Base.OneTo(2))
    A = colorview(RGB, o, zeroarray, o)
    @test A[1,1] == RGB(1,0,1)
    @test size(A) == (2,2)
    @test axes(A) == (Base.OneTo(2),Base.OneTo(2))
    A = colorview(RGB, zeroarray, o, o)
    @test A[1,1] == RGB(0,1,1)
    @test size(A) == (2,2)
    @test axes(A) == (Base.OneTo(2),Base.OneTo(2))
end

@testset "mosaicviews" begin
    # only test image cases
    @testset "2D inputs" begin
        A1 = fill(Gray(1.), 2, 2)
        A2 = fill(RGB(1., 0., 0.), 3, 3)
        A3 = fill(RGB(0., 1., 0.), 3, 3)
        out = mosaic(A1, A2, A3) |> collect
        @test_reference "references/mosaicviews/2d_opaque_1.png" out by=isequal
        out = mosaic(A1, A2, A3; npad=2, fillvalue=Gray(0.), nrow=2) |> collect
        @test_reference "references/mosaicviews/2d_opaque_2.png" out by=isequal
        out = mosaic(A1, A2, A3; npad=2, fillvalue=Gray(0.), nrow=2, rowmajor=true) |> collect
        @test_reference "references/mosaicviews/2d_opaque_3.png" out by=isequal

        A1 = fill(GrayA(1.), 2, 2)
        A2 = fill(RGBA(1., 0., 0.), 3, 3)
        A3 = fill(RGBA(0., 1., 0.), 3, 3)
        out = mosaic(A1, A2, A3) |> collect
        @test_reference "references/mosaicviews/2d_transparent_1.png" out by=isequal
        out = mosaic(A1, A2, A3; npad=2, fillvalue=GrayA(0., 0.), nrow=2) |> collect
        @test_reference "references/mosaicviews/2d_transparent_2.png" out by=isequal
        out = mosaic(A1, A2, A3; npad=2, fillvalue=GrayA(0., 0.), nrow=2, rowmajor=true) |> collect
        @test_reference "references/mosaicviews/2d_transparent_3.png" out by=isequal
    end

    @testset "3D inputs" begin
        A = fill(RGB(0., 0., 0.), 2, 2, 3)
        A[:, :, 1] .= RGB(1., 0., 0.)
        A[:, :, 2] .= RGB(0., 1., 0.)
        A[:, :, 3] .= RGB(0., 0., 1.)
        out = mosaic(A) |> collect
        @test_reference "references/mosaicviews/3d_opaque_1.png" out by=isequal
        out = mosaic(A; npad=2, fillvalue=Gray(0.), nrow=2) |> collect
        @test_reference "references/mosaicviews/3d_opaque_2.png" out by=isequal
        out = mosaic(A; npad=2, fillvalue=Gray(0.), nrow=2, rowmajor=true) |> collect
        @test_reference "references/mosaicviews/3d_opaque_3.png" out by=isequal
        out = mosaic(A, A; npad=2, fillvalue=Gray(0.), nrow=2) |> collect
        @test_reference "references/mosaicviews/3d_opaque_4.png" out by=isequal

        A = fill(RGBA(0., 0., 0.), 2, 2, 3)
        A[:, :, 1] .= RGBA(1., 0., 0.)
        A[:, :, 2] .= RGBA(0., 1., 0.)
        A[:, :, 3] .= RGBA(0., 0., 1.)
        out = mosaic(A) |> collect
        @test_reference "references/mosaicviews/3d_transparent_1.png" out by=isequal
        out = mosaic(A; npad=2, fillvalue=GrayA(0., 0.), nrow=2) |> collect
        @test_reference "references/mosaicviews/3d_transparent_2.png" out by=isequal
        out = mosaic(A; npad=2, fillvalue=GrayA(0.), nrow=2, rowmajor=true) |> collect
        @test_reference "references/mosaicviews/3d_transparent_3.png" out by=isequal
        out = mosaic(A, A; npad=2, fillvalue=GrayA(0.), nrow=2) |> collect
        @test_reference "references/mosaicviews/3d_transparent_4.png" out by=isequal
    end

    @testset "4D inputs" begin
        A = fill(RGB(0., 0., 0.), 2, 2, 2, 2)
        A[1, :, 1, 1] .= RGB(1., 0., 0.)
        A[:, :, 1, 2] .= RGB(0., 1., 0.)
        A[:, :, 2, 1] .= RGB(0., 0., 1.)
        out = mosaic(A) |> collect
        @test_reference "references/mosaicviews/4d_opaque_1.png" out by=isequal
        out = mosaic(A; npad=2, fillvalue=Gray(0.), nrow=2) |> collect
        @test_reference "references/mosaicviews/4d_opaque_2.png" out by=isequal
        out = mosaic(A; npad=2, fillvalue=Gray(0.), nrow=2, rowmajor=true) |> collect
        @test_reference "references/mosaicviews/4d_opaque_3.png" out by=isequal
        out = mosaic(A, A; npad=2, fillvalue=Gray(0.), nrow=2) |> collect
        @test_reference "references/mosaicviews/4d_opaque_4.png" out by=isequal

        A = fill(RGBA(0., 0., 0.), 2, 2, 2, 2)
        A[1, :, 1, 1] .= RGBA(1., 0., 0.)
        A[:, :, 1, 2] .= RGBA(0., 1., 0.)
        A[:, :, 2, 1] .= RGBA(0., 0., 1.)
        out = mosaic(A) |> collect
        @test_reference "references/mosaicviews/4d_transparent_1.png" out by=isequal
        out = mosaic(A; npad=2, fillvalue=GrayA(0., 0.), nrow=2) |> collect
        @test_reference "references/mosaicviews/4d_transparent_2.png" out by=isequal
        out = mosaic(A; npad=2, fillvalue=GrayA(0., 0.), nrow=2, rowmajor=true) |> collect
        @test_reference "references/mosaicviews/4d_transparent_3.png" out by=isequal
        out = mosaic(A, A; npad=2, fillvalue=GrayA(0.), nrow=2) |> collect
        @test_reference "references/mosaicviews/4d_transparent_4.png" out by=isequal
    end
end

@testset "PaddedViews" begin
    # cases that don't promote array eltype:
    #   * (Number, Colorant)
    #   * (Gray, Gray)
    #   * (Gray, Color3)
    #   * (Color3, Color3)
    #   * (Color3, TransparentColor)
    #   * (TransparentColor, TransparentColor)
    @test @inferred(filltype(Float32, Gray{N0f8})) == Gray{N0f8}
    @test @inferred(filltype(Float32, RGB{N0f8})) == RGB{N0f8}
    @test @inferred(filltype(Gray{Float32}, Gray{N0f8})) === Gray{N0f8}
    @test @inferred(filltype(Gray{N0f8}, Gray{Float32})) === Gray{Float32}
    @test @inferred(filltype(Gray{Float32}, RGB{N0f8})) === RGB{N0f8}
    @test @inferred(filltype(Gray{N0f8}, RGB{Float32})) === RGB{Float32}
    @test @inferred(filltype(RGB{Float32}, RGB{N0f8})) === RGB{N0f8}
    @test @inferred(filltype(RGB{N0f8}, RGB{Float32})) === RGB{Float32}
    @test @inferred(filltype(Lab{Float32}, RGB{N0f8})) === RGB{N0f8}
    @test @inferred(filltype(RGB{N0f8}, Lab{Float32})) === Lab{Float32}

    @test @inferred(filltype(Gray{N0f8}, AGray{Float32})) === AGray{Float32}
    @test @inferred(filltype(Gray{Float32}, AGray{N0f8})) === AGray{N0f8}
    @test @inferred(filltype(Gray{N0f8}, ARGB{Float32})) === ARGB{Float32}
    @test @inferred(filltype(Gray{Float32}, ARGB{N0f8})) === ARGB{N0f8}
    @test @inferred(filltype(BGR{N0f8}, ARGB{Float32})) === ARGB{Float32}
    @test @inferred(filltype(BGR{Float32}, ARGB{N0f8})) === ARGB{N0f8}
    @test @inferred(filltype(AGray{N0f8}, ARGB{Float32})) === ARGB{Float32}
    @test @inferred(filltype(AGray{Float32}, ARGB{N0f8})) === ARGB{N0f8}

    # cases that promote both colorant type and storage type
    #   * (Colorant, Number) -- numeric arrays are treated like gray images
    #   * (Color3, Gray)
    #   * (TransparentColor, Colorant)
    @test @inferred(filltype(Gray{N0f8}, Float32)) === Gray{Float32}
    @test @inferred(filltype(Gray{Float64}, Float32)) === Gray{Float32}
    @test @inferred(filltype(RGB{N0f8}, Float32)) === RGB{Float32}
    @test @inferred(filltype(RGB{Float64}, Float32)) === RGB{Float32}
    for T in (Gray24, AGray32, RGB24, ARGB32)
        @test @inferred(filltype(T, Float32)) === @inferred(filltype(Float32, T)) === T
    end

    @test @inferred(filltype(RGB{N0f8}, Gray{Float32})) === RGB{Float32}
    @test @inferred(filltype(RGB{Float32}, Gray{N0f8})) === RGB{Float32}
    @test @inferred(filltype(Lab{Float32}, Gray{N0f8})) === Lab{Float32}
    @test @inferred(filltype(Lab{Float32}, Gray{Float64})) === Lab{Float64}

    @test @inferred(filltype(AGray{N0f8}, Gray{Float32})) === AGray{Float32}
    @test @inferred(filltype(AGray{Float32}, Gray{N0f8})) === AGray{Float32}
    @test @inferred(filltype(AGray{N0f8}, RGB{Float32})) === ARGB{Float32}
    @test @inferred(filltype(AGray{Float32}, RGB{N0f8})) === ARGB{Float32}
    @test @inferred(filltype(AGray{N0f8}, Lab{Float32})) === ALab{Float32}
    @test @inferred(filltype(AGray{Float64}, Lab{Float32})) === ALab{Float64}

    @test @inferred(filltype(ARGB{N0f8}, Gray{Float32})) === ARGB{Float32}
    @test @inferred(filltype(ARGB{Float32}, Gray{N0f8})) === ARGB{Float32}
    @test @inferred(filltype(ARGB{N0f8}, BGR{Float32})) === ABGR{Float32}
    @test @inferred(filltype(ARGB{Float32}, BGR{N0f8})) === ABGR{Float32}
    @test @inferred(filltype(ARGB{N0f8}, Lab{Float32})) === ALab{Float32}
    @test @inferred(filltype(ARGB{Float64}, Lab{Float32})) === ALab{Float64}


    A = Gray{N0f8}[Gray(0.) Gray(0.3); Gray(0.6) Gray(1.)]
    fv = Gray{N0f8}(0)
    Ap = PaddedView(fv, A, (-1:4, -1:4))
    @test eltype(Ap) == Gray{N0f8}
    @test @inferred(getindex(Ap, -1, -1)) === fv
    @test @inferred(getindex(Ap, 2, 2)) === Gray{N0f8}(1.)
    @test Ap[axes(A)...] == A

    fv = RGB{Float32}(1., 0., 0.)
    Ap = PaddedView(fv, A, (-1:4, -1:4))
    @test eltype(Ap) == RGB{Float32}
    @test @inferred(getindex(Ap, -1, -1)) === fv
    @test @inferred(getindex(Ap, 2, 2)) === RGB{Float32}(1., 1., 1.)
    @test Ap[axes(A)...] == RGB{Float32}.(A)

    fv = AGray{Float32}(0.5, 0)
    Ap = PaddedView(fv, A, (-1:4, -1:4))
    @test eltype(Ap) == AGray{Float32}
    @test @inferred(getindex(Ap, -1, -1)) === fv
    @test @inferred(getindex(Ap, 2, 2)) === AGray{Float32}(1., 1.)
    @test Ap[axes(A)...] == AGray{Float32}.(A)

    fv = ARGB{Float32}(0.5, 0., 0., 0.)
    Ap = PaddedView(fv, A, (-1:4, -1:4))
    @test eltype(Ap) == ARGB{Float32}
    @test @inferred(getindex(Ap, -1, -1)) === fv
    @test @inferred(getindex(Ap, 2, 2)) === ARGB{Float32}(1., 1., 1., 1.)
    @test Ap[axes(A)...] == ARGB{Float32}.(A)

    # order irrelevant
    A = mosaic(rand(Float32, 4, 4), rand(RGB{N0f8}, 4, 4))
    @test eltype(A) == RGB{Float32}
    A = mosaic(rand(RGB{N0f8}, 4, 4), rand(Float32, 4, 4))
    @test eltype(A) == RGB{Float32}

    A = mosaic(rand(Float32, 4, 4), rand(Gray{Float64}, 4, 4))
    @test eltype(A) == Gray{Float64}
    A = mosaic(rand(Gray{Float64}, 4, 4), rand(Float32, 4, 4))
    @test eltype(A) == Gray{Float64}

    A = mosaic(rand(Float32, 4, 4), rand(Gray{Float64}, 4, 4), rand(RGB{N0f8}, 4, 4))
    @test eltype(A) == RGB{Float64}
    A = mosaic(rand(Float32, 4, 4), rand(RGB{N0f8}, 4, 4), rand(Gray{Float64}, 4, 4))
    @test eltype(A) == RGB{Float64}

    A = mosaic(rand(Gray, 4, 4), rand(RGB, 4, 4))
    @test eltype(A) == RGB{Float64} # the filltype is always a concrete type
end

nothing
