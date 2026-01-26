using Colors, ImageCore, OffsetArrays, FixedPointNumbers, Test
using OffsetArrays: IdentityUnitRange
using BlockArrays

# backward-compatibility to ColorTypes < v0.9 or Colors < v0.11
using ImageCore: XRGB, RGBX

struct ArrayLF{T,N} <: AbstractArray{T,N}
    A::Array{T,N}
end
Base.IndexStyle(::Type{A}) where {A<:ArrayLF} = IndexLinear()
Base.size(A::ArrayLF) = size(A.A)
Base.getindex(A::ArrayLF, i::Int) = A.A[i]
Base.setindex!(A::ArrayLF, val, i::Int) = A.A[i] = val

struct ArrayLS{T,N} <: AbstractArray{T,N}
    A::Array{T,N}
end
Base.IndexStyle(::Type{A}) where {A<:ArrayLS} = IndexCartesian()
Base.size(A::ArrayLS) = size(A.A)
Base.getindex(A::ArrayLS{T,N}, i::Vararg{Int,N}) where {T,N} = A.A[i...]
Base.setindex!(A::ArrayLS{T,N}, val, i::Vararg{Int,N}) where {T,N} = A.A[i...] = val

@testset "ChannelView" begin

@testset "grayscale" begin
    a = rand(2,3)
    @test channelview(a) === a

    a0 = [Gray(N0f8(0.2)), Gray(N0f8(0.4))]
    for (a, LI) in ((copy(a0), IndexLinear()),
                    (ArrayLF(copy(a0)), IndexLinear()),
                    (ArrayLS(copy(a0)), IndexCartesian()))
        v = @inferred(channelview(a))
        @test @inferred(IndexStyle(v)) == LI
        @test ndims(v) == 1
        @test size(v) == (2,)
        @test eltype(v) == N0f8
        @test @inferred(colorview(Gray, v)) === a
        @test parent(parent(v)) === a
        @test v[1] == N0f8(0.2)
        @test v[2] == N0f8(0.4)
        @test_throws BoundsError v[0]
        @test_throws BoundsError v[3]
        v[1] = 0.8
        @test a[1] === Gray(N0f8(0.8))
        @test_throws BoundsError (v[0] = 0.6)
        @test_throws BoundsError (v[3] = 0.6)
        c = similar(v)
        @test eltype(c) == N0f8 && ndims(c) == 1
        @test length(c) == 2
        c = similar(v, 3)
        @test eltype(c) == N0f8 && ndims(c) == 1
        @test length(c) == 3
        c = similar(v, Float32)
        @test eltype(c) == Float32 && ndims(c) == 1
        @test length(c) == 2
        c = similar(v, Float16, (5,5))
        @test eltype(c) == Float16 && ndims(c) == 2
        @test size(c) == (5,5)
    end
end

@testset "RGB, HSV, etc" begin
    for T in (RGB, BGR, XRGB, RGBX, HSV, Lab, XYZ)
        a0 = [T(0.1,0.2,0.3), T(0.4, 0.5, 0.6)]
        for a in (copy(a0),
                  ArrayLS(copy(a0)))
            v = @inferred(channelview(a))
            @test ndims(v) == 2
            @test size(v) == (3,2)
            @test eltype(v) == Float64
            if T in (RGB, HSV, Lab, XYZ)
                @test @inferred(colorview(T, v)) == a && colorview(T, v) isa typeof(a)
            else
                @test @inferred(colorview(T, v)) == a
            end
            @test v[1] == v[1,1] == 0.1
            @test v[2] == v[2,1] == 0.2
            @test v[3] == v[3,1] == 0.3
            @test v[4] == v[1,2] == 0.4
            @test v[5] == v[2,2] == 0.5
            @test v[6] == v[3,2] == 0.6
            @test_throws BoundsError v[0,1]
            @test_throws BoundsError v[4,1]
            @test_throws BoundsError v[2,0]
            @test_throws BoundsError v[2,3]
            v[2] = 0.8
            @test a[1] == T(0.1,0.8,0.3)
            v[2,1] = 0.7
            @test a[1] == T(0.1,0.7,0.3)
            @test_throws BoundsError (v[0,1] = 0.7)
            @test_throws BoundsError (v[4,1] = 0.7)
            @test_throws BoundsError (v[2,0] = 0.7)
            @test_throws BoundsError (v[2,3] = 0.7)
            c = similar(v)
            @test size(c) == (3,2) && eltype(c) == Float64
            c = similar(v, (3,4))
            @test size(c) == (3,4) && eltype(c) == Float64
            c = similar(v, Float32)
            @test size(c) == (3,2) && eltype(c) == Float32
            c = similar(v, Float16, (3,5,5))
            @test size(c) == (3,5,5) && eltype(c) == Float16
        end
    end
    a = reshape([RGB(1,0,0)])  # 0-dimensional
    v = @inferred(channelview(a))
    @test axes(v) === (Base.OneTo(3),)
    v = @inferred(channelview(a))
    @test axes(v) === (Base.OneTo(3),)
end

@testset "Gray+Alpha" begin
    for T in (AGray, GrayA)
        a = [T(0.1f0,0.2f0), T(0.3f0,0.4f0), T(0.5f0,0.6f0)]
        v = @inferred(channelview(a))
        @test @inferred(colorview(T, v)) == a
        @test ndims(v) == 2
        @test size(v) == (2,3)
        @test eltype(v) == Float32
        @test v[1] == v[1,1] == 0.1f0
        @test v[2] == v[2,1] == 0.2f0
        @test v[3] == v[1,2] == 0.3f0
        @test v[4] == v[2,2] == 0.4f0
        @test v[5] == v[1,3] == 0.5f0
        @test v[6] == v[2,3] == 0.6f0
        @test_throws BoundsError v[0,1]
        @test_throws BoundsError v[3,1]
        @test_throws BoundsError v[2,0]
        @test_throws BoundsError v[2,4]
        v[2] = 0.8
        @test a[1] == T(0.1f0,0.8f0)
        v[2,1] = 0.7
        @test a[1] == T(0.1f0,0.7f0)
        @test_throws BoundsError (v[0,1] = 0.7)
        @test_throws BoundsError (v[3,1] = 0.7)
        @test_throws BoundsError (v[2,0] = 0.7)
        @test_throws BoundsError (v[2,4] = 0.7)
        c = similar(v)
        @test eltype(c) == Float32
        @test size(c) == (2,3)
        c = similar(v, (2,4))
        @test eltype(c) == Float32
        @test size(c) == (2,4)
        c = similar(v, Float64)
        @test eltype(c) == Float64
        @test size(c) == (2,3)
        c = similar(v, Float16, (2,5,5))
        @test eltype(c) == Float16
        @test size(c) == (2,5,5)
    end
end

@testset "Alpha+RGB, HSV, etc" begin
    for T in (ARGB,
              ABGR,
              AHSV,
              ALab,
              AXYZ,
              RGBA,
              BGRA,
              LabA,
              XYZA)
        a = [T(0.1,0.2,0.3,0.4), T(0.5,0.6,0.7,0.8)]
        v = @inferred(channelview(a))
        @test ndims(v) == 2
        @test size(v) == (4,2)
        @test eltype(v) == Float64
        @test v[1] == v[1,1] == 0.1
        @test v[2] == v[2,1] == 0.2
        @test v[3] == v[3,1] == 0.3
        @test v[4] == v[4,1] == 0.4
        @test v[5] == v[1,2] == 0.5
        @test v[6] == v[2,2] == 0.6
        @test v[7] == v[3,2] == 0.7
        @test v[8] == v[4,2] == 0.8
        @test_throws BoundsError v[0,1]
        @test_throws BoundsError v[5,1]
        @test_throws BoundsError v[2,0]
        @test_throws BoundsError v[2,3]
        v[2] = 0.9
        @test a[1] == T(0.1,0.9,0.3,0.4)
        v[2,1] = 0.7
        @test a[1] == T(0.1,0.7,0.3,0.4)
        @test_throws BoundsError (v[0,1] = 0.7)
        @test_throws BoundsError (v[5,1] = 0.7)
        @test_throws BoundsError (v[2,0] = 0.7)
        @test_throws BoundsError (v[2,3] = 0.7)
        c = similar(v)
        @test eltype(c) == Float64
        @test size(c) == (4,2)
        c = similar(v, (4,4))
        @test eltype(c) == Float64
        @test size(c) == (4,4)
        c = similar(v, Float32)
        @test eltype(c) == Float32
        @test size(c) == (4,2)
        c = similar(v, Float16, (4,5,5))
        @test eltype(c) == Float16
        @test size(c) == (4,5,5)
    end

    @testset "Non-1 indices" begin
        a = OffsetArray(rand(RGB{N0f8}, 3, 5), -1:1, -2:2)
        v = @inferred(channelview(a))
        @test @inferred(axes(v)) == (IdentityUnitRange(1:3), IdentityUnitRange(-1:1), IdentityUnitRange(-2:2))
        @test @inferred(v[1,0,0]) === a[0,0].r
        a = OffsetArray(rand(Gray{Float32}, 3, 5), -1:1, -2:2)
        v = @inferred(channelview(a))
        @test @inferred(axes(v)) == (IdentityUnitRange(-1:1), IdentityUnitRange(-2:2))
        @test @inferred(v[0,0]) === gray(a[0,0])
        @test @inferred(v[5]) === gray(a[5])
        v[5] = -1
        @test v[5] === -1.0f0
    end
end

end

@testset "ColorView" begin

@testset "grayscale" begin
    a0 = [N0f8(0.2), N0f8(0.4)]
    for (a, LI) in ((copy(a0), IndexLinear()),
                    (ArrayLF(copy(a0)), IndexLinear()),
                    (ArrayLS(copy(a0)), IndexCartesian()))
        @test_throws MethodError colorview(a)
        v = @inferred(colorview(Gray, a))
        @test colorview(Gray)(a) === colorview(Gray, a)
        @test @inferred(IndexStyle(v)) == LI
        @test ndims(v) == 1
        @test size(v) == (2,)
        @test eltype(v) == Gray{N0f8}
        @test @inferred(channelview(v)) === a
        @test parent(parent(v)) === a
        @test v[1] == Gray(N0f8(0.2))
        @test v[2] == Gray(N0f8(0.4))
        @test_throws BoundsError v[0]
        @test_throws BoundsError v[3]
        v[1] = 0.8
        @test a[1] === N0f8(0.8)
        @test_throws BoundsError (v[0] = 0.6)
        @test_throws BoundsError (v[3] = 0.6)
        c = similar(v)
        @test eltype(c) == Gray{N0f8} && ndims(c) == 1
        @test length(c) == 2
        c = similar(v, 3)
        @test eltype(c) == Gray{N0f8} && ndims(c) == 1
        @test length(c) == 3
        c = similar(v, Gray{Float32})
        @test eltype(c) == Gray{Float32}
        @test length(c) == 2
        c = similar(v, Gray{Float16}, (5,5))
        @test eltype(c) == Gray{Float16}
        @test size(c) == (5,5)
        c = similar(v, Float32)
        @test isa(c, Array{Float32, 1})
    end
    # two dimensional images and linear indexing
    a0 = N0f8[0.2 0.4; 0.6 0.8]
    for (a, LI) in ((copy(a0), IndexLinear()),
                    (ArrayLF(copy(a0)), IndexLinear()),
                    (ArrayLS(copy(a0)), IndexCartesian()))
        @test_throws MethodError colorview(a)
        v = @inferred(colorview(Gray, a))
        @test colorview(Gray)(a) === colorview(Gray, a)
        @test @inferred(IndexStyle(v)) == LI
        @test ndims(v) == 2
        @test size(v) == (2,2)
        @test eltype(v) == Gray{N0f8}
        @test @inferred(channelview(v)) === a
        @test parent(parent(v)) === a
        @test v[1] == Gray(N0f8(0.2))
        @test v[2] == Gray(N0f8(0.6))
        @test_throws BoundsError v[0]
        @test_throws BoundsError v[5]
        v[1] = 0.9
        @test a[1] === N0f8(0.9)
        @test_throws BoundsError (v[0] = 0.6)
        @test_throws BoundsError (v[5] = 0.6)
    end
end

@testset "RGB, HSV, etc" begin
    for T in (RGB, BGR, XRGB, RGBX, HSV, Lab, XYZ)
        a0 = [0.1 0.2 0.3; 0.4 0.5 0.6]'
        for a in (copy(a0),
                  ArrayLS(copy(a0)))
            v = @inferred(colorview(T,a))
            @test v === colorview(T)(a)
            @test @inferred(channelview(v)) === a
            @test ndims(v) == 1
            @test size(v) == (2,)
            @test eltype(v) == T{Float64}
            @test v[1] == T(0.1,0.2,0.3)
            @test v[2] == T(0.4,0.5,0.6)
            @test_throws BoundsError v[0]
            @test_throws BoundsError v[3]
            v[2] = T(0.8, 0.7, 0.6)
            @test a == [0.1 0.2 0.3; 0.8 0.7 0.6]'
            @test_throws BoundsError (v[0] = T(0.8, 0.7, 0.6))
            @test_throws BoundsError (v[3] = T(0.8, 0.7, 0.6))
            c = similar(v)
            @test size(c) == (2,)
            c = similar(v, 4)
            @test size(c) == (4,)
            c = similar(v, T{Float32})
            @test size(c) == (2,)
            c = similar(v, T)
            @test size(c) == (2,)
            c = similar(v, T{Float16}, (5,5))
            @test size(c) == (5,5)
            c = similar(v, RGB24)
            @test eltype(c) == RGB24
            @test size(c) == size(v)
        end
    end
end

@testset "Gray+Alpha" begin
    for T in (AGray, GrayA)
        a = [0.1f0 0.2f0; 0.3f0 0.4f0; 0.5f0 0.6f0]'
        v = @inferred(colorview(T, a))
        @test colorview(T, a) === colorview(T)(a)
        @test ndims(v) == 1
        @test size(v) == (3,)
        @test eltype(v) == T{Float32}
        @test @inferred(channelview(v)) === a
        @test v[1] == T(0.1f0, 0.2f0)
        @test v[2] == T(0.3f0, 0.4f0)
        @test v[3] == T(0.5f0, 0.6f0)
        @test_throws BoundsError v[0]
        @test_throws BoundsError v[4]
        v[2] = T(0.8, 0.7)
        @test a[1,2] === 0.8f0
        @test a[2,2] === 0.7f0
        @test_throws BoundsError (v[0] = T(0.8,0.7))
        @test_throws BoundsError (v[4] = T(0.8,0.7))
        c = similar(v)
        @test eltype(c) == T{Float32}
        @test size(c) == (3,)
        c = similar(v, (4,))
        @test eltype(c) == T{Float32}
        @test size(c) == (4,)
        c = similar(v, T{Float64})
        @test eltype(c) == T{Float64}
        @test size(c) == (3,)
        c = similar(v, T{Float16}, (5,5))
        @test eltype(c) == T{Float16}
        @test size(c) == (5,5)
    end
end

@testset "Alpha+RGB, HSV, etc" begin
    a = rand(ARGB{N0f8}, 5, 5)
    vc = @inferred(channelview(a))
    @test eltype(@inferred(colorview(ARGB, vc))) == ARGB{N0f8}
    cvc = @inferred(colorview(RGBA, vc))
    @test all(cvc .== a)

    for T in (ARGB,
              ABGR,
              AHSV,
              ALab,
              AXYZ,
              RGBA,
              BGRA,
              HSVA,
              LabA,
              XYZA)
        a = [0.1 0.2 0.3 0.4; 0.5 0.6 0.7 0.8]'
        v = @inferred(colorview(T, a))
        @test colorview(T, a) === colorview(T)(a)
        @test eltype(v) == T{Float64}
        @test @inferred(channelview(v)) === a
        @test ndims(v) == 1
        @test size(v) == (2,)
        @test eltype(v) == T{Float64}
        @test v[1] == T(0.1,0.2,0.3,0.4)
        @test v[2] == T(0.5,0.6,0.7,0.8)
        @test_throws BoundsError v[0]
        @test_throws BoundsError v[3]
        v[2] = T(0.9,0.8,0.7,0.6)
        @test a[1,2] == 0.9
        @test a[2,2] == 0.8
        @test a[3,2] == 0.7
        @test a[4,2] == 0.6
        @test_throws BoundsError (v[0] = T(0.9,0.8,0.7,0.6))
        @test_throws BoundsError (v[3] = T(0.9,0.8,0.7,0.6))
        c = similar(v)
        @test eltype(c) == T{Float64}
        @test size(c) == (2,)
        c = similar(v, 4)
        @test eltype(c) == T{Float64}
        @test size(c) == (4,)
        c = similar(v, T{Float32})
        @test eltype(c) == T{Float32}
        @test size(c) == (2,)
        c = similar(v, T{Float16}, (5,5))
        @test eltype(c) == T{Float16}
        @test size(c) == (5,5)
    end

    @testset "Non-1 indices" begin
        a = OffsetArray(rand(3, 3, 5), 1:3, -1:1, -2:2)
        v = @inferred(colorview(RGB, a))
        @test colorview(RGB, a) === colorview(RGB)(a)
        @test @inferred(axes(v)) == (IdentityUnitRange(-1:1), IdentityUnitRange(-2:2))
        @test @inferred(v[0,0]) === RGB(a[1,0,0], a[2,0,0], a[3,0,0])
        a = OffsetArray(rand(3, 3, 5), 0:2, -1:1, -2:2)
        @test_throws ArgumentError colorview(RGB, a)
    end

    @testset "Custom/divergent axis types" begin
        img1 = rand(5, 4, 2)
        img2_2 = mortar(reshape([rand(5, 4, 1), rand(5, 4, 1)], 1, 1, 2))
        img2_all = mortar(reshape([rand(5, 4, 1), rand(5, 4, 1), rand(5, 4, 1)], 1, 1, 3))
        img2_odd = img2_all[:,:,1:2:end]
        for img2 in (img2_2, img2_odd)
            for imgrgb in (colorview(RGB, img1, img2, zeroarray),
                        colorview(RGB, img2, img1, zeroarray))
                @test eltype(imgrgb) === RGB{Float64}
                @test size(imgrgb) == size(img1)
            end
        end
    end
end

end

nothing
