@testset "restrict" begin
    @testset "interfaces" begin
        A = rand(N0f8, 4, 5, 3)

        Ar = @inferred restrict(A, 1)
        @test typeof(Ar) <: Array
        @test size(Ar) == (3, 5, 3)

        Ar = @inferred restrict(A, (1, ))
        @test typeof(Ar) <: Array
        @test size(Ar) == (3, 5, 3)

        Ar = @inferred restrict(A, (1, 2, 3))
        @test typeof(Ar) <: Array
        @test size(Ar) == (3, 3, 2)
        @test Ar == restrict(A)

        @test_throws MethodError restrict(A, 1, 2, 3)
    end

    @testset "numerical test" begin
        A = reshape([UInt16(i) for i = 1:60], 4, 5, 3)
        B = restrict(A, (1,2))
        Btarget = cat([0.96875   4.625   5.96875;
                       2.875    10.5    12.875;
                       1.90625   5.875   6.90625],
                      [8.46875  14.625  13.46875;
                       17.875    30.5    27.875;
                       9.40625  15.875  14.40625],
                      [15.96875  24.625  20.96875;
                       32.875    50.5    42.875;
                       16.90625  25.875  21.90625], dims=3)
        @test B ≈ Btarget

        Argb = reinterpretc(RGB, reinterpret(N0f16, permutedims(A, (3,1,2))))
        B = restrict(Argb)
        Bf = permutedims(reinterpretc(eltype(eltype(B)), B), (2,3,1))
        # isapprox no longer lies, so atol is now serious
        @test isapprox(Bf, Btarget/reinterpret(one(N0f16)), atol=1e-10)

        Argba = reinterpretc(RGBA{N0f16}, reinterpret(N0f16, A))
        B = restrict(Argba)
        @test isapprox(reinterpretc(eltype(eltype(B)), B), restrict(A, (2,3))/reinterpret(one(N0f16)), atol=1e-10)

        A = reshape(1:60, 5, 4, 3)
        B = restrict(A, (1,2,3))
        @test cat([2.6015625   8.71875   6.1171875;
                   4.09375    12.875     8.78125;
                   3.5390625  10.59375   7.0546875],
                  [10.1015625  23.71875  13.6171875;
                   14.09375    32.875    18.78125;
                   11.0390625  25.59375  14.5546875], dims=3) ≈ B
    end

    @testset "singleton dimension" begin
        # https://github.com/JuliaImages/ImageTransformations.jl/issues/47
        @test restrict([1]) == [1.0]
        @test restrict([0 1 0]) == [0.25 0.25]
        @test restrict([0, 1, 0], 5) == [0, 1, 0]
    end
   
    @testset "OffsetArray" begin
        A = rand(5, 4, 3)
        Ao = OffsetArray(A, (-2,1,0))

        for (dims, offsets) in [
                (1,      (-1, 1, 0)),
                (2,      (-2, 0, 0)),
                ((1, 2), (-1, 0, 0))
            ]
            Ar = @inferred restrict(Ao, dims)
            @test typeof(Ar) <: OffsetArray
            @test Ar.offsets == offsets
            @test parent(Ar) == restrict(A, dims)
        end
    end

    @testset "SubArray" begin
        x0 = [0, 0, 1, 1, 0, 0]
        xv = @view x0[1:2:6]
        x = x0[1:2:6]
        xr = restrict(xv)
        @test xr == restrict(x)
        @test xr isa Array

        x = view(x0, IdentityUnitRange(3:5))
        xr = restrict(x)
        @test xr == restrict(collect(x))
        @test xr isa Array
    end

    @testset "FixedPoint overflow" begin
        # issue https://github.com/JuliaImages/Images.jl/issues/395
        img1 = colorview(RGB, fill(0.9, 3, 5, 5))
        img2 = colorview(RGB, fill(N0f8(0.9), 3, 5, 5))
        @test isapprox(channelview(restrict(img1)), channelview(restrict(img2)), rtol=0.01)
    end

    @testset "Array of arrays" begin
        a = Vector{Int}[[3,3,3], [2,1,7],[-11,4,2]]
        @test restrict(a) == Vector{Float64}[[2,3.5/2,6.5/2], [-5,4.5/2,5.5/2]]
    end

    @testset "various colorant" begin
        # issue https://github.com/JuliaImages/Images.jl/issues/652
        img = testimage("cameraman")
        @test eltype(@inferred(restrict(img))) == Gray{Float32}
        img = testimage("mandrill")
        @test eltype(@inferred(restrict(img))) == RGB{Float32}
        @test eltype(@inferred(restrict(Lab.(img)))) == RGB{Float32}
        img = rand(RGBA{N0f8}, 11, 11)
        @test eltype(@inferred(restrict(img))) == RGBA{Float32}
        @test eltype(@inferred(restrict(LabA.(img)))) == ARGB{Float32}

        ori = repeat(distinguishable_colors(10), inner=(1, 10))
        for T in (
            RGB, BGR, RGBX, XRGB,
            ARGB, RGBA,
            RGB24, ARGB32,
        )
            img = T.(ori)
            out = @inferred restrict(img)
            if T == RGB24
                @test eltype(out) == RGB{Float32}
            elseif T == ARGB32
                @test eltype(out) == ARGB{Float32}
            else
                @test eltype(out) <: T
            end
            ref = restrict(ori)
            @test ref ≈ RGB.(out)
        end
        for T in (Gray, AGray, GrayA, Gray24)
            img = T.(ori)
            out = @inferred restrict(img)
            if T == Gray24
                @test eltype(out) == Gray{Float32}
            else
                @test eltype(out) <: T
            end
            ref = restrict(Gray.(ori))
            @test ref ≈ Gray.(out)
        end
    end
end
