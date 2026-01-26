using ImageCore, Colors, FixedPointNumbers, OffsetArrays
using Test, Random

@testset "reinterpret" begin
    # Gray
    for sz in ((4,), (4,5))
        a = rand(Gray{N0f8}, sz)
        for T in (Gray{N0f8}, Gray{Float32}, Gray{Float64})
            b = @inferred(map(T, a))
            rb = @inferred(reinterpretc(eltype(T), b))
            @test eltype(rb) == eltype(T) && ndims(rb) == length(sz)
            @test size(rb) == sz
            c = copy(rb)
            rc = @inferred(reinterpretc(T, c))
            @test eltype(rc) == T && ndims(rc) == length(sz)
            @test size(rc) == sz
        end
    end
    for sz in ((4,), (4,5))
        # Bool/Gray{Bool}
        b = rand(Bool, sz)
        rb = @inferred(reinterpretc(Gray{Bool}, b))
        @test eltype(rb) == Gray{Bool} && ndims(rb) == length(sz)
        @test size(rb) == sz
        c = copy(rb)
        rc = @inferred(reinterpretc(Bool, c))
        @test eltype(rc) == Bool && ndims(rc) == length(sz)
        @test size(rc) == sz
    end
    for sz in ((4,), (4,5))
        b = Gray24.(reinterpret(N0f8, rand(UInt8, sz)))
        for T in (UInt32, RGB24)
            rb = @inferred(reinterpretc(T, b))
            @test eltype(rb) == T && ndims(rb) == length(sz)
            @test size(rb) == sz
            c = copy(rb)
            rc = @inferred(reinterpretc(Gray24, c))
            @test eltype(rc) == Gray24 && ndims(rc) == length(sz)
            @test size(rc) == sz
        end
    end
    # TransparentGray
    a = rand(AGray{N0f8}, (4,5))
    for T in (AGray{N0f8}, GrayA{Float32}, AGray{Float64})
        b = @inferred(broadcast(T, a))
        rb = @inferred(reinterpretc(eltype(T), b))
        @test eltype(rb) == eltype(T) && ndims(rb) == 3
        @test size(rb) == (2,4,5)
        c = copy(rb)
        rc = @inferred(reinterpretc(T, c))
        @test eltype(rc) == T && ndims(rc) == 2
        @test size(rc) == (4,5)
    end
    # Color3
    a = rand(RGB{N0f8}, (4,5))
    for T in (RGB{N0f8}, HSV{Float32}, XYZ{Float64})
        b = @inferred(broadcast(T, a))
        rb = @inferred(reinterpretc(eltype(T), b))
        @test eltype(rb) == eltype(T) && ndims(rb) == 3
        @test size(rb) == (3,4,5)
        c = copy(rb)
        rc = @inferred(reinterpretc(T, c))
        @test eltype(rc) == T && ndims(rc) == 2
        @test size(rc) == (4,5)
    end
    for a in (rand(RGB{N0f8}, 4), rand(RGB{N0f8}, (4,5)))
        b = @inferred(reinterpretc(HSV{Float32}, float32.(a)))
        @test eltype(b) == HSV{Float32}
        @test ndims(b) == ndims(a)
    end
    # Transparent color
    a = rand(ARGB{N0f8}, (4,5))
    for T in (ARGB{N0f8}, AHSV{Float32}, AXYZ{Float64})
        b = @inferred(broadcast(T, a))
        rb = @inferred(reinterpretc(eltype(T), b))
        @test eltype(rb) == eltype(T) && ndims(rb) == 3
        @test size(rb) == (4,4,5)
        c = copy(rb)
        rc = @inferred(reinterpretc(T, c))
        @test eltype(rc) == T && ndims(rc) == 2
        @test size(rc) == (4,5)
    end
    # XRGB/RGBX
    a = rand(RGB{N0f8}, (4,5))
    for T in (XRGB{N0f8},RGBX{Float32})
        b = @inferred(broadcast(T, a))
        rb = @inferred(reinterpretc(eltype(T), b))
        @test eltype(rb) == eltype(T) && ndims(rb) == 3
        @test size(rb) == (4,4,5)
        c = copy(rb)
        rc = @inferred(reinterpretc(T, c))
        @test eltype(rc) == T && ndims(rc) == 2
        @test size(rc) == (4,5)
    end
    a = [RGB(1,0,0) RGB(0,0,1);
         RGB(0,1,0) RGB(1,1,1)]
    @test @inferred(reinterpretc(N0f8, a)) == cat([1 0; 0 1; 0 0], [0 1; 0 1; 1 1]; dims=3)
    b = BGR{N0f8}.(a)
    @test @inferred(reinterpretc(N0f8, b)) == cat([0 0; 0 1; 1 0], [1 1; 0 1; 0 1]; dims=3)
    # RGB24, ARGB32
    for sz in ((4,), (4,5))
        a = rand(UInt32, sz)
        for T in (RGB24, ARGB32)
            b = @inferred(reinterpretc(T, a))
            @test eltype(b) == T && ndims(b) == length(sz)
            @test size(b) == sz
            @test eltype(b) == T
            @test @inferred(reinterpret(UInt32, b)) == a
        end
    end

    # 1d
    a = RGB{Float64}[RGB(1,1,0)]
    af = @inferred(reinterpretc(Float64, a))
    anew = @inferred(reinterpretc(RGB, vec(af)))
    @test anew[1] == a[1]
    @test ndims(anew) == 0

    # #33 and its converse
    a = reinterpretc(BGRA{N0f8}, [0xf0884422])
    @test isa(a, AbstractVector) && a == [BGRA{N0f8}(0.533,0.267,0.133,0.941)]
    @test reinterpretc(UInt32, a) == [0xf0884422]
    @test size(reinterpretc(BGRA{N0f8}, rand(UInt32, 5, 5))) == (5,5)
    @test size(colorview(ARGB32, rand(BGRA{N0f8}, 5, 5))) == (5,5)
    a = reinterpretc(BGRA{N0f8}, reshape([0x22, 0x44, 0x88, 0xf0, 0x01, 0x02, 0x03, 0x04], 4, 2))
    @test a == [BGRA{N0f8}(0.533,0.267,0.133,0.941), BGRA{N0f8}(0.012, 0.008, 0.004, 0.016)]
    @test reinterpret(UInt8, a) == [0x22, 0x44, 0x88, 0xf0, 0x01, 0x02, 0x03, 0x04]
    @test colorview(ARGB32, a) == reinterpretc(ARGB32, [0xf0884422,0x04030201])

    # indeterminate type tests
    a = Array{RGB{AbstractFloat}}(undef, 3)
    @test_throws Union{ArgumentError,ErrorException} reinterpretc(Float64, a)
    a = Vector{RGB}(undef, 3)
    @test_throws Union{ArgumentError,ErrorException} reinterpretc(Float64, a)

    # Invalid conversions
    a = rand(UInt8, 4,5)
    ret = @test_throws TypeError reinterpretc(Gray, a)
    a = rand(Int8, 4,5)
    ret = @test_throws TypeError reinterpretc(Gray, a)
end

@testset "eltype conversion" begin
    @test float32(Float64) == Float32
    @test float32(N0f8)      == Float32
    @test float64(RGB{N0f8}) == RGB{Float64}

    a = [RGB(1,0,0) RGB(0,0,1);
         RGB(0,1,0) RGB(1,1,1)]
    @test eltype(a) == RGB{N0f8}
    @test eltype(n0f8.(a))       == RGB{N0f8}
    @test eltype(n6f10.(a)) == RGB{N6f10}
    @test eltype(n4f12.(a)) == RGB{N4f12}
    @test eltype(n2f14.(a)) == RGB{N2f14}
    @test eltype(n0f16.(a)) == RGB{N0f16}
#    @test eltype(float16.(a)) == RGB{Float16}
    @test eltype(float32.(a)) == RGB{Float32}
    @test eltype(float64.(a)) == RGB{Float64}

    a = N0f8[0.1,0.2,0.3]
    @test eltype(a) == N0f8
    @test eltype(n0f8.(a))       == N0f8
    @test eltype(n6f10.(a)) == N6f10
    @test eltype(n4f12.(a)) == N4f12
    @test eltype(n2f14.(a)) == N2f14
    @test eltype(n0f16.(a)) == N0f16
#    @test eltype(float16.(a)) == Float16
    @test eltype(float32.(a)) == Float32
    @test eltype(float64.(a)) == Float64

    a = OffsetArray(N0f8[0.1,0.2,0.3], -1:1)
    @test eltype(a) == N0f8
    @test eltype(n0f8.(a))       == N0f8
    @test eltype(n6f10.(a)) == N6f10
    @test eltype(n4f12.(a)) == N4f12
    @test eltype(n2f14.(a)) == N2f14
    @test eltype(n0f16.(a)) == N0f16
#    @test eltype(float16.(a)) == Float16
    @test eltype(float32.(a)) == Float32
    @test eltype(float64.(a)) == Float64
    @test axes(float32.(a)) == (-1:1,)

    @testset "float" begin
        # float(::Type) is equivalent to floattype
        @test @inferred(float(RGBA{Float32})) == RGBA{Float32}
        @test @inferred(float(BGR{N0f8})    ) == BGR{Float32}
        @test @inferred(float(Gray{N0f8})   ) == Gray{Float32}

        # `floattype(::Number)` is inconsistent to `float(::NUmber)`, see issue:
        # https://github.com/JuliaMath/FixedPointNumbers.jl/issues/127
        # @test @inferred(float(N0f8)         ) == Float32
        # @test @inferred(float(Bool)         ) == Float32
        # @test @inferred(float(Float32)      ) == Float32
        # @test @inferred(float(Float64)      ) == Float64

        @test @inferred(float(ARGB32)) == ARGB{Float32}
        @test @inferred(float(AGray32)) == AGray{Float32}
        @test @inferred(float(RGB24)) == RGB{Float32}
        @test @inferred(float(Gray24)) == Gray{Float32}

        # float(x::Colorant)
        @test float(oneunit(Gray{N0f8})) == oneunit(Gray{Float32})
        @test float(RGB(oneunit(Gray{N0f8}))) == RGB(oneunit(Gray{Float32}))
    end
end

@testset "convert and broadcasting" begin
    a = [RGB(1,0,0) RGB(0,0,1);
         RGB(0,1,0) RGB(1,1,1)]
    c = @inferred(convert(Array{BGR}, a))
    @test eltype(c) === BGR
    c = @inferred(convert(Array{BGR{Float32}}, a))
    @test eltype(c) === BGR{Float32}
    c = @inferred(convert(Array{Lab}, a))
    @test eltype(c) === Lab
    for a in (rand(Float32, (4,5)),
              bitrand(4,5))
        b = @inferred(convert(Array{Gray}, a))
        @test eltype(b) === Gray
        b = @inferred(convert(Array{Gray{N0f8}}, a))
        @test eltype(b) === Gray{N0f8}
    end

    # Gray images wrapped by an OffsetArray.
    A = rand(8,8)
    for img in ( Gray.(A),
                 Gray.(N0f8.(A)),
                 Gray.(N0f16.(A)) )
        imgo = OffsetArray(img, -2, -1)
        for T in (Gray{N0f8}, Gray{N0f16}, Gray{Float32})
            s = @inferred(broadcast(T, imgo))
            @test eltype(s) == T
            @test s isa OffsetArray{T,2,Array{T,2}}
            @test permutedims(permutedims(s,(2,1)),(2,1)) == s
            @test axes(s) === axes(imgo)
        end
    end

    # Color images wrapped by an OffsetArray.
    A = rand(RGB{Float32},8,8)
    for img in ( A,
                 n0f8.(A),
                 n6f10.(A),
                 n4f12.(A),
                 n2f14.(A),
                 n0f16.(A))
        imgo = OffsetArray(img, -2, -1)
        for T in (RGB{N0f8}, RGB{Float32})
            s = @inferred(broadcast(T, imgo))
            @test eltype(s) == T
            @test s isa OffsetArray{T,2,Array{T,2}}
            @test permutedims(permutedims(s,(2,1)),(2,1)) == s
            @test axes(s) === axes(imgo)
        end
    end
end


@testset "ambiguities" begin
    # issue #40
    d = Dict(:a=>1, :b=>2.0)
    @test isa(Dict(k=>v for (k,v) in d), Dict{Symbol,Real})
end

nothing
