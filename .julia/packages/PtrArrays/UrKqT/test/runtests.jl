using PtrArrays
using Test
using Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(PtrArrays, deps_compat=false)
    Aqua.test_deps_compat(PtrArrays, check_extras=false)
end

@testset "Basics" begin
    x = malloc(Int, 10)
    @test x isa AbstractVector{Int}
    @test x isa DenseVector{Int}
    @test x isa PtrArray

    x .= 1:10
    @test x == 1:10

    @test_throws BoundsError x[11]
    @test_throws BoundsError x[0]
    @test_throws BoundsError x[0] = 0
    @test_throws BoundsError x[11] = 0

    @test free(x) === nothing

    y = malloc(Complex{Float64}, 4, 10)
    @test length(y) == 40
    @test size(y) == (4, 10)
    @test strides(y) == (1, 4)
    @test y isa DenseMatrix{Complex{Float64}}
    @test y isa PtrArray

    fill!(y, im)
    @test all(z -> z === 1.0im, y)
    @test count(!iszero, y) == 40
    y[4, 10] = 0
    @test y[40] == y[4, 10] == 0
    @test_throws BoundsError y[41]
    @test_throws BoundsError y[10, 4]

    @test free(y) === nothing

    @test_throws ArgumentError malloc(Vector{Int}, 10)

    # Strided array API
    z = malloc(Int, 4, 6, 10)
    @test length(z) == 240
    @test size(z) == (4, 6, 10)
    @test strides(z) == (1, 4, 24)
    @test Base.elsize(z) == sizeof(Int)
    @test z isa PtrArray
    @test z isa DenseArray{Int, 3}
    free(z)
end

@testset "zero-dimensional (#16)" begin
    s = @inferred malloc(Int32)
    @test s isa PtrArray{Int32, 0}
    s[] = 7
    @test Int32(7) === @inferred s[]
    free(s)
end

function f(x, y)
    z = malloc(Int, x)
    z .= y
    res = sum(z)
    free(z)
    res
end
@testset "Allocations" begin
    @test f(10, 1:10) == 55
    @test 0 == @allocated f(10, 1:10)
end

@testset "Invalid dimensions" begin
    @test_throws ArgumentError("invalid malloc dimensions") malloc(Int, -10, 10, 0)
    @test_throws ArgumentError("invalid malloc dimensions") malloc(Int, 1000000, 1000000, 1000000, 1000000)
    @test_throws ArgumentError("invalid malloc dimensions") malloc(Nothing, 1000000, 1000000, 1000000, 1000000)
    words_in_word_size = Sys.WORD_SIZE - Int(log2(sizeof(Int)))
    @test_throws ArgumentError("invalid malloc dimensions") malloc(Int, 2^(words_in_word_size-1))
    @test_throws ArgumentError("invalid malloc dimensions") malloc(Int, 2^(words_in_word_size+1)-1)
    @test_throws ArgumentError("invalid malloc dimensions") malloc(Int, 2^(words_in_word_size)-1)
    @test_throws ArgumentError("invalid malloc dimensions") malloc(Int, 2^(words_in_word_size+1))
    @test_throws ArgumentError("invalid malloc dimensions") malloc(Int, 2^(words_in_word_size))
    @test_throws ArgumentError("invalid malloc dimensions") malloc(UInt128, 2^3, 2^(Sys.WORD_SIZE-5))
    Sys.WORD_SIZE == 64 && @test_throws OutOfMemoryError() malloc(Int, 2^(Sys.WORD_SIZE-5)) # We could actually have enough memory on 32-bit systems
end
