module OneZero

using Test
using VectorInterface

const Z = Zero()
const I = One()
const typelist = (Int32, Int64, Float16, Float32, Float64, ComplexF16, ComplexF32,
                  ComplexF64, BigFloat)

@testset "equalities" begin
    @test I == 1
    @test Z == 0
    @test I != Z
    @test I == I
    @test Z == Z

    @test iszero(Z) == true
    @test isone(I) == true
    @test isone(Z) == false
    @test iszero(I) == false
end

@testset "arithmetic" begin
    @test Z + Z === Z
    @test I + Z === I
    @test Z + I === I
    @test I + I == 2

    @test @inferred(-Z) == Z
    @test @inferred(-I) == -1
    @test @inferred(I - I) === Z
    @test @inferred(I - Z) === I
    @test @inferred(Z - I) == -I

    @test @inferred(Z * Z) === Z
    @test @inferred(I * Z) === Z
    @test @inferred(Z * I) === Z
    @test @inferred(I * I) === I

    @test @inferred(Z / I) === Z
    @test @inferred(I / I) === I
    @test_throws DivideError @inferred(I / Z)
    @test_throws DivideError @inferred(Z / Z)

    for T in typelist
        x = rand(T)
        while iszero(x)
            x = rand(T)
        end

        @test @inferred(Z + x) == x
        @test @inferred(x + Z) == x
        @test @inferred(I + x) == x + 1
        @test @inferred(x + I) == x + 1

        @test @inferred(Z - x) == -x
        @test @inferred(x - Z) == x
        @test @inferred(I - x) == 1 - x
        @test @inferred(x - I) == x - 1

        @test @inferred(Z * x) == zero(x)
        @test @inferred(x * Z) == zero(x)
        @test @inferred(I * x) == x
        @test @inferred(x * I) == x

        @test_throws DivideError @inferred(x / Z)
        @test @inferred(x / I) == x
        @test @inferred(Z / x) == zero(x)
        @test @inferred(I / x) == inv(x)
    end
end

@testset "promotion" begin
    for T in typelist
        @test @inferred(promote_type(typeof(I), T)) == T
        @test @inferred(promote_type(typeof(Z), T)) == T
        @test @inferred(promote_type(T, typeof(I))) == T
        @test @inferred(promote_type(T, typeof(Z))) == T

        @test @inferred(promote_type(One, Zero, T)) == T
        @test @inferred(promote_type(One, T, Zero)) == T
        @test @inferred(promote_type(T, One, Zero)) == T
        @test @inferred(promote_type(Zero, One, T)) == T
        @test @inferred(promote_type(Zero, T, One)) == T
        @test @inferred(promote_type(T, Zero, One)) == T

        @test @inferred(T(I)) == one(T)
        @test @inferred(T(Z)) == zero(T)
        @test @inferred(convert(T, I)) == one(T)
        @test @inferred(convert(T, Z)) == zero(T)
    end
end

end
