# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

using Test
import MutableArithmetics as MA
import LinearAlgebra

struct DummyMutable end

function MA.promote_operation(
    ::typeof(+),
    ::Type{DummyMutable},
    ::Type{DummyMutable},
)
    return DummyMutable
end
# Test that this does not triggers any ambiguity even if there is a fallback specific to `/`.
function MA.promote_operation(
    ::Union{typeof(-),typeof(/)},
    ::Type{DummyMutable},
    ::Type{DummyMutable},
)
    return DummyMutable
end
MA.mutability(::Type{DummyMutable}) = MA.IsMutable()

# Theodorus' constant √3 (deined in global scope due to const)
Base.@irrational theodorus 1.73205080756887729353 sqrt(big(3))
@testset "promote_operation" begin
    @test MA.promote_operation(/, Rational{Int}, Rational{Int}) == Rational{Int}
    @test MA.promote_operation(-, DummyMutable, DummyMutable) == DummyMutable
    @test MA.promote_operation(/, DummyMutable, DummyMutable) == DummyMutable
    @testset "Issue #164" begin
        iπ() = MA.promote_operation(+, Int, typeof(π))
        @test iπ() == Float64
        alloc_test(iπ, 0)
        ℯbf() = MA.promote_operation(+, typeof(ℯ), BigFloat)
        @test ℯbf() == BigFloat
        # TODO this allocates as it creates the `BigFloat`
        #alloc_test(ℯbf, 0)
        bγ() = MA.promote_operation(/, Bool, typeof(Base.MathConstants.γ))
        @test bγ() == Float64
        alloc_test(bγ, 0)
        φf32() = MA.promote_operation(*, typeof(Base.MathConstants.φ), Float32)
        @test φf32() == Float32
        alloc_test(φf32, 0)
        # test user-defined Irrational
        i_theodorus() = MA.promote_operation(+, Int, typeof(theodorus))
        @test i_theodorus() == Float64
        alloc_test(i_theodorus, 0)
        # test _instantiate(::Type{S}) where {S<:Irrational} return value
        @test MA._instantiate(typeof(π)) == π
        @test MA._instantiate(typeof(MathConstants.catalan)) ==
              MathConstants.catalan
        @test MA._instantiate(typeof(theodorus)) == theodorus
    end

    for op in [real, imag]
        @test MA.promote_operation(op, ComplexF64) == Float64
        @test MA.promote_operation(op, Complex{Real}) == Real
    end
end

@testset "Errors" begin
    err = ArgumentError(
        "Cannot call `operate_to!(::$Int, +, ::$Int, ::$Int)` as objects of type `$Int` cannot be modifed to equal the result of the operation. Use `operate_to!!` instead which returns the value of the result (possibly modifying the first argument) to write generic code that also works when the type cannot be modified.",
    )
    @test_throws err MA.operate_to!(0, +, 0, 0)
    err = ArgumentError(
        "Cannot call `operate!(+, ::$Int, ::$Int)` as objects of type `$Int` cannot be modifed to equal the result of the operation. Use `operate!!` instead which returns the value of the result (possibly modifying the first argument) to write generic code that also works when the type cannot be modified.",
    )
    @test_throws err MA.operate!(+, 0, 0)
    err = ArgumentError(
        "Cannot call `buffered_operate_to!(::$Int, ::$Int, +, ::$Int, ::$Int)` as objects of type `$Int` cannot be modifed to equal the result of the operation. Use `buffered_operate_to!!` instead which returns the value of the result (possibly modifying the first argument) to write generic code that also works when the type cannot be modified.",
    )
    @test_throws err MA.buffered_operate_to!(0, 0, +, 0, 0)
    err = ArgumentError(
        "Cannot call `buffered_operate!(::$Int, +, ::$Int, ::$Int)` as objects of type `$Int` cannot be modifed to equal the result of the operation. Use `buffered_operate!!` instead which returns the value of the result (possibly modifying the first argument) to write generic code that also works when the type cannot be modified.",
    )
    @test_throws err MA.buffered_operate!(0, +, 0, 0)
    x = DummyMutable()
    err = ErrorException(
        "`operate_to!(::DummyMutable, +, ::DummyMutable, ::DummyMutable)` is not implemented yet.",
    )
    @test_throws err MA.operate_to!(x, +, x, x)
    err = ErrorException(
        "`operate!(+, ::DummyMutable, ::DummyMutable)` is not implemented yet.",
    )
    @test_throws err MA.operate!(+, x, x)
    err = ErrorException(
        "`buffered_operate_to!(::DummyMutable, ::DummyMutable, +, ::DummyMutable, ::DummyMutable)` is not implemented.",
    )
    @test_throws err MA.buffered_operate_to!(x, x, +, x, x)
    err = ErrorException(
        "`buffered_operate!(::DummyMutable, +, ::DummyMutable, ::DummyMutable)` is not implemented.",
    )
    @test_throws err MA.buffered_operate!(x, +, x, x)
end

@testset "operate" begin
    @testset "$T" for T in (Int, BigInt, Rational{Int})
        x = T(7)
        @testset "1-ary $op" for op in [+, *, gcd, lcm, copy, abs]
            a = op(x)
            b = MA.operate(op, x)
            @test a == b
            if MA.mutability(T, op, T) == MA.IsMutable()
                @test a !== b
            end
        end
        ops = [+, *, MA.add_mul, MA.sub_mul, MA.add_dot, gcd, lcm]
        @testset "4-ary $op" for op in ops
            a = op(x, x, x, x)
            b = MA.operate(op, x, x, x, x)
            @test a == b
            if MA.mutability(T, op, T, T, T, T) == MA.IsMutable()
                @test a !== b
            end
        end
        @testset "2-ary $op" for op in [-, /, div]
            a = op(x, x)
            b = MA.operate(op, x, x)
            @test a == b
            if MA.mutability(T, op, T, T) == MA.IsMutable()
                @test a !== b
            end
        end
        @testset "dot" for U in (Int, BigInt, Rational{Int}, Rational{BigInt})
            y = U(5)
            a = LinearAlgebra.dot(x, y)
            b = MA.operate(LinearAlgebra.dot, [x], [y])
            @test a == b
            if MA.mutability(
                promote_type(T, U),
                LinearAlgebra.dot,
                Vector{T},
                Vector{U},
            ) == MA.IsMutable()
                @test a !== b
            end
        end
    end
end

@testset "add_mul for BitArray" begin
    x = BigInt[0, 0]
    MA.operate!!(MA.add_mul, x, big(2), trues(2))
    @test x == BigInt[2, 2]
    MA.operate!!(MA.add_mul, x, big(3), BitVector([true, false]))
    @test x == BigInt[5, 2]
    x = BigInt[0 0; 0 0]
    MA.operate!!(MA.add_mul, x, big(2), trues(2, 2))
    @test x == BigInt[2 2; 2 2]
    MA.operate!!(MA.add_mul, x, big(3), BitArray([true false; true true]))
    @test x == BigInt[5 2; 5 5]
end

@testset "similar_array_type" begin
    @test MA.similar_array_type(BitArray{2}, Int) == Array{Int,2}
    @test MA.similar_array_type(BitArray{2}, Bool) == BitArray{2}
end

@testset "similar_array_type_Diagonal" begin
    z = zeros(2, 2)
    y = MA.operate!!(MA.add_mul, z, big(1), LinearAlgebra.I(2))
    @test y == BigFloat[1 0; 0 1]
    y = MA.operate!!(MA.add_mul, z, 2.4, LinearAlgebra.I(2))
    @test y === z
    @test y == Float64[2.4 0; 0 2.4]
    z = zeros(2, 2)
    y = MA.operate!!(MA.add_mul, z, 2.4, LinearAlgebra.Diagonal(1:2))
    @test y == LinearAlgebra.Diagonal(2.4 * (1:2))
end

@testset "unary op(::$T)" for T in (
    Float64,
    BigFloat,
    Int,
    BigInt,
    Rational{Int},
    Rational{BigInt},
)
    @test MA.operate!!(+, T(7)) == 7
    @test MA.operate!!(*, T(7)) == 7
    @test MA.operate!!(-, T(7)) == -7

    @test MA.operate_to!!(T(6), +, T(7)) == 7
    @test MA.operate_to!!(T(6), *, T(7)) == 7
    @test MA.operate_to!!(T(6), -, T(7)) == -7

    @test MA.operate!!(abs, T(7)) == 7
    @test MA.operate!!(abs, T(-7)) == 7

    @test MA.operate_to!!(T(6), abs, T(7)) == 7
    @test MA.operate_to!!(T(6), abs, T(-7)) == 7
end

@testset "Error-free mutability (issue #240)" begin
    for op in (+, -, *, /, div)
        for T in
            (Float64, BigFloat, Int, BigInt, Rational{Int}, Rational{BigInt})
            @test_nowarn MA.mutability(T, op, T, T) # should run without error
        end
    end
end

@testset "issue_271_mutability" begin
    a = 1
    x = [1; 2;;]
    y = [1 2; 3 4]
    z = [1 2; 3 4; 5 6]
    @test MA.mutability(x, *, x, x') == MA.IsNotMutable()
    @test MA.mutability(x, *, x', x) == MA.IsNotMutable()
    @test MA.mutability(x, *, x, a, x') == MA.IsNotMutable()
    @test MA.mutability(x, *, x', a, x) == MA.IsNotMutable()
    @test MA.mutability(y, *, y, y) == MA.IsMutable()
    @test MA.mutability(y, *, y, y, y) == MA.IsMutable()
    @test MA.mutability(y, *, y, a, y') == MA.IsMutable()
    @test MA.mutability(y, *, y', a, y) == MA.IsMutable()
    @test MA.mutability(y, *, a, a, y) == MA.IsMutable()
    @test MA.mutability(y, *, y, z', z) == MA.IsMutable()
    @test MA.mutability(z, *, z, z) == MA.IsNotMutable()
    @test MA.mutability(z, *, z, z, y) == MA.IsNotMutable()
end

@testset "issue_316_SubArray" begin
    y = reshape([1.0], 1, 1, 1)
    Y = view(y, :, :, 1)
    ret = reshape([1.0], 1, 1)
    ret = MA.operate!!(MA.add_mul, ret, 2.0, Y)
    @test ret == reshape([3.0], 1, 1)
    @test y == reshape([1.0], 1, 1, 1)
end

@testset "issue_318_neutral_element" begin
    a = rand(3)
    A = [rand(2, 2) for _ in 1:3]
    @test_throws DimensionMismatch MA.operate(LinearAlgebra.dot, a, A)
    y = a' * A
    @test isapprox(MA.fused_map_reduce(MA.add_mul, a', A), y)
    z = MA.operate(LinearAlgebra.dot, Int[], Int[])
    @test iszero(z) && z isa Int
    z = MA.operate(LinearAlgebra.dot, BigInt[], Int[])
    @test iszero(z) && z isa BigInt
    z = MA.operate(LinearAlgebra.dot, Int[], Float64[])
    @test iszero(z) && z isa Float64
    z = MA.operate(LinearAlgebra.dot, Matrix{Int}[], Matrix{Float64}[])
    @test iszero(z) && z isa Float64
    @test MA.fused_map_reduce(MA.add_mul, Matrix{Int}[], Float64[]) isa MA.Zero
    @test MA.fused_map_reduce(MA.add_mul, Float64[], Matrix{Int}[]) isa MA.Zero
end
