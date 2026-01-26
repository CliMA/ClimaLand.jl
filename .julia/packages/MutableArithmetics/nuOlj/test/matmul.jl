# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

using Test
import MutableArithmetics as MA

struct CustomArray{T,N} <: AbstractArray{T,N} end

import LinearAlgebra, SparseArrays

function dot_test(x, y)
    @test MA.operate(LinearAlgebra.dot, x, y) == LinearAlgebra.dot(x, y)
    @test MA.operate(LinearAlgebra.dot, y, x) == LinearAlgebra.dot(y, x)
    @test MA.operate(*, x', y) == x' * y
    @test MA.operate(*, y', x) == y' * x
    @test MA.operate(*, LinearAlgebra.transpose(x), y) ==
          LinearAlgebra.transpose(x) * y
    @test MA.operate(*, LinearAlgebra.transpose(y), x) ==
          LinearAlgebra.transpose(y) * x
end

@testset "dot" begin
    x = [1im]
    y = [1]
    A = reshape(x, 1, 1)
    B = reshape(y, 1, 1)
    dot_test(x, x)
    dot_test(y, y)
    dot_test(A, A)
    dot_test(B, B)
    dot_test(x, y)
    dot_test(x, A)
    dot_test(x, B)
    dot_test(y, A)
    dot_test(y, B)
    dot_test(A, B)
end

@testset "promote_operation" begin
    x = [1]
    @test MA.promote_operation(*, typeof(x'), typeof(x)) == Int
    @test MA.promote_operation(*, typeof(transpose(x)), typeof(x)) == Int
end

@testset "convert" begin
    @test MA.scaling_convert(
        LinearAlgebra.UniformScaling{Int},
        LinearAlgebra.I,
    ) isa LinearAlgebra.UniformScaling
    @test MA.scaling_convert(Int, LinearAlgebra.I) === 1
    @test MA.scaling_convert(Int, 1) === 1
    @test MA.operate(
        convert,
        LinearAlgebra.UniformScaling{Int},
        LinearAlgebra.I,
    ) isa LinearAlgebra.UniformScaling
    @test MA.operate(convert, Int, LinearAlgebra.I) === 1
    @test MA.operate(convert, Int, 1) === 1
end

const EXPECTED_ERROR = string(
    "Cannot multiply a `Matrix{NoProdMutable}` with a ",
    "`Matrix{NoProdMutable}` because the sum of the product of a ",
    "`NoProdMutable` and a `NoProdMutable` could not be inferred so a ",
    "`Matrix{Union{}}` allocated to store the output of the ",
    "multiplication instead of a `Matrix{$(Int)}`.",
)

struct NoProdMutable <: MA.AbstractMutable end
function MA.promote_operation(
    ::typeof(*),
    ::Type{NoProdMutable},
    ::Type{NoProdMutable},
)
    return Int # Dummy result just to test error message
end

function unsupported_product()
    A = [NoProdMutable() for i in 1:2, j in 1:2]
    err = ErrorException(EXPECTED_ERROR)
    @test_throws err A * A
end

@testset "Errors" begin
    @testset "`promote_op` error" begin
        AT = CustomArray{Int,3}
        f(x, y) = nothing
        err = ErrorException(
            "`promote_operation($(f), $(CustomArray{Int,3}), $(CustomArray{Int,3}))` not implemented yet, please report this.",
        )
        @test_throws err MA.promote_operation(f, AT, AT)
    end

    @testset "Dimension mismatch" begin
        A = zeros(1, 1)
        B = zeros(2, 2)
        @test_throws DimensionMismatch MA.@rewrite A + B
        x = ones(1)
        y = ones(2)
        err = DimensionMismatch(
            "one array has length 1 which does not match the length of the next one, 2.",
        )
        @test_throws err MA.operate(*, x', y)
        @test_throws err MA.operate(*, LinearAlgebra.transpose(x), y)
        err = DimensionMismatch(
            "matrix A has dimensions (2,2), vector B has length 1",
        )
        @test_throws err MA.operate(*, x', B)
        a = zeros(0)
        @test iszero(@inferred MA.operate(LinearAlgebra.dot, a, a))
        @test iszero(@inferred MA.operate(*, a', a))
        @test iszero(@inferred MA.operate(*, LinearAlgebra.transpose(a), a))
        A = zeros(2)
        B = zeros(2, 1)
        err = DimensionMismatch(
            "Cannot sum or substract a matrix of axes `$(axes(B))` into matrix of axes `$(axes(A))`, expected axes `$(axes(B))`.",
        )
        @test_throws err MA.operate!(+, A, B)
        A = SparseArrays.spzeros(2)
        B = SparseArrays.spzeros(2, 1)
        err = DimensionMismatch(
            "Cannot sum or substract a matrix of axes `$(axes(B))` into matrix of axes `$(axes(A))`, expected axes `$(axes(B))`.",
        )
        @test_throws err MA.operate!(+, A, B)
        output = zeros(2)
        A = zeros(2, 1)
        B = zeros(2, 1)
        err = DimensionMismatch(
            "Cannot sum or substract matrices of axes `$(axes(A))` and `$(axes(B))` into a matrix of axes `$(axes(output))`, expected axes `$(axes(B))`.",
        )
        @test_throws err MA.operate_to!(output, +, A, B)
        output = SparseArrays.spzeros(2)
        A = SparseArrays.spzeros(2, 1)
        B = SparseArrays.spzeros(2, 1)
        err = DimensionMismatch(
            "Cannot sum or substract matrices of axes `$(axes(A))` and `$(axes(B))` into a matrix of axes `$(axes(output))`, expected axes `$(axes(B))`.",
        )
        @test_throws err MA.operate_to!(output, +, A, B)
        err = DimensionMismatch(
            "Cannot sum or substract a matrix of axes `$(axes(A))` into a matrix of axes `$(axes(output))`, expected axes `$(axes(A))`.",
        )
        @test_throws err MA.operate_to!(output, +, A)
        @test_throws err MA.operate_to!(output, -, A)
    end
    @testset "unsupported_product" begin
        unsupported_product()
    end
end

@testset "Matrix multiplication" begin
    @testset "matrix-vector product" begin
        A = [1 1 1; 1 1 1; 1 1 1]
        x = [1; 1; 1]
        y = [0; 0; 0]

        @test MA.mul(A, x) == [3; 3; 3]
        @test MA.mul_to!!(y, A, x) == [3; 3; 3] && y == [3; 3; 3]

        A = BigInt[1 1 1; 1 1 1; 1 1 1]
        x = BigInt[1; 1; 1]
        y = BigInt[0; 0; 0]

        @test MA.mutability(y, *, A, x) isa MA.IsMutable
        @test MA.mul(A, x) == BigInt[3; 3; 3]
        @test MA.mul_to!!(y, A, x) == BigInt[3; 3; 3] && y == BigInt[3; 3; 3]
        @test_throws DimensionMismatch MA.mul(BigInt[1 1; 1 1], BigInt[])
        @test MA.mul_to!!(BigInt[], BigInt[1 1; 1 1], BigInt[1; 1]) ==
              BigInt[2, 2]
        z = BigInt[0, 0]
        @test MA.mul_to!!(z, BigInt[1 1; 1 1], BigInt[1; 1]) === z
        @test z == BigInt[2, 2]

        @testset "mutability" begin
            alloc_test(() -> MA.promote_operation(*, typeof(A), typeof(x)), 0)
            alloc_test(
                () -> MA.promote_operation(
                    +,
                    typeof(y),
                    MA.promote_operation(*, typeof(A), typeof(x)),
                ),
                0,
            )
            alloc_test(
                () -> MA.promote_operation(
                    MA.add_mul,
                    typeof(y),
                    typeof(A),
                    typeof(x),
                ),
                0,
            )
            alloc_test(
                () -> MA.mutability(
                    typeof(y),
                    MA.add_mul,
                    typeof(y),
                    typeof(A),
                    typeof(x),
                ),
                0,
            )
            alloc_test(() -> MA.mutability(y, MA.add_mul, y, A, x), 0)
        end

        alloc_test(() -> MA.add_mul!!(y, A, x), BIGINT_ALLOC)
        alloc_test(
            () -> MA.operate_fallback!!(MA.IsMutable(), MA.add_mul, y, A, x),
            BIGINT_ALLOC,
        )
        alloc_test(() -> MA.operate!!(MA.add_mul, y, A, x), BIGINT_ALLOC)
        alloc_test(() -> MA.operate!(MA.add_mul, y, A, x), BIGINT_ALLOC)
        # Apparently, all allocations were on creating the buffer since this is allocation free:
        buffer = MA.buffer_for(MA.add_mul, typeof(y), typeof(A), typeof(x))
        alloc_test(() -> MA.buffered_operate!(buffer, MA.add_mul, y, A, x), 0)
    end
    @testset "matrix-matrix product" begin
        A = [1 2 3; 4 5 6; 6 8 9]
        B = [1 -1 2; -2 3 1; 2 -3 1]
        C = [one(Int) for i in 1:3, j in 1:3]

        D = [3 -4 7; 6 -7 19; 8 -9 29]
        @test MA.mul(A, B) == D
        @test MA.mul_to!!(C, A, B) == D
        @test C == D

        A = BigInt[1 2 3; 4 5 6; 6 8 9]
        B = BigInt[1 -1 2; -2 3 1; 2 -3 1]
        C = [one(BigInt) for i in 1:3, j in 1:3]

        D = BigInt[3 -4 7; 6 -7 19; 8 -9 29]
        @test MA.mul(A, B) == D
        @test MA.mul_to!!(C, A, B) == D
        @test C == D

        @test MA.mutability(C, *, A, B) isa MA.IsMutable
        @test_throws DimensionMismatch MA.mul(
            BigInt[1 1; 1 1],
            zeros(BigInt, 1, 1),
        )
        @test MA.mul_to!!(
            zeros(BigInt, 1, 1),
            BigInt[1 1; 1 1],
            zeros(BigInt, 2, 1),
        ) == zeros(BigInt, 2, 1)

        @testset "mutability" begin
            alloc_test(() -> MA.promote_operation(*, typeof(A), typeof(B)), 0)
            alloc_test(
                () -> MA.promote_operation(
                    +,
                    typeof(C),
                    MA.promote_operation(*, typeof(A), typeof(B)),
                ),
                0,
            )
            alloc_test(
                () -> MA.promote_operation(
                    MA.add_mul,
                    typeof(C),
                    typeof(A),
                    typeof(B),
                ),
                0,
            )
            alloc_test(
                () -> MA.mutability(
                    typeof(C),
                    MA.add_mul,
                    typeof(C),
                    typeof(A),
                    typeof(B),
                ),
                0,
            )
            alloc_test(() -> MA.mutability(C, MA.add_mul, C, A, B), 0)
        end

        alloc_test(() -> MA.add_mul!!(C, A, B), BIGINT_ALLOC)
        alloc_test(() -> MA.operate!!(MA.add_mul, C, A, B), BIGINT_ALLOC)
        alloc_test(() -> MA.operate!(MA.add_mul, C, A, B), BIGINT_ALLOC)
    end
end

@testset "matrix multiplication" begin
    X = ones(BigInt, 1, 1)
    M = ones(1, 1)
    C = X * M
    D = MA.operate!!(MA.add_mul, C, X, M)
    @test D == X * M + X * M
end

@testset "sub_mul" begin
    x = BigFloat[1, 1]
    A = BigFloat[2 2; 2 2]
    y = BigFloat[3, 3]
    MA.operate!!(MA.sub_mul, x, A, y)
    @test x == [-11, -11]
end

@testset "add_mul" begin
    x = BigFloat[1, 1]
    A = BigFloat[2 2; 2 2]
    y = BigFloat[3, 3]
    MA.operate!!(MA.add_mul, x, A, y)
    @test x == [13, 13]
end

struct Issue65Matrix <: AbstractMatrix{Float64}
    x::Matrix{Float64}
end

struct Issue65OneTo
    N::Int
end

Base.size(x::Issue65Matrix) = size(x.x)
Base.getindex(x::Issue65Matrix, args...) = getindex(x.x, args...)
Base.axes(x::Issue65Matrix, n) = Issue65OneTo(size(x.x, n))
Base.convert(::Type{Base.OneTo}, x::Issue65OneTo) = Base.OneTo(x.N)
Base.iterate(x::Issue65OneTo) = iterate(Base.OneTo(x.N))
Base.iterate(x::Issue65OneTo, arg) = iterate(Base.OneTo(x.N), arg)

@testset "Issue #65" begin
    x = [1.0 2.0; 3.0 4.0]
    A = Issue65Matrix(x)
    @test MA.operate(*, A, x[:, 1]) == x * x[:, 1]
    @test MA.operate(*, A, x) == x * x
end

@testset "Issue 154" begin
    X = big.([1 2; 3 4])
    c = big.([5, 6])
    MA.operate!!(MA.add_mul, X, c, c')
    @test X == [26 32; 33 40]
end

@testset "Issue 153-vector" begin
    A = big.([1 2; 3 4])
    b = big.([5, 6])
    ret = big.([0, 0])
    LinearAlgebra.mul!(ret, A, b)
    @test ret == A * b
end

@testset "Issue 153-matrix" begin
    A = big.([1 2; 3 4])
    B = big.([5 6; 7 8])
    ret = big.([0 0; 0 0])
    LinearAlgebra.mul!(ret, A, B)
    @test ret == A * B
end

@testset "Abstract eltype in matmul" begin
    # Test that we don't initialize the output with zero(T), which might not
    # exist.
    for M in (Matrix, LinearAlgebra.Diagonal)
        for T in (Any, Union{String,Int})
            x, x12, x22 = T[1, 2], T[1 2], M([1 2; 3 4])
            @test MA.operate(*, x, x') ≈ x * x'
            @test MA.operate(*, x', x) ≈ x' * x
            @test MA.operate(*, x12, x) ≈ x12 * x
            @test MA.operate(*, x22, x) ≈ x22 * x
            @test MA.operate(*, x', x22) ≈ x' * x22
            @test MA.operate(*, x12, x22) ≈ x12 * x22
            @test MA.operate(*, x22, x22) ≈ x22 * x22
            y = M([1.1 1.2; 1.3 1.4])
            @test MA.operate(*, y, x) ≈ y * x
            @test MA.operate(*, x', y) ≈ x' * y
            @test MA.operate(*, y, x12') ≈ y * x12'
            @test MA.operate(*, x12, y) ≈ x12 * y
            @test MA.operate(*, x22, y) ≈ x22 * y
            @test MA.operate(*, y, x22) ≈ y * x22
        end
    end
    for T in (Any, Union{String,Int})
        x, x12, x22 = T[1, 2], T[1 2], LinearAlgebra.LowerTriangular([1 2; 3 4])
        @test MA.operate(*, x, x') ≈ x * x'
        @test MA.operate(*, x', x) ≈ x' * x
        @test MA.operate(*, x12, x) ≈ x12 * x
        @test MA.operate(*, x22, x22) ≈ x22 * x22
        y = LinearAlgebra.LowerTriangular([1.1 1.2; 1.3 1.4])
        @test MA.operate(*, x22, y) ≈ x22 * y
        @test MA.operate(*, y, x22) ≈ y * x22
        # TODO(odow): These tests are broken because `Base` is also broken.
        # Although it fixed y * x12' in Julia v1.9.0.
        # @test_broken MA.operate(*, x22, x) ≈ x22 * x
        # @test_broken MA.operate(*, x', x22) ≈ x' * x22
        # @test_broken MA.operate(*, x12, x22) ≈ x12 * x22
        # @test_broken MA.operate(*, y, x) ≈ y * x
        # @test_broken MA.operate(*, x', y) ≈ x' * y
        # @test_broken MA.operate(*, y, x12') ≈ y * x12'
        # @test_broken MA.operate(*, x12, y) ≈ x12 * y
    end
end

@testset "Union{Int,Float64} eltype in matmul" begin
    # Test that we don't initialize the output with zero(Int), either by taking
    # the first available type in the union, or by looking at the first element
    # in the array.
    T = Union{Int,Float64}
    x, x12, x22 = T[1, 2.5], T[1 2.5], T[1 2.5; 3.5 4]
    @test MA.operate(*, x, x') == x * x'
    @test MA.operate(*, x', x) == x' * x
    @test MA.operate(*, x12, x) == x12 * x
    @test MA.operate(*, x22, x) == x22 * x
    @test MA.operate(*, x', x22) == x' * x22
    @test MA.operate(*, x12, x22) == x12 * x22
    @test MA.operate(*, x22, x22) == x22 * x22
    y = [1.1 1.2; 1.3 1.4]
    @test MA.operate(*, y, x) == y * x
    @test MA.operate(*, x', y) == x' * y
    @test MA.operate(*, y, x12') == y * x12'
    @test MA.operate(*, x12, y) == x12 * y
    @test MA.operate(*, x22, y) == x22 * y
    @test MA.operate(*, y, x22) == y * x22
end

@testset "Vector*Transpose{Vector}_issue_256" begin
    x = BigInt[1 2; 3 4]
    A = [1 2; 3 4]
    y = MA.@rewrite sum(A[i, :] * LinearAlgebra.transpose(x[i, :]) for i in 1:2)
    @test y == BigInt[10 14; 14 20]
end

struct Monomial end
LinearAlgebra.transpose(m::Monomial) = m
LinearAlgebra.adjoint(m::Monomial) = m
MA.promote_operation(::typeof(*), ::Type{Monomial}, ::Type{Monomial}) = Monomial
Base.:*(m::Monomial, ::Monomial) = m

# `Monomial` does not implement `+`, we should check that it does not prevent
# to do outer products of vectors
@testset "Vector*Transpose{Vector}_issue_256 with Monomial" begin
    m = Monomial()
    a = [m, m]
    for f in [LinearAlgebra.transpose, LinearAlgebra.adjoint]
        b = f(a)
        T = MA.promote_operation(*, typeof(a), typeof(b))
        @test T == typeof(a * b)
        @test T == typeof(MA.operate(*, a, b))
    end
end

@testset "Issue_271" begin
    A = reshape([1, 2], (2, 1))
    B = [1 2]
    C = MA.operate!!(*, A, B)
    @test A == reshape([1, 2], (2, 1))
    @test B == [1 2]
    @test C == A * B
    D = MA.operate!!(*, B, A)
    @test A == reshape([1, 2], (2, 1))
    @test B == [1 2]
    @test D == B * A
end

function test_array_sum(::Type{T}) where {T}
    x = zeros(T, 2)
    y = copy(x)
    z = copy(y)
    alloc_test(() -> MA.operate!(+, y, z), 0)
    alloc_test(() -> MA.add!!(y, z), 0)
    alloc_test(() -> MA.operate_to!(x, +, y, z), 0)
    alloc_test(() -> MA.add_to!!(x, y, z), 0)
    return
end

function test_sparse_vector_sum(::Type{T}) where {T}
    x = SparseArrays.sparsevec([1, 3], T[5, 7])
    y = copy(x)
    z = copy(y)
    # FIXME not sure what is allocating
    alloc_test_le(() -> MA.operate!(+, y, z), 200)
    alloc_test_le(() -> MA.operate!(-, y, z), 200)
    alloc_test(() -> MA.operate_to!(x, +, y, z), 0)
    alloc_test(() -> MA.operate_to!(x, -, y, z), 0)
    alloc_test(() -> MA.operate_to!(x, +, y), 0)
    alloc_test(() -> MA.operate_to!(x, -, y), 0)
    return
end

@testset "Array sum" begin
    test_array_sum(Int)
    test_sparse_vector_sum(Int)
end
