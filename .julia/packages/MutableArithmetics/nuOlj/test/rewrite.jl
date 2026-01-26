# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

using LinearAlgebra, SparseArrays, Test
import MutableArithmetics as MA

@testset "Zero" begin
    z = MA.Zero()
    @test zero(MA.Zero) isa MA.Zero
    @test z + z isa MA.Zero
    @test z + 1 == 1
    @test 1 + z == 1
    @test z * z isa MA.Zero
    @test z * 1 isa MA.Zero
    @test 1 * z isa MA.Zero
    @test -z isa MA.Zero
    @test +z isa MA.Zero
    @test *(z) isa MA.Zero
    @test iszero(z)
    @test !isone(z)
end

# Test that the macro call `m` throws an error exception during pre-compilation
macro test_macro_throws(error, m)
    # See https://discourse.julialang.org/t/test-throws-with-macros-after-pr-23533/5878
    return quote
        @test_throws $(esc(error)) try
            @eval $m
        catch catched_error
            throw(catched_error.error)
        end
    end
end

function is_supported_test(T)
    @test MA.Test._is_supported(*, T, T)
    @test MA.Test._is_supported(+, T, T)
end

function error_test(x, y, z)
    err = ErrorException(
        "The curly syntax (sum{},prod{},norm2{}) is no longer supported. Expression: `sum{(i for i = 1:2)}`.",
    )
    @test_macro_throws err MA.@rewrite sum{i for i in 1:2}
    err = ErrorException("Unexpected comparison in expression `z >= y`.")
    @test_macro_throws err MA.@rewrite z >= y
    err = ErrorException("Unexpected comparison in expression `y <= z`.")
    @test_macro_throws err MA.@rewrite y <= z
    err = ErrorException("Unexpected comparison in expression `x <= y <= z`.")
    @test_macro_throws err MA.@rewrite x <= y <= z
end

using LinearAlgebra
using OffsetArrays

@testset "@rewrite with $T" for (T, supports_float) in [
    (Int, true),
    (Float64, true),
    (BigInt, true),
    (BigFloat, true),
    (Rational{Int}, true),
    (Rational{BigInt}, true),
    (DummyBigInt, false),
]
    @testset "is_supported_test" begin
        is_supported_test(T)
    end
    @testset "Error" begin
        error_test(T(1), T(2), T(3))
    end
    @testset "Scalar" begin
        MA.Test.scalar_test(T(2))
        MA.Test.scalar_test(T(3))
    end
    @testset "Quadratic" begin
        exclude = T == DummyBigInt ? ["quadratic_division"] : String[]
        MA.Test.quadratic_test(T(1), T(2), T(3), T(4), exclude = exclude)
    end
    @testset "Sparse" begin
        MA.Test.sparse_test(T(4), T(5), T[8 1 9; 4 3 1; 2 0 8])
    end
    exclude = if T == DummyBigInt
        ["broadcast_division", "matrix_vector_division"]
    elseif T <: Rational
        ["broadcast_division"]
    else
        String[]
    end
    @testset "Vector" begin
        MA.Test.array_test(T[-7, 1, 4], exclude = exclude)
        MA.Test.array_test(T[3, 2, 6], exclude = exclude)
    end
    @testset "Square Matrix" begin
        MA.Test.array_test(T[1 2; 2 4], exclude = exclude)
        MA.Test.array_test(T[2 4; -1 3], exclude = exclude)
        MA.Test.array_test(T[0 -4; 6 -5], exclude = exclude)
        MA.Test.array_test(T[5 1 9; -7 2 4; -2 -7 5], exclude = exclude)
    end
    @testset "Non-square matrix" begin
        MA.Test.array_test(T[5 1 9; -7 2 4], exclude = exclude)
    end
    @testset "Tensor" begin
        S = zeros(T, 2, 2, 2)
        S[1, :, :] = T[5 -8; 3 -7]
        S[2, :, :] = T[-2 8; 8 -1]
        MA.Test.array_test(S, exclude = exclude)
    end
    @testset "Non-array" begin
        x = T[2, 4, 3]
        MA.Test.non_array_test(x, OffsetArray(x, -length(x)))
    end
end

@testset "generator with not sum" begin
    @test MA.@rewrite(minimum(j^2 for j in 2:3)) == 4
    @test MA.@rewrite(maximum(j^2 for j in 2:3)) == 9
end

@testset "issue_76_trailing_dimensions" begin
    @testset "(4,) + (4,1)" begin
        X = [1, 2, 3, 4]
        Y = reshape([1, 2, 3, 4], (4, 1))
        @test MA.@rewrite(X + Y) == X + Y
        @test MA.@rewrite(Y + X) == Y + X
    end
    @testset "(4,) + (4,1,1)" begin
        X = [1, 2, 3, 4]
        Y = reshape([1, 2, 3, 4], (4, 1, 1))
        @test MA.@rewrite(X + Y) == X + Y
        @test MA.@rewrite(Y + X) == Y + X
    end
    @testset "(2,2) + (2,2,1)" begin
        X = [1 2; 3 4]
        Y = reshape([1, 2, 3, 4], (2, 2, 1))
        @test MA.@rewrite(X + Y) == X + Y
        @test MA.@rewrite(Y + X) == Y + X
    end
    @testset "(4,) - (4,1)" begin
        X = [1, 3, 5, 7]
        Y = reshape([1, 2, 3, 4], (4, 1))
        @test MA.@rewrite(X - Y) == X - Y
        @test MA.@rewrite(Y - X) == Y - X
    end
    @testset "(4,) - (4,1,1)" begin
        X = [1, 3, 5, 7]
        Y = reshape([1, 2, 3, 4], (4, 1, 1))
        @test MA.@rewrite(X - Y) == X - Y
        @test MA.@rewrite(Y - X) == Y - X
    end
    @testset "(2,2) - (2,2,1)" begin
        X = [1 3; 5 7]
        Y = reshape([1, 2, 3, 4], (2, 2, 1))
        @test MA.@rewrite(X - Y) == X - Y
        @test MA.@rewrite(Y - X) == Y - X
    end
end

struct _KwargRef{K,V}
    data::Dict{K,V}
end

Base.getindex(x::_KwargRef; i) = x.data[i]

@testset "test_rewrite_kw_in_ref" begin
    x = _KwargRef(Dict(i => i + 1 for i in 2:4))
    @test MA.@rewrite(sum(x[i = j] for j in 2:4)) == 12
end

@testset "dispatch_dot" begin
    # Symmetric
    x = DummyBigInt[1 2; 2 3]
    y = LinearAlgebra.Symmetric(x)
    @test MA.isequal_canonical(
        LinearAlgebra.dot(x, y),
        MA.operate(LinearAlgebra.dot, x, y),
    )
    a = @allocated LinearAlgebra.dot(x, y)
    b = @allocated MA.operate(LinearAlgebra.dot, x, y)
    @test a == b
    # Symmetric
    x = DummyBigInt[1 2; 2 3]
    y = LinearAlgebra.Hermitian(x)
    @test MA.isequal_canonical(
        LinearAlgebra.dot(x, y),
        MA.operate(LinearAlgebra.dot, x, y),
    )
    a = @allocated LinearAlgebra.dot(x, y)
    b = @allocated MA.operate(LinearAlgebra.dot, x, y)
    @test a == b
    # AbstractArray
    x = DummyBigInt[1 2; 2 3]
    y = x
    @test MA.isequal_canonical(
        LinearAlgebra.dot(x, y),
        MA.operate(LinearAlgebra.dot, x, y),
    )
    a = @allocated LinearAlgebra.dot(x, y)
    b = @allocated MA.operate(LinearAlgebra.dot, x, y)
    @test a == b
end

@testset "test_multiply_expr_MA_Zero" begin
    x = DummyBigInt(1)
    f = DummyBigInt(2)
    @test MA.@rewrite(
        f * sum(x for i in 1:0),
        move_factors_into_sums = false
    ) == MA.Zero()
    @test MA.@rewrite(
        sum(x for i in 1:0) * f,
        move_factors_into_sums = false
    ) == MA.Zero()
    @test MA.@rewrite(
        -f * sum(x for i in 1:0),
        move_factors_into_sums = false
    ) == MA.Zero()
    @test MA.@rewrite(
        sum(x for i in 1:0) * -f,
        move_factors_into_sums = false
    ) == MA.Zero()
    @test MA.@rewrite(
        (f + f) * sum(x for i in 1:0),
        move_factors_into_sums = false
    ) == MA.Zero()
    @test MA.@rewrite(
        sum(x for i in 1:0) * (f + f),
        move_factors_into_sums = false
    ) == MA.Zero()
    @test MA.@rewrite(
        -[f] * sum(x for i in 1:0),
        move_factors_into_sums = false
    ) == MA.Zero()
    @test MA.@rewrite(
        sum(x for i in 1:0) * -[f],
        move_factors_into_sums = false
    ) == MA.Zero()
    @test MA.isequal_canonical(
        MA.@rewrite(f + sum(x for i in 1:0), move_factors_into_sums = false),
        f,
    )
    @test MA.isequal_canonical(
        MA.@rewrite(sum(x for i in 1:0) + f, move_factors_into_sums = false),
        f,
    )
    @test MA.isequal_canonical(
        MA.@rewrite(-f + sum(x for i in 1:0), move_factors_into_sums = false),
        -f,
    )
    @test MA.isequal_canonical(
        MA.@rewrite(sum(x for i in 1:0) + -f, move_factors_into_sums = false),
        -f,
    )
    @test MA.isequal_canonical(
        MA.@rewrite(
            (f + f) + sum(x for i in 1:0),
            move_factors_into_sums = false
        ),
        f + f,
    )
    @test MA.isequal_canonical(
        MA.@rewrite(
            sum(x for i in 1:0) + (f + f),
            move_factors_into_sums = false
        ),
        f + f,
    )
end
