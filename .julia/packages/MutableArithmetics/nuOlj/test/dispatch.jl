# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

using LinearAlgebra

# Tests that the calls are correctly redirected to the mutable calls
# by checking allocations
function dispatch_tests(::Type{T}) where {T}
    buffer = zero(T)
    a = one(T)
    b = one(T)
    c = one(T)
    x = convert.(T, [1, 2, 3])
    # Need to allocate 1 BigInt for the result and one for the buffer
    alloc_test(() -> MA.fused_map_reduce(MA.add_mul, x, x), 2BIGINT_ALLOC)
    alloc_test(() -> MA.fused_map_reduce(MA.add_dot, x, x), 2BIGINT_ALLOC)
    if T <: MA.AbstractMutable
        alloc_test(() -> x'x, 2BIGINT_ALLOC)
        alloc_test(() -> transpose(x) * x, 2BIGINT_ALLOC)
        alloc_test(() -> LinearAlgebra.dot(x, x), 2BIGINT_ALLOC)
    end
end

@testset "Dispatch tests" begin
    dispatch_tests(BigInt)
    # On `DummyBigInt` allocates more on previous releases of Julia
    # as it's dynamically allocated
    dispatch_tests(DummyBigInt)

    @testset "dot non-concrete vector" begin
        x = [5.0, 6.0]
        y = Vector{Union{Float64,String}}(x)
        @test MA.operate(LinearAlgebra.dot, x, y) == LinearAlgebra.dot(x, y)
        @test MA.operate(*, x', y) == x' * y
    end

    @testset "dot vector of vectors" begin
        x = [5.0, 6.0]
        z = [x, x]
        @test MA.operate(LinearAlgebra.dot, z, z) == LinearAlgebra.dot(z, z)
    end
end

@testset "*(::Real, ::Hermitian)" begin
    A = DummyBigInt[1 2; 2 3]
    B = DummyBigInt[2 4; 4 6]
    D = LinearAlgebra.Hermitian(B)
    for s in (:L, :U)
        Ah = LinearAlgebra.Hermitian(A, s)
        @test all(MA.isequal_canonical.(2 * Ah, D))
        @test all(MA.isequal_canonical.(Ah * 2, D))
    end
end

@testset "*(::AbstractMutable, ::Symmetric)" begin
    for A in ([1 2; 5 3], DummyBigInt[1 2; 5 3])
        for x in (2, DummyBigInt(2))
            for s in (:L, :U)
                # *(::AbstractMutable, ::Symmetric)
                B = LinearAlgebra.Symmetric(A, s)
                C = LinearAlgebra.Symmetric(x * A, s)
                D = x * B
                @test D isa LinearAlgebra.Symmetric
                @test MA.isequal_canonical(D, C)
                @test MA.isequal_canonical(D + D, 2 * D)
                # *(::Symmetric, ::AbstractMutable)
                B = LinearAlgebra.Symmetric(A, s)
                C = LinearAlgebra.Symmetric(A * x, s)
                D = B * x
                @test D isa LinearAlgebra.Symmetric
                @test MA.isequal_canonical(D, C)
                @test MA.isequal_canonical(D + D, 2 * D)
            end
        end
    end
end

@testset "*(::Complex, ::Hermitian)" begin
    A = BigInt[1 2; 2 3]
    B = LinearAlgebra.Hermitian(DummyBigInt.(A))
    C = 2im * A
    @test 2im * B == C
    @test C isa Matrix{Complex{BigInt}}
end

@testset "operate_to!(::Array, ::typeof(*), ::AbstractMutable, ::Array)" begin
    for A in ([1 2; 5 3], DummyBigInt[1 2; 5 3])
        for x in (2, DummyBigInt(2))
            # operate_to!(::Array, *, ::AbstractMutable, ::Array)
            B = x * A
            C = zero(B)
            D = MA.operate_to!(C, *, x, A)
            @test C === D
            @test typeof(B) == typeof(C)
            @test MA.isequal_canonical(B, C)
            # operate_to!(::Array, *, ::Array, ::AbstractMutable)
            B = A * x
            C = zero(B)
            D = MA.operate_to!(C, *, A, x)
            @test C === D
            @test typeof(B) == typeof(C)
            @test MA.isequal_canonical(B, C)
        end
    end
end

function non_mutable_sum_pr306(x)
    y = zero(eltype(x))
    for xi in x
        y += xi
    end
    return y
end

@testset "sum_with_init" begin
    x = convert(Vector{DummyBigInt}, 1:100)
    # compilation
    @allocated sum(x)
    @allocated sum(x; init = DummyBigInt(0))
    @allocated non_mutable_sum_pr306(x)
    # now test actual allocations
    no_init = @allocated sum(x)
    with_init = @allocated sum(x; init = DummyBigInt(0))
    no_ma = @allocated non_mutable_sum_pr306(x)
    # There's an additional 16 bytes for kwarg version. Upper bound by 40 to be
    # safe between Julia versions
    @test with_init <= no_init + 40
    # MA is at least 10-times better than no MA for this example
    @test 10 * with_init < no_ma
end

@testset "sum_with_init_and_dims" begin
    x = reshape(convert(Vector{DummyBigInt}, 1:12), 3, 4)
    X = reshape(1:12, 3, 4)
    for dims in (1, 2, :, 1:2, (1, 2))
        # Without (; init)
        @test MA.isequal_canonical(sum(x; dims), DummyBigInt.(sum(X; dims)))
        # With (; init)
        y = sum(x; init = DummyBigInt(0), dims)
        @test MA.isequal_canonical(y, DummyBigInt.(sum(X; dims)))
    end
end

@testset "issue_336" begin
    A = DummyBigInt[1 2; 3 4]
    B = Any[A[1], A[2]]
    @test MA.isequal_canonical(A * B, A * A[1:2])
    @test_throws DimensionMismatch A * Any[A[1]]
    @test_throws DimensionMismatch MA.operate(*, A, Any[A[1]])
end
