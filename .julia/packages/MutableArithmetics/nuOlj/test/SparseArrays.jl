# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

module TestInterfaceSparseArrays

using Test

import LinearAlgebra
import MutableArithmetics as MA
import Random
import SparseArrays

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_spmatmul()
    Random.seed!(1234)
    for m in [1, 2, 3, 5, 11]
        for n in [1, 2, 3, 5, 11]
            A = SparseArrays.sprand(Float64, m, n, 0.5)
            B = SparseArrays.sprand(Float64, n, m, 0.5)
            ret = SparseArrays.spzeros(Float64, m, m)
            MA.operate!(MA.add_mul, ret, A, B)
            @test ret ≈ A * B
            ret = SparseArrays.spzeros(Float64, m, m)
            MA.operate!(MA.add_mul, ret, A, A')
            @test ret ≈ A * A'
            ret = SparseArrays.spzeros(Float64, m, m)
            MA.operate!(MA.add_mul, ret, A, B, 2.0)
            @test ret ≈ A * B * 2.0
            ret = SparseArrays.spzeros(Float64, m, m)
            MA.operate!(MA.add_mul, ret, A, B, 2.0, 1.5)
            @test ret ≈ A * B * 2.0 * 1.5
        end
    end
    return
end

function test_spmatmul_prefer_sort()
    Random.seed!(1234)
    m = n = 100
    p = 0.01
    A = SparseArrays.sprand(Float64, m, n, p)
    B = SparseArrays.sprand(Float64, n, m, p)
    ret = SparseArrays.spzeros(Float64, m, m)
    MA.operate!(MA.add_mul, ret, A, B)
    @test ret ≈ A * B
    ret = SparseArrays.spzeros(Float64, m, m)
    MA.operate!(MA.add_mul, ret, A, B, 2.0)
    @test ret ≈ A * B * 2.0
    ret = SparseArrays.spzeros(Float64, m, m)
    MA.operate!(MA.add_mul, ret, A, B, 2.0, 1.5)
    @test ret ≈ A * B * 2.0 * 1.5
    return
end

function test_sparse_addmul_matmat()
    Random.seed!(1234)
    m = 10
    ret = SparseArrays.spzeros(Float64, m, m)
    A, B = rand(m, m), rand(m, m)
    MA.operate!(MA.add_mul, ret, A, B)
    @test ret ≈ A * B
    MA.operate!(MA.add_mul, ret, A, B)
    @test ret ≈ A * B + A * B
    return
end

function test_sparse_addmul_mat_adjspmat()
    Random.seed!(1234)
    m = 10
    ret = zeros(m, m)
    A = rand(m, m)
    B = SparseArrays.sprand(Float64, m, m, 0.5)'
    MA.operate!(MA.add_mul, ret, A, B)
    @test ret ≈ A * B
    MA.operate!(MA.add_mul, ret, A, B, 3.0)
    @test ret ≈ A * B + A * B * 3.0
    return
end

function test_spmatmul_adj()
    Random.seed!(1234)
    for m in [1, 2, 3, 5, 11]
        for f in (identity, LinearAlgebra.adjoint)
            for g in (identity, LinearAlgebra.adjoint)
                A = f(SparseArrays.sprand(Float64, m, m, 0.5))
                B = g(SparseArrays.sprand(Float64, m, m, 0.5))
                ret = SparseArrays.spzeros(Float64, m, m)
                MA.operate!(MA.add_mul, ret, A, B)
                @test ret ≈ A * B
                ret = SparseArrays.spzeros(Float64, m, m)
                MA.operate!(MA.add_mul, ret, A, A')
                @test ret ≈ A * A'
                ret = SparseArrays.spzeros(Float64, m, m)
                MA.operate!(MA.add_mul, ret, A, B, 2.0)
                @test ret ≈ A * B * 2.0
                ret = SparseArrays.spzeros(Float64, m, m)
                MA.operate!(MA.add_mul, ret, A, B, 2.0, 1.5)
                @test ret ≈ A * B * 2.0 * 1.5
            end
        end
    end
    return
end

function test_sparse_broadcast()
    Random.seed!(1234)
    x = [BigInt(1) for _ in 1:10]
    y = SparseArrays.sprand(Float64, 10, 0.5)
    @test MA.broadcast!!(+, x, y) ≈ x .+ y
    Random.seed!(1234)
    X = [BigInt(1) for _ in 1:10, _ in 1:10]
    Y = SparseArrays.sprand(Float64, 10, 10, 0.5)
    @test MA.broadcast!!(+, X, Y) ≈ X .+ Y
    return
end

function test_zero()
    v = SparseArrays.sparsevec([1, 2], [3.0, 4.0])
    MA.operate!(zero, v)
    @test isempty(SparseArrays.nonzeroinds(v))
    @test isempty(SparseArrays.nonzeros(v))
    A = SparseArrays.sprand(Float64, 2, 2, 0.5)
    MA.operate!(zero, A)
    @test iszero(SparseArrays.nnz(A))
    @test isempty(SparseArrays.nonzeros(A))
    return
end

end  # module

TestInterfaceSparseArrays.runtests()
