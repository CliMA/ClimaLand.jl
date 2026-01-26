# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

module TestRewriteGeneric

using Test

import MutableArithmetics as MA

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

macro test_rewrite(expr)
    return quote
        esc(
            @test MA.isequal_canonical(
                MA.@rewrite($expr, move_factors_into_sums = false),
                $expr,
            )
        )
    end
end

function test_rewrite()
    x, expr = MA.rewrite(1 + 1, move_factors_into_sums = false)
    @test x isa Symbol
    @test Meta.isexpr(expr, :(=), 2)
    return
end

function test_rewrite_basic_ops()
    @test_rewrite 1 + 2
    @test_rewrite 1 - 2
    @test_rewrite 1 * 2
    @test_rewrite 1 / 2
    @test_rewrite 1 + 2 - 3
    @test_rewrite 1 - 2 + 4
    @test_rewrite 1 * 2 / 3
    @test_rewrite 1 / 2 * 3
    return
end

function test_rewrite_not_call()
    x = [1, 2, 3]
    for i in 1:3
        @test MA.@rewrite(x[i], move_factors_into_sums = false) == i
    end
    return
end

function test_rewrite_zero_argument()
    @test_rewrite Dict()
    return
end

function test_rewrite_sum_to_add_mul()
    @test_rewrite +1
    @test_rewrite +(-2)
    x = [1.2]
    @test MA.@rewrite(+x, move_factors_into_sums = false) == x
    @test MA.@rewrite(x + x, move_factors_into_sums = false) == 2 * x
    @test MA.@rewrite(+(x, x, x), move_factors_into_sums = false) == 3 * x
    @test MA.@rewrite(+(x, x, x, x), move_factors_into_sums = false) == 4 * x
    return
end

function test_rewrite_prod_to_add_mul()
    @test_rewrite 1 * 1
    @test_rewrite 2 * -2
    A = [1.0 2.0; 3.0 4.0]
    x = [5.0, 6.0]
    @test MA.@rewrite(A * x, move_factors_into_sums = false) == A * x
    return
end

function test_rewrite_nonconcrete_vector()
    x = [5.0, 6.0]
    y = Vector{Union{Float64,String}}(x)
    @test MA.@rewrite(x' * y, move_factors_into_sums = false) == x' * y
    @test MA.@rewrite(x .+ y, move_factors_into_sums = false) == x .+ y
    @test MA.@rewrite(x + y, move_factors_into_sums = false) == x + x
    return
end

function test_rewrite_minus_to_add_mul()
    @test_rewrite -1
    @test_rewrite -(+2)
    @test_rewrite -(-2)
    x = [1.2]
    @test MA.@rewrite(-x, move_factors_into_sums = false) == -x
    @test MA.@rewrite(x - x, move_factors_into_sums = false) == [0.0]
    @test MA.@rewrite(-(x, x, x), move_factors_into_sums = false) == -1 * x
    @test MA.@rewrite(-(x, x, x, x), move_factors_into_sums = false) == -2 * x
    return
end

function test_rewrite_sum()
    @test MA.@rewrite(sum(i for i in 1:0), move_factors_into_sums = false) ==
          MA.Zero()
    @test_rewrite sum(i for i in 1:10)
    @test_rewrite sum(i + i^2 for i in 1:10)
    @test_rewrite sum(i * i for i in 1:10)
    return
end

function test_rewrite_unknown_generators()
    @test_rewrite prod(i for i in 1:10)
    @test_rewrite minimum(i for i in 1:10)
    @test_rewrite maximum(i for i in 1:10)
    return
end

function test_rewrite_generator()
    # Unfiltered generators
    # Univariate
    @test_rewrite sum(i for i in 1:2)
    @test_rewrite sum(j^2 for j in 2:3)
    # Binary
    @test_rewrite sum(i + j^2 for i in 1:2, j in 2:3)
    @test_rewrite sum(i + j^2 for i in 1:2 for j in 2:3)
    # Tertiary
    @test_rewrite sum(i + j^2 for i in 1:2, j in 2:3, k in 3:4)
    @test_rewrite sum(i + j^2 for i in 1:2, j in 2:3 for k in i:j)
    @test_rewrite sum(i + j^2 for i in 1:2 for j in 2:3, k in 3:4)
    @test_rewrite sum(i + j^2 for i in 1:2 for j in 2:3 for k in 3:4)
    # Generators with dependent variables
    # This syntax is unsupported by Julia!
    @test MA.@rewrite(
        sum(i + j^2 for i in 1:2, j in i:3),
        move_factors_into_sums = false,
    ) == 34
    @test MA.@rewrite(
        sum(i + j^2 for i in 1:2, j in 2:3, k in i:j),
        move_factors_into_sums = false,
    ) == 68
    # Unnivariate generators with an if statement
    @test_rewrite sum(i for i in 1:2 if i >= 1)
    @test_rewrite sum(j^2 for j in 2:3 if j <= 3)
    @test_rewrite sum(i for i in 1:2 if isodd(i))
    @test_rewrite sum(j^2 for j in 2:3 if isodd(j))
    # Binary generators with an if statement
    @test_rewrite sum(i + j^2 for i in 1:2, j in 2:3 if i + j >= 0)
    @test_rewrite sum(i + j^2 for i in 1:2, j in 2:3 if i <= 1)
    @test_rewrite sum(i + j^2 for i in 1:2, j in 2:3 if j <= 2)
    # Flatten with generator on first key
    @test_rewrite sum(i for i in 1:2 if isodd(i) for j in 1:3)
    @test_rewrite sum(i for i in 1:3 if i <= 2 for j in 1:2)
    # Flatten with generator on second key
    @test_rewrite sum(i for i in 1:2 for j in 1:i if isodd(j))
    @test_rewrite sum(i for i in 1:2 for j in 1:i if i != j)
    # Generator with two filters
    @test_rewrite sum(i for i in 1:3 if i <= 2 for j in 1:2 if j >= 1)
    return
end

function test_rewrite_linear_algebra()
    @test_rewrite [1 2; 3 4] * [5.0, 6.0]
    @test_rewrite [1 2; 3 4] * [5, 6]
    @test_rewrite [1 2; 3 4] * [1 2; 3 4]
    @test_rewrite [1 2; 3 4] * Float64[1 2; 3 4]
    @test_rewrite [5, 6]' * [1 2; 3 4] * [5, 6]
    @test_rewrite 2 * [1 2; 3 4]
    @test_rewrite [1 2; 3 4] * 2
    @test_rewrite [1, 2] * 2.0
    @test_rewrite 2.0 * [1, 2]
    x = [1, 2]
    A = [1 2; 3 4]
    y = reshape(x, (1, length(x))) * A * x .- 1
    @test MA.@rewrite(
        reshape(x, (1, length(x))) * A * x .- 1,
        move_factors_into_sums = false,
    ) == y
    return
end

function test_rewrite_broadcast_add_mul()
    @test_rewrite 1 .+ 2
    @test_rewrite [1, 2] .+ [2, 3]
    @test_rewrite [5.0, 6.0] .+ [1 2; 3 4] * [5.0, 6.0]
    @test_rewrite .+(BigInt(1))

    return
end

function test_rewrite_broadcast_sub_mul()
    @test_rewrite 1 .* 2
    @test_rewrite [1, 2] .* [2, 3]
    @test_rewrite [5.0, 6.0] .- [1 2; 3 4] * [5.0, 6.0]
    @test_rewrite .-(BigInt(1))
    return
end

function test_rewrite_repeated_sums()
    n = 10
    x = BigInt[2 for _ in 1:n]
    @test MA.@rewrite(sum(x[i] for i in 1:n) + sum(x[i] for i in 1:n)) == 4 * n
    @test MA.@rewrite(sum(x[i] for i in 1:n) - sum(x[i] for i in 1:n)) == 0
    @test MA.@rewrite(sum(x[i]^2 for i in 1:n) - sum(x[i] for i in 1:n)) ==
          sum(x .^ 2) - sum(x)
    @test MA.@rewrite(sum(x[i]^2 for i in 1:n) + sum(x[i] for i in 1:n)) ==
          sum(x .^ 2) + sum(x)
    @test MA.@rewrite(sum(x[i] for i in 1:n) - sum(x[i]^2 for i in 1:n)) ==
          sum(x) - sum(x .^ 2)
    @test MA.@rewrite(sum(x[i] for i in 1:n) + sum(x[i]^2 for i in 1:n)) ==
          sum(x) + sum(x .^ 2)
    return
end

function test_multi_prod()
    @test_rewrite *(sum([BigInt(5)] for _ in 1:2), 2)
    @test_rewrite *(sum([BigInt(5)] for _ in 1:2), BigInt(2))
    @test_rewrite *(sum([BigInt(5)] for _ in 1:2)', [2])
    @test_rewrite *(sum(BigInt(5) for _ in 1:2), BigInt(2), 3, BigInt(4))
    return
end

function test_rewrite_sums_init()
    @test_rewrite sum(i for i in 1:2; init = 2)
    @test_rewrite sum(i for i in 1:2; init = -1)
    @test_rewrite sum(i * j for i in 1:2 for j in 1:3; init = -1)
    @test_rewrite sum(i * j for i in 1:2 for j in 1:3; init = 2 * 3)
    @test_rewrite sum(i * j for i in 1:2 for j in i:3; init = 2 * 3)
    return
end

function test_rewrite_sums_init_positional()
    # Needed so we don't replace `, init` with `; init`.
    #! format:off
    @test MA.@rewrite(
        sum(i for i in 1:2, init = 2),
        move_factors_into_sums = false,
    ) == 5
    @test MA.@rewrite(
        sum(i for i in 1:2, init = -1),
        move_factors_into_sums = false,
    ) == 2
    @test MA.@rewrite(
        sum(i * j for i in 1:2 for j in 1:3, init = -1),
        move_factors_into_sums = false,
    ) == 17
    @test MA.@rewrite(
        sum(i * j for i in 1:2 for j in 1:3, init = 2 * 3),
        move_factors_into_sums = false,
    ) == 24
    @test MA.@rewrite(
        sum(i * j for i in 1:2 for j in i:3, init = 2 * 3),
        move_factors_into_sums = false,
    ) == 22
    #! format:on
    return
end

function test_rewrite_kwarg()
    f(; x) = x^2
    @test MA.@rewrite(f(; x = 2), move_factors_into_sums = false) == 4
    g(x; y) = x + y^2
    @test MA.@rewrite(g(2; y = 2), move_factors_into_sums = false) == 6
    return
end

function test_rewrite_sum_function()
    @test_rewrite sum(i -> i^2, 1:2)
    @test_rewrite sum(i -> i^2, 1:2; init = 3)
    f(x) = x^2
    y = 1:2
    @test MA.@rewrite(sum(f, y), move_factors_into_sums = false) == sum(f, y)
    return
end

struct NoProduct <: MA.AbstractMutable
    x::Int
end

MA.operate(::typeof(*), x::NoProduct, y::NoProduct) = NoProduct(x.x * y.x)

function test_no_product()
    x, y = NoProduct(2), NoProduct(3)
    @test MA.@rewrite(x * y, move_factors_into_sums = false) == NoProduct(6)
    return
end

function test_splatting()
    x = [1, 2, 3]
    @test MA.@rewrite(+(x...), move_factors_into_sums = false) == 6
    @test MA.@rewrite(+(4, x...), move_factors_into_sums = false) == 10
    @test MA.@rewrite(+(x..., 4), move_factors_into_sums = false) == 10
    @test MA.@rewrite(+(4, x..., 5), move_factors_into_sums = false) == 15
    @test MA.@rewrite(*(x...), move_factors_into_sums = false) == 6
    @test MA.@rewrite(*(4, x...), move_factors_into_sums = false) == 24
    @test MA.@rewrite(*(x..., 4), move_factors_into_sums = false) == 24
    @test MA.@rewrite(*(4, x..., 5), move_factors_into_sums = false) == 120
    @test MA.@rewrite(
        +(4, x..., *(4, x..., 5)),
        move_factors_into_sums = false,
    ) == 130
    @test MA.@rewrite(vcat(x...), move_factors_into_sums = false) == x
    return
end

struct _KwargRef{K,V}
    data::Dict{K,V}
end

Base.getindex(x::_KwargRef; i) = x.data[i]

function test_rewrite_kw_in_ref()
    x = _KwargRef(Dict(i => i + 1 for i in 2:4))
    @test MA.@rewrite(
        sum(x[i = j] for j in 2:4),
        move_factors_into_sums = false,
    ) == 12
    return
end

function test_rewrite_expression()
    x = [1.2]
    @test MA.@rewrite(x + 2 * x, move_factors_into_sums = false) == 3x
    @test MA.@rewrite(x + *(2, x, 3), move_factors_into_sums = false) ==
          x + *(2, x, 3)
    y = 1.2
    @test MA.@rewrite(
        sum(y for i in 1:2) + 2y,
        move_factors_into_sums = false
    ) == 4y
    @test MA.@rewrite(
        sum(y for i in 1:2) + y * y,
        move_factors_into_sums = false
    ) == 2y + y^2
    @test MA.@rewrite(y + 2 * y, move_factors_into_sums = false) == 3y
    @test MA.@rewrite(y + *(2, y, 3), move_factors_into_sums = false) ==
          y + *(2, y, 3)
    return
end

function test_rewrite_generic_sum_dims()
    @test_rewrite sum([1 2; 3 4]; dims = 1)
    @test_rewrite sum([1 2; 3 4]; dims = 2)
    @test_rewrite sum([1 2; 3 4]; dims = 1, init = 0)
    @test_rewrite sum([1 2; 3 4]; dims = 2, init = 0)
    @test_rewrite sum([1 2; 3 4]; init = 0, dims = 1)
    @test_rewrite sum([1 2; 3 4]; init = 0, dims = 2)
    @test_rewrite sum([1 2; 3 4], dims = 1)
    @test_rewrite sum([1 2; 3 4], dims = 2)
    @test_rewrite sum([1 2; 3 4], dims = 1, init = 0)
    @test_rewrite sum([1 2; 3 4], dims = 2, init = 0)
    @test_rewrite sum([1 2; 3 4], init = 0, dims = 1)
    @test_rewrite sum([1 2; 3 4], init = 0, dims = 2)
    return
end

function test_rewrite_block()
    @test_rewrite begin
        x = 1
        y = x + 2
        z = 3 * y
    end
    @test_rewrite begin
        x = [1]
        y = x + [2]
        z = 3 * y
    end
    return
end

function test_rewrite_ifelse()
    @test_rewrite begin
        x = -1
        y = [3.0]
        if x < 1
            y .+ x
        else
            2 * x
        end
    end
    @test_rewrite begin
        x = 2
        y = [3.0]
        if x < 1
            y .+ x
        else
            2 * x
        end
    end
    @test_rewrite begin
        x = 2
        y = [3.0, 4.0]
        if x < 1
            y .+ x
        elseif length(y) == 2
            0.0
        else
            2 * x
        end
    end
    @test_rewrite begin
        x = 2
        y = Float64[]
        if x < 1
            y .+ x
        elseif length(y) == 2
            0.0
        elseif isempty(y)
            -1.0
        else
            2 * x
        end
    end
    @test_rewrite begin
        x = 2
        y = Float64[1.0]
        if x < 1
            1.0
        elseif length(y) == 2
            2.0
        elseif isempty(y)
            3.0
        else
            4.0
        end
    end
    return
end

function test_return_is_mutable()
    function _rewrite(expr)
        return MA.rewrite(
            expr;
            move_factors_into_sums = false,
            return_is_mutable = true,
        )
    end
    x, expr, is_mutable = _rewrite(1)
    @test x isa Symbol
    @test Meta.isexpr(expr, :(=), 2)
    @test is_mutable
    y = 1
    x, expr, is_mutable = _rewrite(:(y))
    @test x isa Symbol
    @test Meta.isexpr(expr, :(=), 2)
    @test !is_mutable
    x, expr, is_mutable = _rewrite(:(1 + 1))
    @test x isa Symbol
    @test Meta.isexpr(expr, :(=), 2)
    @test is_mutable
    @test_throws(
        AssertionError,
        MA.rewrite(
            :(1 + 1);
            move_factors_into_sums = true,
            return_is_mutable = true,
        ),
    )
    return
end

function test_rewrite_sum_unary_minus()
    x = big.(1:3)
    y = MA.@rewrite(sum(-x[i] for i in 1:3), move_factors_into_sums = false)
    @test y == big(-6)
    return
end

function test_allocations_rewrite_unary_minus()
    N = 100
    x = big.(1:N)
    a = big(-1)
    # precompilaton
    MA.@rewrite(sum(-x[i] for i in 1:N), move_factors_into_sums = false)
    MA.@rewrite(-sum(x[i] for i in 1:N), move_factors_into_sums = false)
    MA.@rewrite(sum(a * x[i] for i in 1:N), move_factors_into_sums = false)
    sum(-x[i] for i in 1:N)
    total = @allocated sum(-x[i] for i in 1:N)
    # actual
    value = @allocated(
        MA.@rewrite(sum(-x[i] for i in 1:N), move_factors_into_sums = false),
    )
    @test value < total
    value = @allocated(
        MA.@rewrite(-sum(x[i] for i in 1:N), move_factors_into_sums = false),
    )
    @test value < total
    value = @allocated(
        MA.@rewrite(sum(a * x[i] for i in 1:N), move_factors_into_sums = false),
    )
    @test value < total
    return
end

function test_allocations_rewrite_mult()
    x = BigFloat[1, 2, 3]
    y = MA.@rewrite(x[1] * x[2] * x[3], move_factors_into_sums = false)
    @test y == BigFloat(6)
    @test x == BigFloat[1, 2, 3]
    total = @allocated (x[1] * x[2]) * x[3]
    @test @allocated(
        MA.@rewrite(x[1] * x[2] * x[3], move_factors_into_sums = false),
    ) < total
    return
end

function test_rewrite_init_symbol()
    x = Int[]
    y = MA.@rewrite(sum(x; init = 0), move_factors_into_sums = false)
    @test y == 0
    y = MA.@rewrite(sum(x, init = 0), move_factors_into_sums = false)
    @test y == 0
    return
end

end  # module

TestRewriteGeneric.runtests()
