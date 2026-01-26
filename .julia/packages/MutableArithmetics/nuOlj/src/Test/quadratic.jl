# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

function quadratic_division_test(w, x, y, z)
    @test_rewrite w / 2
    a = 7
    b = 2
    aff = @inferred a * x + b
    @test_rewrite aff / 2
    q = @inferred b * y * z + aff
    @test_rewrite q / 2
end

function quadratic_iszero_test(w, x, y, z)
    @test !iszero(x)
    a = 7
    b = 2
    aff = @inferred a * x + b
    @test !iszero(aff)
    @test iszero(zero(typeof(aff)))
    q = @inferred b * y * z + aff
    @test !iszero(q)
    @test iszero(zero(typeof(q)))
end

function quadratic_isequal_canonical_test(w, x, y, z)
    a = 7
    b = 2
    c = 1
    @test MA.isequal_canonical((@inferred 3w + 2y), @inferred 2y + 3w)
    @test !MA.isequal_canonical((@inferred 3w + 2y + 1), @inferred 3w + 2y)
    if !MA.isequal_canonical(w, y)
        @test !MA.isequal_canonical((@inferred 3w + 2y), @inferred 3y + 2w)
    end
    if !MA.isequal_canonical(w, y)
        @test !MA.isequal_canonical((@inferred 3w + 2y), @inferred 3w + y)
    end

    aff = @inferred a * x + b
    aff2 = @inferred c * y + c
    @test !MA.isequal_canonical(aff, aff2)
    @test !MA.isequal_canonical(aff2, aff)

    @test MA.isequal_canonical(
        MA.sub_mul(MA.add_mul(aff, c, y), a, x) + c - b,
        aff2,
    )
    @test MA.isequal_canonical(
        MA.sub_mul(MA.add_mul(aff2, a, x), c, y) + b - c,
        aff,
    )

    q = @inferred b * y * z + aff
    if MA.isequal_canonical(y * z, z * y)
        @test MA.isequal_canonical(q, @inferred b * z * y + aff)
    end
    @test !MA.isequal_canonical(q, @inferred b * y * z + aff2)
    if !MA.isequal_canonical(x, y)
        @test !MA.isequal_canonical(q, @inferred b * x * z + aff)
    end
    if !MA.isequal_canonical(x, z)
        @test !MA.isequal_canonical(q, @inferred b * y * x + aff)
    end
    @test !MA.isequal_canonical(q, @inferred (b - 1) * y * z + aff)

    q2 = @inferred 8 * x * z + aff2
    if MA.isequal_canonical(x * z, z * x)
        @test MA.isequal_canonical(q2, @inferred 8z * x + aff2)
    end
    @test !MA.isequal_canonical(q2, @inferred 8x * z + aff)
    @test !MA.isequal_canonical(q2, @inferred 7x * z + aff2)
    if !MA.isequal_canonical(y, z)
        @test !MA.isequal_canonical(q2, @inferred 8x * y + aff2)
    end
    if !MA.isequal_canonical(x, y)
        @test !MA.isequal_canonical(q2, @inferred 8y * z + aff2)
    end
end

function quadratic_add_test(w, x, y, z)
    w_copy = MA.copy_if_mutable(w)
    x_copy = MA.copy_if_mutable(x)
    y_copy = MA.copy_if_mutable(y)
    z_copy = MA.copy_if_mutable(z)
    a = 7
    b = 2
    c = 1

    @testset "Variable--???" begin
        # 2-0 Variable unary
        @test (+x) === x
        @test w == w
        unary_test(x)
        add_test(x, x)
        add_test(x, w)
        add_test(w, 4)
        add_test(w, 3)
        @test_rewrite w * 4
        @test_rewrite x * y - 1
        @test_rewrite x^2
        @test_rewrite x^1
        @test_rewrite x^0
        @test_rewrite w * x
        @test_rewrite y * z - x
    end

    @test MA.isequal_canonical(w, w_copy)
    @test MA.isequal_canonical(x, x_copy)
    @test MA.isequal_canonical(y, y_copy)
    @test MA.isequal_canonical(z, z_copy)

    aff = @inferred a * x + b
    @test_rewrite a * x + b
    @test aff == aff
    aff2 = @inferred c * y + c
    @test_rewrite c * y + c

    @test MA.isequal_canonical(w, w_copy)
    @test MA.isequal_canonical(x, x_copy)
    @test MA.isequal_canonical(y, y_copy)
    @test MA.isequal_canonical(z, z_copy)

    @testset "Affine" begin
        unary_test(aff)
        add_test(aff, aff)
        add_test(aff, 1)
        add_test(aff, z)
        add_test(aff, aff2)
        @test_rewrite 7 * x - aff
        @test_rewrite aff^2
        @test_rewrite (a * x + b)^2
        @test_rewrite aff^1
        @test_rewrite (a * x + b)^1
        @test_rewrite aff^0
        @test_rewrite (a * x + b)^0
        # 3-2 AffExpr--Variable
        @test_rewrite aff * z
        @test_rewrite aff - a * x
        # 3-3 AffExpr--AffExpr
        @test_rewrite aff * aff2
    end

    @test MA.isequal_canonical(w, w_copy)
    @test MA.isequal_canonical(x, x_copy)
    @test MA.isequal_canonical(y, y_copy)
    @test MA.isequal_canonical(z, z_copy)

    @testset "Quadratic" begin
        @test_rewrite 2 * x * x + 1 * y * y + z + 3

        @test_rewrite b * y * z + aff
        q = @inferred b * y * z + aff
        @test q == q
        unary_test(q)
        add_test(q, q)
        add_test(q, 1)
        add_test(q, w)
        add_test(q, aff2)

        @test_rewrite 8 * x * z + aff2
        q2 = @inferred 8 * x * z + aff2
        add_test(q, q2)

        vx = [x]
        vy = [z]
        v1 = [1]
        @test_rewrite (1 + vx'vy) * 1
        @test_rewrite (1 + vx'v1) * x
        @test_rewrite (1 + v1'vy) * y
    end

    @test MA.isequal_canonical(w, w_copy)
    @test MA.isequal_canonical(x, x_copy)
    @test MA.isequal_canonical(y, y_copy)
    @test MA.isequal_canonical(z, z_copy)
end

const quadratic_tests = Dict(
    "quadratic_division" => quadratic_division_test,
    "quadratic_iszero" => quadratic_iszero_test,
    "quadratic_isequal_canonical" => quadratic_isequal_canonical_test,
    "quadratic_add_canonical" => quadratic_add_test,
)

@test_suite quadratic
