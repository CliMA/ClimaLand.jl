# Copyright (c) 2023 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

function test_fma_output_values(x::F, y::F, z::F) where {F<:BigFloat}
    two_roundings_reference = x * y + z
    one_rounding_reference = fma(x, y, z)
    @test one_rounding_reference != two_roundings_reference

    @testset "fma $op output values" for op in (MA.operate!, MA.operate!!)
        (a, b, c) = map(MA.copy_if_mutable, (x, y, z))
        @inferred op(fma, a, b, c)
        @test one_rounding_reference == a
        @test y == b
        @test z == c
    end

    @testset "fma $op output values" for op in (MA.operate_to!, MA.operate_to!!)
        (a, b, c) = map(MA.copy_if_mutable, (x, y, z))
        out = F()
        @inferred op(out, fma, a, b, c)
        @test one_rounding_reference == out
        @test x == a
        @test y == b
        @test z == c
    end

    return nothing
end

function test_fma_output_values(x::F, y::F, z::F) where {F<:Float64}
    return test_fma_output_values(map(BigFloat, (x, y, z))...)
end

function test_fma_output_values_func(x::F, y::F, z::F) where {F<:Float64}
    return let x = x, y = y, z = z
        () -> test_fma_output_values(x, y, z)
    end
end

@testset "fma output values: $exp_x $exp_y $sign_x $sign_y" for exp_x in (-3):3,
    exp_y in (-3):3,
    sign_x in (-1, 1),
    sign_y in (-1, 1)

    # Assuming a two-bit mantissa
    significand_length = 2

    sign_z = -sign_x * sign_y
    exp_z = exp_x + exp_y

    base = 2.0
    bit_0 = 2.0^0
    bit_1 = 2.0^-1

    x = sign_x * base^exp_x * (bit_0 + bit_1)
    y = sign_y * base^exp_y * (bit_0 + bit_1)
    z = sign_z * base^exp_z * (bit_0 + bit_1)

    setprecision(
        test_fma_output_values_func(x, y, z),
        BigFloat,
        significand_length,
    )
end

@testset "muladd operate_to!! type inferred" begin
    m1 = BigFloat(-1.0)
    out = BigFloat()
    @test iszero(@inferred MA.operate_to!!(out, Base.muladd, m1, m1, m1))
end

@testset "muladd operate!! type inferred" begin
    x = BigFloat(-1.0)
    y = BigFloat(-1.0)
    z = BigFloat(-1.0)
    @test iszero(@inferred MA.operate!!(Base.muladd, x, y, z))
end

@testset "fma $op doesn't allocate" for op in (MA.operate_to!, MA.operate_to!!)
    alloc_test(let op = op, o = big"1.3", x = big"1.3"
        () -> op(o, Base.fma, x, x, x)
    end, 0)
end

@testset "fma $op doesn't allocate" for op in (MA.operate!, MA.operate!!)
    alloc_test(let op = op, x = big"1.3", y = big"1.3", z = big"1.3"
        () -> op(Base.fma, x, y, z)
    end, 0)
end
