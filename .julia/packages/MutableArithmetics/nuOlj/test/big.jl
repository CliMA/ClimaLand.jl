# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

function allocation_test(
    op,
    T,
    short,
    short_to,
    n;
    a = T(2),
    b = T(3),
    c = T(4),
)
    @test MA.promote_operation(op, T, T) == T
    alloc_test(() -> MA.promote_operation(op, T, T), 0)
    if op != div && op != -
        @test MA.promote_operation(op, T, T, T) == T
        alloc_test(() -> MA.promote_operation(op, T, T, T), 0)
    end
    g = op(a, b)
    @test c === short_to(c, a, b)
    @test g == c
    A = MA.copy_if_mutable(a)
    @test A === short(A, b)
    @test g == A
    alloc_test_le(() -> short(A, b), n)
    alloc_test_le(() -> short_to(c, a, b), n)
    @test g == MA.buffered_operate!(nothing, op, MA.copy_if_mutable(a), b)
    @test g == MA.buffered_operate_to!(nothing, c, op, a, b)
    buffer = MA.buffer_for(op, typeof(a), typeof(b))
    @test g == MA.buffered_operate_to!(buffer, c, op, a, b)
    alloc_test(() -> MA.buffered_operate_to!(buffer, c, op, a, b), 0)
    return
end

function add_sub_mul_test(op, T; a = T(2), b = T(3), c = T(4))
    g = op(a, b, c)
    @test g == MA.buffered_operate!(nothing, op, MA.copy_if_mutable(a), b, c)
    buffer = MA.buffer_for(op, typeof(a), typeof(b), typeof(c))
    return alloc_test(() -> MA.buffered_operate!(buffer, op, a, b, c), 0)
end
@testset "$T" for T in [BigInt, BigFloat, Rational{BigInt}]
    MA.Test.int_test(T)
    @testset "Allocation" begin
        MAX_ALLOC = T <: Rational ? 280 : 0
        allocation_test(+, T, MA.add!!, MA.add_to!!, MAX_ALLOC)
        allocation_test(-, T, MA.sub!!, MA.sub_to!!, MAX_ALLOC)
        allocation_test(*, T, MA.mul!!, MA.mul_to!!, MAX_ALLOC)
        add_sub_mul_test(MA.add_mul, T)
        add_sub_mul_test(MA.sub_mul, T)
        if T <: Rational # https://github.com/jump-dev/MutableArithmetics.jl/issues/167
            allocation_test(
                +,
                T,
                MA.add!!,
                MA.add_to!!,
                MAX_ALLOC,
                a = T(1 // 2),
                b = T(3 // 2),
                c = T(5 // 2),
            )
            allocation_test(
                -,
                T,
                MA.sub!!,
                MA.sub_to!!,
                MAX_ALLOC,
                a = T(1 // 2),
                b = T(3 // 2),
                c = T(5 // 2),
            )
        end
        # Requires https://github.com/JuliaLang/julia/commit/3f92832df042198b2daefc1f7ca609db38cb8173
        # for `gcd` to be defined on `Rational`.
        if T == BigInt
            allocation_test(div, T, MA.div!!, MA.div_to!!, 0)
        end
        if T == BigInt || T == Rational{BigInt}
            allocation_test(gcd, T, MA.gcd!!, MA.gcd_to!!, 0)
            allocation_test(lcm, T, MA.lcm!!, MA.lcm_to!!, 0)
        end
    end
end
