# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

using Test
import MutableArithmetics as MA

@testset "promote_operation" begin
    @test MA.promote_operation(MA.zero, Int) == Int
    @test MA.promote_operation(MA.one, Int) == Int
    @test MA.promote_operation(+, Int, Int) == Int
    @test MA.promote_operation(+, Int, Int, Int) == Int
    @test MA.promote_operation(-, Int, Int) == Int
    @test MA.promote_operation(*, Int, Int) == Int
    @test MA.promote_operation(*, Int, Int, Int) == Int
    @test MA.promote_operation(gcd, Int, Int) == Int
    @test MA.promote_operation(gcd, Int, Int, Int) == Int
    @test MA.promote_operation(MA.add_mul, Int, Int, Int) == Int
    for op in [+, -, MA.add_mul, MA.sub_mul]
        err = ErrorException(
            "Operation `$op` between `$(Vector{Int})` and `$Int` is not allowed. You should use broadcast.",
        )
        @test_throws err MA.promote_operation(op, Vector{Int}, Int)
        err = ErrorException(
            "Operation `$op` between `$Int` and `$(Vector{Int})` is not allowed. You should use broadcast.",
        )
        @test_throws err MA.promote_operation(op, Int, Vector{Int})
    end
    for op in [+, -, *, /, div]
        @test MA.promote_operation(op, Int, Number) == Number
        @test MA.promote_operation(op, Number, Int) == Number
    end
    @test MA.promote_operation(/, Int, Integer) == Float64
    @test MA.promote_operation(/, Integer, Integer) == Float64
    @test MA.promote_operation(/, Integer, Int) == Float64
    @test MA.promote_operation(gcd, Int, Integer) == Integer
    @test MA.promote_operation(gcd, Integer, Integer) == Integer
    @test MA.promote_operation(gcd, Integer, Int) == Integer
    @test MA.promote_operation(&, Integer, Integer, Integer) == Integer
    @test MA.promote_operation(&, Integer, Integer, Int) == Integer
end
@testset "add_to!! / add!!" begin
    @test MA.mutability(Int, MA.add_to!!, Int, Int) isa MA.IsNotMutable
    @test MA.mutability(Int, MA.add!!, Int) isa MA.IsNotMutable
    a = 5
    b = 28
    c = 41
    @test MA.add_to!!(a, b, c) == 69
    @test MA.add!!(b, c) == 69
    a = 165
    b = 255
    @test MA.add!!(a, b) == 420
end

@testset "mul_to!! / mul!!" begin
    @test MA.mutability(Int, MA.mul_to!!, Int, Int) isa MA.IsNotMutable
    @test MA.mutability(Int, MA.mul!!, Int) isa MA.IsNotMutable
    a = 5
    b = 23
    c = 3
    @test MA.mul_to!!(a, b, c) == 69
    @test MA.mul!!(b, c) == 69
    a = 15
    b = 28
    @test MA.mul!!(a, b) == 420
end

@testset "add_mul_to!! / add_mul!! / add_mul_buf_to!! /add_mul_buf!!" begin
    @test MA.mutability(Int, MA.add_mul_to!!, Int, Int, Int) isa MA.IsNotMutable
    @test MA.mutability(Int, MA.add_mul!!, Int, Int) isa MA.IsNotMutable
    @test MA.mutability(Int, MA.add_mul_buf_to!!, Int, Int, Int, Int) isa
          MA.IsNotMutable
    @test MA.mutability(Int, MA.add_mul_buf!!, Int, Int, Int) isa
          MA.IsNotMutable
    a = 5
    b = 9
    c = 3
    d = 20
    buf = 24
    @test MA.add_mul_to!!(a, b, c, d) == 69
    @test MA.add_mul!!(b, c, d) == 69
    @test MA.add_mul_buf_to!!(buf, a, b, c, d) == 69
    @test MA.add_mul_buf!!(buf, b, c, d) == 69
    a = 148
    b = 16
    c = 17
    @test MA.add_mul!!(a, b, c) == 420
    a = 148
    b = 16
    c = 17
    buf = 56
    @test MA.add_mul_buf!!(buf, a, b, c) == 420
    a = 148
    b = 16
    c = 17
    d = 42
    buf = 56
    @test MA.add_mul_buf_to!!(buf, d, a, b, c) == 420
end

@testset "zero!!" begin
    @test MA.mutability(Int, MA.zero!!) isa MA.IsNotMutable
    a = 5
    @test MA.zero!!(a) == 0
end

@testset "one!!" begin
    @test MA.mutability(Int, MA.one!!) isa MA.IsNotMutable
    a = 5
    @test MA.one!!(a) == 1
end

@testset "Zero / Int" begin
    @test MA.Zero() / 1 == MA.Zero()
    @test_throws DivideError MA.Zero() / 0
end

@testset "Division" begin
    @test 1 / 2 == MA.operate!!(/, 1, 2)
end
