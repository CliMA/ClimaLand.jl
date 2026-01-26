# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

function int_add_test(::Type{T}) where {T}
    @testset "add_to! / add!" begin
        @test MA.mutability(T, +, T, T) isa MA.IsMutable

        t(n) = convert(T, n)
        a = t(5)
        b = t(28)
        c = t(41)
        expected = t(69)
        @test MA.isequal_canonical(MA.add_to!!(a, b, c), expected)
        @test MA.isequal_canonical(a, expected)
        @test MA.isequal_canonical(MA.add!!(b, c), expected)
        @test MA.isequal_canonical(b, expected)

        a = t(165)
        b = t(255)
        expected = t(420)
        @test MA.isequal_canonical(MA.add!!(a, b), expected)
        @test MA.isequal_canonical(a, expected)
    end
end
function int_sub_test(::Type{T}) where {T}
    @testset "sub_to! / sub!" begin
        @test MA.mutability(T, -, T, T) isa MA.IsMutable

        t(n) = convert(T, n)
        a = t(5)
        b = t(28)
        c = t(41)
        expected = t(-13)
        @test MA.isequal_canonical(MA.sub_to!!(a, b, c), expected)
        @test MA.isequal_canonical(a, expected)
        @test MA.isequal_canonical(MA.sub!!(b, c), expected)
        @test MA.isequal_canonical(b, expected)

        a = t(165)
        b = t(255)
        expected = t(-90)
        @test MA.isequal_canonical(MA.sub!!(a, b), expected)
        @test MA.isequal_canonical(a, expected)
    end
end
function int_mul_test(::Type{T}) where {T}
    @testset "mul_to! / mul!" begin
        @test MA.mutability(T, *, T, T) isa MA.IsMutable

        t(n) = convert(T, n)
        a = t(5)
        b = t(23)
        c = t(3)
        @test MA.isequal_canonical(MA.mul_to!!(a, b, c), t(69))
        @test MA.isequal_canonical(a, t(69))
        @test MA.isequal_canonical(MA.mul!!(b, c), t(69))
        @test MA.isequal_canonical(b, t(69))

        a = t(15)
        b = t(28)
        @test MA.isequal_canonical(MA.mul!!(a, b), t(420))
        @test MA.isequal_canonical(a, t(420))
    end
end
function int_add_mul_test(::Type{T}) where {T}
    @testset "add_mul_to! / add_mul! / add_mul_buf_to! / add_mul_buf!" begin
        @test MA.mutability(T, MA.add_mul, T, T) isa MA.IsMutable
        @test MA.mutability(T, MA.add_mul, T, T, T) isa MA.IsMutable
        @test MA.mutability(T, MA.add_mul, T, T, T, T) isa MA.IsMutable

        t(n) = convert(T, n)
        a = t(5)
        b = t(9)
        c = t(3)
        d = t(20)
        buf = MA.buffer_for(MA.add_mul, typeof(a), typeof(b), typeof(c))

        @test MA.isequal_canonical(MA.add_mul_to!!(a, b, c, d), t(69))
        @test MA.isequal_canonical(a, t(69))
        a = t(5)
        @test MA.isequal_canonical(MA.add_mul!!(b, c, d), t(69))
        @test MA.isequal_canonical(b, t(69))
        b = t(9)

        @test MA.isequal_canonical(MA.add_mul_buf_to!!(buf, a, b, c, d), t(69))
        @test MA.isequal_canonical(a, t(69))
        @test MA.isequal_canonical(MA.add_mul_buf!!(buf, b, c, d), t(69))
        @test MA.isequal_canonical(b, t(69))

        a = t(148)
        b = t(16)
        c = t(17)
        d = t(42)
        @test MA.isequal_canonical(MA.add_mul!!(a, b, c), t(420))
        @test MA.isequal_canonical(a, t(420))
        a = t(148)
        @test MA.isequal_canonical(MA.add_mul_buf_to!!(buf, d, a, b, c), t(420))
        @test MA.isequal_canonical(d, t(420))
        @test MA.isequal_canonical(MA.add_mul_buf!!(buf, a, b, c), t(420))
        @test MA.isequal_canonical(a, t(420))
    end
end
function int_sub_mul_test(::Type{T}) where {T}
    @testset "sub_mul_to! / sub_mul! / sub_mul_buf_to! / sub_mul_buf!" begin
        @test MA.mutability(T, MA.sub_mul, T, T) isa MA.IsMutable
        @test MA.mutability(T, MA.sub_mul, T, T, T) isa MA.IsMutable
        @test MA.mutability(T, MA.sub_mul, T, T, T, T) isa MA.IsMutable

        t(n) = convert(T, n)
        a = t(5)
        b = t(9)
        c = t(3)
        d = t(20)
        buf = t(24)

        expected = t(-51)
        @test MA.isequal_canonical(MA.sub_mul_to!!(a, b, c, d), expected)
        @test MA.isequal_canonical(a, expected)
        a = t(5)
        @test MA.isequal_canonical(MA.sub_mul!!(b, c, d), expected)
        @test MA.isequal_canonical(b, expected)
        b = t(9)

        @test MA.isequal_canonical(
            MA.sub_mul_buf_to!!(buf, a, b, c, d),
            expected,
        )
        @test MA.isequal_canonical(a, expected)
        @test MA.isequal_canonical(MA.sub_mul_buf!!(buf, b, c, d), expected)
        @test MA.isequal_canonical(b, expected)

        a = t(148)
        b = t(16)
        c = t(17)
        d = t(42)
        buf = t(56)
        expected = t(-124)
        @test MA.isequal_canonical(MA.sub_mul!!(a, b, c), expected)
        @test MA.isequal_canonical(a, expected)
        a = t(148)
        @test MA.isequal_canonical(
            MA.sub_mul_buf_to!!(buf, d, a, b, c),
            expected,
        )
        @test MA.isequal_canonical(d, expected)
        @test MA.isequal_canonical(MA.sub_mul_buf!!(buf, a, b, c), expected)
        @test MA.isequal_canonical(a, expected)
    end
end

function int_zero_test(::Type{T}) where {T}
    @testset "zero!" begin
        @test MA.mutability(T, zero, T) isa MA.IsMutable

        t(n) = convert(T, n)
        a = t(5)
        @test MA.isequal_canonical(MA.zero!!(a), t(0))
        @test MA.isequal_canonical(a, t(0))
        @test iszero(a)
    end
end

function int_one_test(::Type{T}) where {T}
    @testset "one!" begin
        @test MA.mutability(T, one, T) isa MA.IsMutable

        t(n) = convert(T, n)
        a = t(5)
        @test MA.isequal_canonical(MA.one!!(a), t(1))
        @test MA.isequal_canonical(a, t(1))
        @test isone(a)
    end
end

const int_tests = Dict(
    "int_add" => int_add_test,
    "int_mul" => int_mul_test,
    "int_add_mul" => int_add_mul_test,
    "int_zero" => int_zero_test,
    "int_one" => int_one_test,
)

@test_suite int
