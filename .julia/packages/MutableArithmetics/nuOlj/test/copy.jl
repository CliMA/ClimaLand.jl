# Copyright (c) 2023 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

@testset "copy: $T" for T in (
    Float64,
    BigFloat,
    Int,
    BigInt,
    Rational{Int},
    Rational{BigInt},
)
    @test MA.operate!!(copy, T(2)) == 2
    @test MA.operate_to!!(T(3), copy, T(2)) == 2
    if MA.mutability(T, copy, T) == MA.IsMutable()
        @testset "mutable" begin
            @testset "correctness" begin
                x = T(2)
                y = T(3)
                @test MA.operate!(copy, x) === x == 2
                @test MA.operate_to!(y, copy, x) === y == 2
            end
            @testset "alloc" begin
                f = let x = T(2)
                    () -> MA.operate!(copy, x)
                end
                g = let x = T(2), y = T(3)
                    () -> MA.operate_to!(y, copy, x)
                end
                alloc_test(f, 0)
                alloc_test(g, 0)
            end
        end
    end
end
