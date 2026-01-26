# Copyright (c) 2023 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

function neg_abs_test(f::F, x::BigFloat) where {F<:Function}
    backup = MA.copy_if_mutable(x)
    output = MA.operate!(f, x)
    @test f(backup) == output == x
    return
end

@testset "$input" for input in map(
    BigFloat,
    (-1.3, -1.0, -0.7, -0.0, 0.0, 0.3, 1.0, 1.7),
)
    @testset "-" begin
        neg_abs_test(-, input)
    end

    @testset "Base.abs" begin
        neg_abs_test(Base.abs, input)
    end
end
