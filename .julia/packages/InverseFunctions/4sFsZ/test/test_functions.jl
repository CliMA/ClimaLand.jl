# This file is a part of InverseFunctions.jl, licensed under the MIT License (MIT).

using Test
using InverseFunctions


@testset "square" begin
    for x in (0.0, 0.73)
        @test InverseFunctions.square(x) ≈ x * x
    end

    @test_throws DomainError InverseFunctions.square(-1)
end
