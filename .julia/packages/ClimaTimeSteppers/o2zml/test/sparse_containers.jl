using ClimaTimeSteppers: SparseContainer

using Test
@testset "SparseContainer" begin
    a1 = ones(3) .* 1
    a2 = ones(3) .* 2
    a3 = ones(3) .* 3
    a4 = ones(3) .* 4
    v = SparseContainer((a1, a2, a3, a4), (1, 3, 5, 7))
    @test v[1] == ones(3) .* 1
    @test v[3] == ones(3) .* 2
    @test v[5] == ones(3) .* 3
    @test v[7] == ones(3) .* 4

    @test parent(v)[1] == ones(3) .* 1
    @test parent(v)[2] == ones(3) .* 2
    @test parent(v)[3] == ones(3) .* 3
    @test parent(v)[4] == ones(3) .* 4

    @test_throws BoundsError v[2]
    @test_throws BoundsError v[8]
    @inferred v[7]

    a1 = ones(3) .* 1
    a2 = ones(3) .* 2
    a3 = ones(3) .* 3
    a4 = ones(3) .* 4
    v = SparseContainer([a1, a2, a3, a4], (1, 3, 5, 7))
    @test v[1] == ones(3) .* 1
    @test v[3] == ones(3) .* 2
    @test v[5] == ones(3) .* 3
    @test v[7] == ones(3) .* 4

    @test parent(v)[1] == ones(3) .* 1
    @test parent(v)[2] == ones(3) .* 2
    @test parent(v)[3] == ones(3) .* 3
    @test parent(v)[4] == ones(3) .* 4

    @test_throws BoundsError v[2]
    @test_throws BoundsError v[8]
    @inferred v[7]
end
