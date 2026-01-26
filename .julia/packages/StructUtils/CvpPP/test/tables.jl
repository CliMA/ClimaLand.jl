using StructUtils, Tables, Test

@testset "tables.jl" begin
    tbl = (a=[1, 2, 3], b=[4, 5, 6])
    tbl2 = [StructUtils.make(@NamedTuple{a::Int, b::Int}, row) for row in Tables.rows(tbl)]
    @test length(tbl2) == 3
    @test tbl2[1].a == 1
    @test tbl2[1].b == 4
    @test tbl2[2].a == 2
    @test tbl2[2].b == 5
    @test tbl2[3].a == 3
    @test tbl2[3].b == 6
end