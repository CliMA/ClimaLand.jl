@testset "insert/deleteat" begin
    @test @inferred(StaticArrayInterface.insert([1,2,3], 2, -2)) == [1, -2, 2, 3]
    @test @inferred(StaticArrayInterface.deleteat([1, 2, 3], 2)) == [1, 3]

    @test @inferred(StaticArrayInterface.deleteat([1, 2, 3], [1, 2])) == [3]
    @test @inferred(StaticArrayInterface.deleteat([1, 2, 3], [1, 3])) == [2]
    @test @inferred(StaticArrayInterface.deleteat([1, 2, 3], [2, 3])) == [1]

    @test @inferred(StaticArrayInterface.insert((2,3,4), 1, -2)) == (-2, 2, 3, 4)
    @test @inferred(StaticArrayInterface.insert((2,3,4), 2, -2)) == (2, -2, 3, 4)
    @test @inferred(StaticArrayInterface.insert((2,3,4), 3, -2)) == (2, 3, -2, 4)

    @test @inferred(StaticArrayInterface.deleteat((2, 3, 4), 1)) == (3, 4)
    @test @inferred(StaticArrayInterface.deleteat((2, 3, 4), 2)) == (2, 4)
    @test @inferred(StaticArrayInterface.deleteat((2, 3, 4), 3)) == (2, 3)
    @test StaticArrayInterface.deleteat((1, 2, 3), [1, 2]) == (3,)
end