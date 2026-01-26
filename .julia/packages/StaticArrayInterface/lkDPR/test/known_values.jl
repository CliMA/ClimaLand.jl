@testset "known values" begin
    CI = CartesianIndices((2, 2))

    @test isnothing(@inferred(StaticArrayInterface.known_first(typeof(1:4))))
    @test isone(@inferred(StaticArrayInterface.known_first(Base.OneTo(4))))
    @test isone(@inferred(StaticArrayInterface.known_first(typeof(Base.OneTo(4)))))
    @test @inferred(StaticArrayInterface.known_first(typeof(CI))) == CartesianIndex(1, 1)
    @test @inferred(StaticArrayInterface.known_first(typeof(CI))) == CartesianIndex(1, 1)

    @test isnothing(@inferred(StaticArrayInterface.known_last(1:4)))
    @test isnothing(@inferred(StaticArrayInterface.known_last(typeof(1:4))))
    @test @inferred(StaticArrayInterface.known_last(typeof(CI))) === nothing

    @test isnothing(@inferred(StaticArrayInterface.known_step(typeof(1:0.2:4))))
    @test isone(@inferred(StaticArrayInterface.known_step(1:4)))
    @test isone(@inferred(StaticArrayInterface.known_step(typeof(1:4))))
    @test isone(@inferred(StaticArrayInterface.known_step(typeof(Base.Slice(1:4)))))
    @test isone(@inferred(StaticArrayInterface.known_step(typeof(view(1:4, 1:2)))))
end