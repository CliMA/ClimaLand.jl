
@testset "Range Constructors" begin
    @test @inferred(static(1):static(10)) == 1:10
    @test @inferred(Static.SUnitRange{1, 10}()) == 1:10
    @test @inferred(static(1):static(2):static(10)) == 1:2:10
    @test @inferred(1:static(2):static(10)) == 1:2:10
    @test @inferred(static(1):static(2):10) == 1:2:10
    @test @inferred(static(1):2:static(10)) == 1:2:10
    @test @inferred(1:2:static(10)) == 1:2:10
    @test @inferred(1:static(2):10) == 1:2:10
    @test @inferred(static(1):2:10) == 1:2:10
    @test @inferred(static(1):UInt(10)) === static(1):10
    @test @inferred(UInt(1):static(1):static(10)) === 1:static(10)
    @test Static.SUnitRange(1, 10) == 1:10
    @test @inferred(Static.OptionallyStaticUnitRange{Int, Int}(1:10)) == 1:10
    @test @inferred(Static.OptionallyStaticUnitRange(1:10)) == 1:10 ==
          @inferred(Static.OptionallyStaticUnitRange(Static.OptionallyStaticUnitRange(1:10)))

    sr = Static.OptionallyStaticStepRange(static(1), static(1), static(1))
    @test @inferred(Static.OptionallyStaticStepRange(sr)) == sr == 1:1:1
    @test @inferred(Static.OptionallyStaticStepRange(static(1), 1, UInt(10))) ==
          static(1):1:10 == Static.SOneTo(10)
    @test @inferred(Static.OptionallyStaticStepRange(UInt(1), 1, static(10))) ==
          static(1):1:10
    @test @inferred(Static.OptionallyStaticStepRange(1:10)) == 1:1:10

    @test_throws ArgumentError Static.OptionallyStaticUnitRange(1:2:10)
    @test_throws ArgumentError Static.OptionallyStaticUnitRange{Int, Int}(1:2:10)
    @test_throws ArgumentError Static.OptionallyStaticStepRange(1, 0, 10)
    @test_throws ArgumentError Static.OptionallyStaticStepRange(1, StaticInt(0), 10)

    @test @inferred(static(1):static(1):static(10)) ===
          Static.OptionallyStaticUnitRange(static(1), static(10))
    @test @inferred(static(1):static(1):10) ===
          Static.OptionallyStaticUnitRange(static(1), 10)
    @test @inferred(1:static(1):10) === Static.OptionallyStaticUnitRange(1, 10)
    @test length(static(-1):static(-1):static(-10)) == 10 ==
          lastindex(static(-1):static(-1):static(-10))

    @test UnitRange(Static.OptionallyStaticUnitRange(static(1), static(10))) ===
          UnitRange(1, 10)
    @test UnitRange{Int}(Static.OptionallyStaticUnitRange(static(1), static(10))) ===
          UnitRange(1, 10)

    @test AbstractUnitRange{Int}(Static.OptionallyStaticUnitRange(
        static(1), static(10))) isa
          Static.OptionallyStaticUnitRange
    @test AbstractUnitRange{UInt}(Static.OptionallyStaticUnitRange(
        static(1), static(10))) isa
          Base.OneTo
    @test AbstractUnitRange{UInt}(Static.OptionallyStaticUnitRange(
        static(2), static(10))) isa
          UnitRange

    @test @inferred((static(1):static(10))[static(2):static(3)]) === static(2):static(3)
    @test @inferred((static(1):static(10))[static(2):3]) === static(2):3
    @test @inferred((static(1):static(10))[2:3]) === 2:3
    @test @inferred((1:static(10))[static(2):static(3)]) === 2:3

    @test Base.checkindex(Bool, static(1):static(10), static(1):static(5))
    @test -(static(1):static(10)) === static(-1):static(-1):static(-10)

    @test reverse(static(1):static(10)) === static(10):static(-1):static(1)
    @test reverse(static(1):static(2):static(9)) === static(9):static(-2):static(1)
end

@testset "range properties" begin
    x = static(1):static(2):static(9)
    @test getproperty(x, :start) === first(x)
    @test getproperty(x, :step) === step(x)
    @test getproperty(x, :stop) === last(x)
    @test_throws ErrorException getproperty(x, :foo)
end

@testset "iterate" begin
    @test iterate(static(1):static(5)) === (1, 1)
    @test iterate(static(1):static(5), 1) === (2, 2)
    @test iterate(static(1):static(5), 5) === nothing
    @test iterate(static(2):static(5), 5) === nothing
    @test iterate(static(1):static(2):static(9), 1) === (3, 3)
    @test iterate(static(1):static(2):static(9), 9) === nothing
    # make sure single length ranges work correctly
    @test iterate(static(2):static(3):static(2))[1] === 2 ===
          (static(2):static(3):static(2))[1]
    @test iterate(static(2):static(3):static(2), 2) === nothing
end

# CartesianIndices
CI = CartesianIndices((static(1):static(2), static(1):static(2)))

@testset "length" begin
    @test @inferred(length(Static.OptionallyStaticUnitRange(1, 0))) == 0
    @test @inferred(length(Static.OptionallyStaticUnitRange(1, 10))) == 10
    @test @inferred(length(Static.OptionallyStaticUnitRange(static(1), 10))) == 10
    @test @inferred(length(Static.OptionallyStaticUnitRange(static(0), 10))) == 11
    @test @inferred(length(Static.OptionallyStaticUnitRange(static(1), static(10)))) == 10
    @test @inferred(length(Static.OptionallyStaticUnitRange(static(0), static(10)))) == 11

    @test @inferred(length(static(1):static(2):static(0))) == 0
    @test @inferred(length(static(0):static(-2):static(1))) == 0

    @test @inferred(length(Static.OptionallyStaticStepRange(static(1), 2, 10))) == 5
    @test @inferred(length(Static.OptionallyStaticStepRange(static(1), static(1),
        static(10)))) == 10
    @test @inferred(length(Static.OptionallyStaticStepRange(static(2), static(1),
        static(10)))) == 9
    @test @inferred(length(Static.OptionallyStaticStepRange(static(2), static(2),
        static(10)))) == 5
end

@test @inferred(getindex(static(1):10, Base.Slice(static(1):10))) === static(1):10
@test @inferred(getindex(Static.OptionallyStaticUnitRange(static(1), 10), 1)) == 1
@test @inferred(getindex(Static.OptionallyStaticUnitRange(static(0), 10), 1)) == 0
@test_throws BoundsError getindex(Static.OptionallyStaticUnitRange(static(1), 10), 0)
@test_throws BoundsError getindex(Static.OptionallyStaticStepRange(static(1), 2, 10), 0)
@test_throws BoundsError getindex(Static.OptionallyStaticUnitRange(static(1), 10), 11)
@test_throws BoundsError getindex(Static.OptionallyStaticStepRange(static(1), 2, 10), 11)

@test Static.static_first(Base.OneTo(one(UInt))) === static(1)
@test Static.static_step(Base.OneTo(one(UInt))) === static(1)

@test @inferred(eachindex(static(-7):static(7))) === static(1):static(15)
@test @inferred((static(-7):static(7))[first(eachindex(static(-7):static(7)))]) == -7

@test @inferred(firstindex(128:static(-1):1)) == 1

@test identity.(static(1):5) isa Vector{Int}
@test (static(1):5) .+ (1:3)' isa Matrix{Int}
@test similar(Array{Int}, (static(1):(4),)) isa Vector{Int}
@test similar(Array{Int}, (static(1):(4), Base.OneTo(4))) isa Matrix{Int}
@test similar(Array{Int}, (Base.OneTo(4), static(1):(4))) isa Matrix{Int}

@test Base.to_shape(static(1):10) == 10
@test Base.to_shape(Base.Slice(static(1):10)) == 10
@test Base.axes1(Base.Slice(static(1):10)) === static(1):10
@test axes(Base.Slice(static(1):10)) === (static(1):10,)
@test isa(axes(Base.Slice(static(0):static(1):10))[1], Base.IdentityUnitRange)
@test isa(Base.axes1(Base.Slice(static(0):static(1):10)), Base.IdentityUnitRange)

@test Base.Broadcast.axistype(static(1):10, static(1):10) === Base.OneTo(10)
@test Base.Broadcast.axistype(Base.OneTo(10), static(1):10) === Base.OneTo(10)

@testset "static_promote(::AbstractRange, ::AbstractRange)" begin
    ur1 = static(1):10
    ur2 = 1:static(10)
    @test @inferred(static_promote(ur1, ur2)) === static(1):static(10)
    @test static_promote(Base.Slice(ur1), Base.Slice(ur2)) ===
          Base.Slice(static(1):static(10))
    sr1 = static(1):2:10
    sr2 = static(1):static(2):static(10)
    @test @inferred(static_promote(sr1, sr2)) === sr2
end

@testset "n-last/first" begin
    ur = static(2):static(10)
    sr = static(2):static(2):static(10)
    n3 = static(3)
    @test @inferred(first(ur, n3)) === static(2):static(4)
    @test @inferred(last(ur, n3)) === static(8):static(10)
    @test @inferred(first(sr, n3)) === static(2):static(2):static(6)
    @test @inferred(last(sr, n3)) === static(5):static(2):static(9)
end
