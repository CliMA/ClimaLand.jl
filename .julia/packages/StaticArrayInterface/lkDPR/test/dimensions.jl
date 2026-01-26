
###
### define wrapper with StaticArrayInterface.dimnames
###

@testset "order_named_inds" begin
    n1 = (static(:x),)
    n2 = (n1..., static(:y))
    n3 = (n2..., static(:z))
    @test @inferred(StaticArrayInterface.find_all_dimnames(n1, (), (), :)) == ()
    @test @inferred(StaticArrayInterface.find_all_dimnames(n1, (static(:x),), (2,), :)) == (2,)
    @test @inferred(StaticArrayInterface.find_all_dimnames(n2, (static(:x),), (2,), :)) == (2, :)
    @test @inferred(StaticArrayInterface.find_all_dimnames(n2, (static(:y),), (2,), :)) == (:, 2)
    @test @inferred(StaticArrayInterface.find_all_dimnames(n2, (static(:y), static(:x)), (20, 30), :)) == (30, 20)
    @test @inferred(StaticArrayInterface.find_all_dimnames(n2, (static(:x), static(:y)), (30, 20), :)) == (30, 20)
    @test @inferred(StaticArrayInterface.find_all_dimnames(n3, (static(:x), static(:y)), (30, 20), :)) == (30, 20, :)

    @test_throws ErrorException StaticArrayInterface.find_all_dimnames(n2, (static(:x), static(:y), static(:z)), (30, 20, 40), :)
end

@testset "StaticArrayInterface.dimnames" begin
    d = (static(:x), static(:y))
    x = NamedDimsWrapper(d, ones(Int64, 2, 2))
    y = NamedDimsWrapper((static(:x),), ones(2))
    z = NamedDimsWrapper((:x, static(:y)), ones(2))
    r1 = reinterpret(Int8, x)
    r2 = reinterpret(reshape, Int8, x)
    r3 = reinterpret(reshape, Complex{Int}, x)
    r4 = reinterpret(reshape, Float64, x)
    w = Wrapper(x)
    dnums = ntuple(+, length(d))
    @test @inferred(StaticArrayInterface.has_dimnames(x)) == true
    @test @inferred(StaticArrayInterface.has_dimnames(z)) == true
    @test @inferred(StaticArrayInterface.has_dimnames(ones(2, 2))) == false
    @test @inferred(StaticArrayInterface.has_dimnames(Array{Int,2})) == false
    @test @inferred(StaticArrayInterface.has_dimnames(typeof(x))) == true
    @test @inferred(StaticArrayInterface.has_dimnames(typeof(view(x, :, 1, :)))) == true
    @test @inferred(StaticArrayInterface.dimnames(x)) === d
    @test @inferred(StaticArrayInterface.dimnames(w)) === d
    @test @inferred(StaticArrayInterface.dimnames(r1)) === d
    @test @inferred(StaticArrayInterface.dimnames(r2)) === (static(:_), d...)
    @test @inferred(StaticArrayInterface.dimnames(r3)) === Base.tail(d)
    @test @inferred(StaticArrayInterface.dimnames(r4)) === d
    @test @inferred(StaticArrayInterface.StaticArrayInterface.dimnames(z)) === (:x, static(:y))
    @test @inferred(StaticArrayInterface.dimnames(parent(x))) === (static(:_), static(:_))
    @test @inferred(StaticArrayInterface.dimnames(reshape(x, (1, 4)))) === d
    @test @inferred(StaticArrayInterface.dimnames(reshape(x, :))) === (static(:_),)
    @test @inferred(StaticArrayInterface.dimnames(x')) === reverse(d)
    @test @inferred(StaticArrayInterface.dimnames(y')) === (static(:_), static(:x))
    @test @inferred(StaticArrayInterface.dimnames(PermutedDimsArray(x, (2, 1)))) === reverse(d)
    @test @inferred(StaticArrayInterface.dimnames(PermutedDimsArray(x', (2, 1)))) === d
    @test @inferred(StaticArrayInterface.dimnames(view(x, :, 1))) === (static(:x),)
    @test @inferred(StaticArrayInterface.dimnames(view(x, :, 1)')) === (static(:_), static(:x))
    @test @inferred(StaticArrayInterface.dimnames(view(x, :, :, :))) === (static(:x), static(:y), static(:_))
    @test @inferred(StaticArrayInterface.dimnames(view(x, :, 1, :))) === (static(:x), static(:_))
    # multidmensional indices
    @test @inferred(StaticArrayInterface.dimnames(view(x, ones(Int, 2, 2), 1))) === (static(:_), static(:_))
    @test @inferred(StaticArrayInterface.dimnames(view(x, [CartesianIndex(1,1), CartesianIndex(1,1)]))) === (static(:_),)

    @test @inferred(StaticArrayInterface.dimnames(x, Static.One())) === static(:x)
    @test @inferred(StaticArrayInterface.dimnames(parent(x), Static.One())) === static(:_)
    @test @inferred(StaticArrayInterface.known_dimnames(Iterators.flatten(1:10))) === (:_,)
    @test @inferred(StaticArrayInterface.known_dimnames(Iterators.flatten(1:10), static(1))) === :_
    # multidmensional indices
    @test @inferred(StaticArrayInterface.known_dimnames(view(x, ones(Int, 2, 2), 1))) === (:_, :_)
    @test @inferred(StaticArrayInterface.known_dimnames(view(x, [CartesianIndex(1,1), CartesianIndex(1,1)]))) === (:_,)

    @test @inferred(StaticArrayInterface.known_dimnames(z)) === (nothing, :y)
    @test @inferred(StaticArrayInterface.known_dimnames(reshape(x, (1, 4)))) === (:x, :y)
    @test @inferred(StaticArrayInterface.known_dimnames(r1)) === (:x, :y)
    @test @inferred(StaticArrayInterface.known_dimnames(r2)) === (:_, :x, :y)
    @test @inferred(StaticArrayInterface.known_dimnames(r3)) === (:y,)
    @test @inferred(StaticArrayInterface.known_dimnames(r4)) === (:x, :y)
    @test @inferred(StaticArrayInterface.known_dimnames(w)) === (:x, :y)
    @test @inferred(StaticArrayInterface.known_dimnames(reshape(x, :))) === (:_,)
    @test @inferred(StaticArrayInterface.known_dimnames(view(x, :, 1)')) === (:_, :x)
end

@testset "to_dims" begin
    x = NamedDimsWrapper(static((:x, :y)), ones(2, 2))
    y = NamedDimsWrapper(static((:x, :y, :a, :b, :c, :d)), ones(6))

    @test @inferred(StaticArrayInterface.to_dims(x, :)) == Colon()
    @test @inferred(StaticArrayInterface.to_dims(x, 1)) == 1
    @testset "small case" begin
        @test @inferred(StaticArrayInterface.to_dims(x, (:x, :y))) == (1, 2)
        @test @inferred(StaticArrayInterface.to_dims(x, (:y, :x))) == (2, 1)
        @test @inferred(StaticArrayInterface.to_dims(x, :x)) == 1
        @test @inferred(StaticArrayInterface.to_dims(x, :y)) == 2
        @test_throws DimensionMismatch StaticArrayInterface.to_dims(x, static(:z))  # not found
        @test_throws DimensionMismatch StaticArrayInterface.to_dims(x, :z)  # not found
    end

    @testset "large case" begin
        @test @inferred(StaticArrayInterface.to_dims(y, :x)) == 1
        @test @inferred(StaticArrayInterface.to_dims(y, :a)) == 3
        @test @inferred(StaticArrayInterface.to_dims(y, :d)) == 6
        @test_throws DimensionMismatch StaticArrayInterface.to_dims(y, :z) # not found
    end
end

@testset "methods accepting StaticArrayInterface.dimnames" begin
    d = (static(:x), static(:y))
    x = NamedDimsWrapper(d, ones(2, 2))
    y = NamedDimsWrapper((static(:x),), ones(2))
    @test @inferred(size(x, first(d))) == size(parent(x), 1)
    @test @inferred(StaticArrayInterface.static_axes(y')) == (static(1):static(1), Base.OneTo(2))
    @test @inferred(axes(x, first(d))) == axes(parent(x), 1)
    @test strides(x, :x) == StaticArrayInterface.static_strides(parent(x))[1]
    @test @inferred(StaticArrayInterface.axes_types(x, static(:x))) <: Base.OneTo{Int}
    @test StaticArrayInterface.axes_types(x, :x) <: Base.OneTo{Int}
    @test @inferred(StaticArrayInterface.axes_types(LinearIndices{2,NTuple{2,Base.OneTo{Int}}})) <: NTuple{2,Base.OneTo{Int}}
    CI = CartesianIndices{2,Tuple{Base.OneTo{Int},UnitRange{Int}}}
    @test @inferred(StaticArrayInterface.axes_types(CI, static(1))) <: Base.OneTo{Int}

    x[x=1] = [2, 3]
    @test @inferred(getindex(x, x=1)) == [2, 3]
    y = NamedDimsWrapper((:x, static(:y)), ones(2, 2))
    # FIXME this doesn't correctly infer the output because it can't infer
    @test getindex(y, x=1) == [1, 1]
end
