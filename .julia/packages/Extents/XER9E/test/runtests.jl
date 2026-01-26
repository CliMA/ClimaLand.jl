using Extents
using Test
using Dates
const E = Extent

ex1 = Extent(X=(1, 2), Y=(3, 4))
ex2 = Extent(Y=(3, 4), X=(1, 2))
ex3 = Extent(X=(1, 2), Y=(3, 4), Z=(5.0, 6.0))

struct HasExtent end
Extents.extent(::HasExtent) = Extent(X=(0, 1), Y=(0, 1))

@testset "getindex" begin
    @test ex3[1] == ex3[:X] == (1, 2)
    @test ex3[[:X, :Z]] == ex3[(:X, :Z)] == Extent{(:X, :Z)}(((1, 2), (5.0, 6.0)))
end

@testset "getproperty" begin
    @test ex3.X == (1, 2)
    @test propertynames(ex3) == (:X, :Y, :Z)
end

@testset "bounds" begin
    @test bounds(ex1) === (X=(1, 2), Y=(3, 4))
    @test bounds(ex2) === (Y=(3, 4), X=(1, 2))
    @test bounds(ex3) === (X=(1, 2), Y=(3, 4), Z=(5.0, 6.0))
end
@testset "NamedTuple" begin
    @test NamedTuple(ex1) === (X=(1, 2), Y=(3, 4))
    @test NamedTuple(ex2) === (Y=(3, 4), X=(1, 2))
    @test NamedTuple(ex3) === (X=(1, 2), Y=(3, 4), Z=(5.0, 6.0))
end

@testset "extent can be called on Extent, returning itself" begin
    @test extent(ex1) === ex1
end

@testset "Base julia equality works in any order for the same dimensions" begin
    # Extents equal themselves
    @test ex1 == ex1
    # Order doesn't matter
    @test ex1 == ex2
    # Extra dimensions are not equal
    @test ex1 != ex3
end

@testset "isapprox" begin
    ex4 = Extent(X=(1.00000000000001, 2), Y=(3, 4))
    ex5 = Extent(X=(1.1, 2), Y=(3, 4))
    @test ex1 ≈ ex1
    @test ex1 ≈ ex2
    @test ex1 ≈ ex4
    @test isapprox(ex1, ex5; atol=0.11)
end

@testset "keys and values are just like a NamedTuple" begin
    @test keys(ex1) == (:X, :Y)
    @test values(ex1) == ((1, 2), (3, 4))
end

@testset "union" begin
    a = E(X=(0.1, 0.5), Y=(1.0, 2.0))
    b = E(X=(2.1, 2.5), Y=(3.0, 4.0), Z=(0.0, 1.0))
    c = E(Z=(0.2, 2.0))
    @test Extents.union(a, b) == Extents.union(a, b, a) == E(X=(0.1, 2.5), Y=(1.0, 4.0))
    @test Extents.union(a, b; strict=true) === nothing
    @test Extents.union(a, c) === nothing

    # If either argument is nothing, return the other
    @test Extents.union(a, nothing) === a
    @test Extents.union(nothing, b) === b
    @test Extents.union(a, nothing; strict=true) === nothing
    @test Extents.union(nothing, b; strict=true) === nothing
    # If both arguments are nothing, return nothing
    @test Extents.union(nothing, nothing) === nothing
    @test Extents.union(nothing, nothing; strict=true) === nothing
end

@testset "covers" begin
    # An extent contains itself
    @test Extents.covers(E(X=(2, 3), Y=(3, 4)), E(X=(2, 3), Y=(3, 4))) == true
    # A larger extent covers a smaller one inside it
    @test Extents.covers(E(X=(0, 4), Y=(1, 5)), E(X=(2, 3), Y=(3, 4))) == true
    # Intersecting but not covers in one dimension
    @test Extents.covers(E(X=(0, 5), Y=(1, 3)), E(X=(2, 3), Y=(3, 4))) == false
    @test Extents.covers(E(X=(3, 5), Y=(1, 5)), E(X=(2, 3), Y=(3, 4))) == false
    # An extent with no interior is still covered
    @test Extents.covers(E(X=(0, 4), Y=(1, 5)), E(X=(2, 2), Y=(3, 4))) == true
    # Not containing in any dimensions
    @test Extents.covers(E(X=(0, 1), Y=(1, 2)), E(X=(2, 3), Y=(3, 4))) == false
    @test Extents.covers(E(X=(4, 5), Y=(5, 6)), E(X=(2, 3), Y=(3, 4))) == false
    # We just ignore missing dimensions
    @test Extents.covers(E(X=(2, 4), Y=(1, 4), Z=(1, 2)), E(X=(2, 3), Y=(3, 4))) == true
    @test Extents.covers(E(X=(2, 4), Y=(1, 4)), E(X=(2, 3), Y=(3, 4), Z=(1, 2))) == true
    # Except when `strict` is true
    @test Extents.covers(E(X=(2, 4), Y=(1, 4)), E(X=(2, 3), Y=(3, 4), Z=(1, 2)); strict=true) == false
    # When they are present *they can change the result*
    @test Extents.covers(E(X=(2, 4), Y=(1, 4), Z=(0, 1)), E(X=(2, 3), Y=(3, 4), Z=(1, 2))) == false
    # Nothing returns false
    @test Extents.covers(E(X=(2, 3), Y=(3, 4)), nothing) == false
    @test Extents.covers(nothing, E(X=(2, 3), Y=(3, 4))) == false
    @test Extents.covers(nothing, nothing) == false
    # Objects that don't have extents also return false
    @test Extents.covers(1, :x) == false
    # Objects that have extents can be used
    @test Extents.covers(HasExtent(), Extents.extent(HasExtent())) == true
end

@testset "coveredby" begin
    # An extent contains itself
    @test Extents.coveredby(E(X=(2, 3), Y=(3, 4)), E(X=(2, 3), Y=(3, 4))) == true
    # A larger extent coveredby a smaller one inside it
    @test Extents.contains(E(X=(0, 4), Y=(1, 5)), E(X=(2, 3), Y=(3, 4))) == true
    @test Extents.coveredby(E(X=(2, 3), Y=(3, 4)), E(X=(0, 4), Y=(1, 5))) == true
    # Intersecting but not coveredby in one dimension
    @test Extents.coveredby(E(X=(2, 3), Y=(3, 4)), E(X=(0, 5), Y=(1, 3))) == false
    @test Extents.coveredby(E(X=(2, 3), Y=(3, 4)), E(X=(3, 5), Y=(1, 5))) == false
    # An extent with no interior is still coveredby
    @test Extents.coveredby(E(X=(2, 2), Y=(3, 4)), E(X=(0, 4), Y=(1, 5))) == true
    # Not coveredby in any dimensions
    @test Extents.coveredby(E(X=(2, 3), Y=(3, 4)), E(X=(0, 1), Y=(1, 2))) == false
    @test Extents.coveredby(E(X=(4, 5), Y=(5, 6)), E(X=(1, 3), Y=(3, 4))) == false
    # We just ignore missing dimensions
    @test Extents.coveredby(E(X=(2, 3), Y=(3, 4)), E(X=(2, 4), Y=(1, 4), Z=(1, 2))) == true
    @test Extents.coveredby(E(X=(2, 3), Y=(3, 4), Z=(1, 2)), E(X=(2, 4), Y=(1, 4))) == true
    # Except when `strict` is true
    @test Extents.coveredby(E(X=(2, 3), Y=(3, 4), Z=(1, 2)), E(X=(2, 4), Y=(1, 4)); strict=true) == false
    # When they are present *they can change the result*
    @test Extents.coveredby(E(X=(2, 4), Y=(1, 4), Z=(0, 1)), E(X=(2, 3), Y=(3, 4), Z=(1, 2))) == false
    # Nothing returns false
    @test Extents.coveredby(E(X=(2, 3), Y=(3, 4)), nothing) == false
    @test Extents.coveredby(nothing, E(X=(2, 3), Y=(3, 4))) == false
    @test Extents.coveredby(nothing, nothing) == false
    # Objects that don't have extents also return false
    @test Extents.coveredby(1, :x) == false
    # Objects that have extents can be used
    @test Extents.coveredby(HasExtent(), Extents.extent(HasExtent())) == true
end

@testset "contains" begin
    # An extent contains itself
    @test Extents.contains(E(X=(2, 3), Y=(3, 4)), E(X=(2, 3), Y=(3, 4))) == true
    # A larger extent contains a smaller one inside it
    @test Extents.contains(E(X=(0, 4), Y=(1, 5)), E(X=(2, 3), Y=(3, 4))) == true
    # Intersecting but not contains in one dimension
    @test Extents.contains(E(X=(0, 5), Y=(1, 3)), E(X=(2, 3), Y=(3, 4))) == false
    @test Extents.contains(E(X=(3, 5), Y=(1, 5)), E(X=(2, 3), Y=(3, 4))) == false
    # An extent with no interior is not contained
    @test Extents.contains(E(X=(0, 4), Y=(1, 5)), E(X=(2, 2), Y=(3, 4))) == false
    # Not containing in any dimensions
    @test Extents.contains(E(X=(0, 1), Y=(1, 2)), E(X=(2, 3), Y=(3, 4))) == false
    @test Extents.contains(E(X=(4, 5), Y=(5, 6)), E(X=(2, 3), Y=(3, 4))) == false
    # We just ignore missing dimensions
    @test Extents.contains(E(X=(2, 4), Y=(1, 4), Z=(1, 2)), E(X=(2, 3), Y=(3, 4))) == true
    @test Extents.contains(E(X=(2, 4), Y=(1, 4)), E(X=(2, 3), Y=(3, 4), Z=(1, 2))) == true
    # Except when `strict` is true
    @test Extents.contains(E(X=(2, 4), Y=(1, 4)), E(X=(2, 3), Y=(3, 4), Z=(1, 2)); strict=true) == false
    # When they are present *they can change the result*
    @test Extents.contains(E(X=(2, 4), Y=(1, 4), Z=(0, 1)), E(X=(2, 3), Y=(3, 4), Z=(1, 2))) == false
    # Nothing returns false
    @test Extents.contains(E(X=(2, 3), Y=(3, 4)), nothing) == false
    @test Extents.contains(nothing, E(X=(2, 3), Y=(3, 4))) == false
    @test Extents.contains(nothing, nothing) == false
    # Objects that don't have extents also return false
    @test Extents.contains(1, :x) == false
    # Objects that have extents can be used
    @test Extents.contains(HasExtent(), Extents.extent(HasExtent())) == true
end

@testset "within" begin
    # An extent contains itself
    @test Extents.within(E(X=(2, 3), Y=(3, 4)), E(X=(2, 3), Y=(3, 4))) == true
    # A larger extent within a smaller one inside it
    @test Extents.contains(E(X=(0, 4), Y=(1, 5)), E(X=(2, 3), Y=(3, 4))) == true
    @test Extents.within(E(X=(2, 3), Y=(3, 4)), E(X=(0, 4), Y=(1, 5))) == true
    # Intersecting but not within in one dimension
    @test Extents.within(E(X=(2, 3), Y=(3, 4)), E(X=(0, 5), Y=(1, 3))) == false
    @test Extents.within(E(X=(2, 3), Y=(3, 4)), E(X=(3, 5), Y=(1, 5))) == false
    # An extent with no interior is not within
    @test Extents.within(E(X=(2, 2), Y=(3, 4)), E(X=(0, 4), Y=(1, 5))) == false
    # Not within in any dimensions
    @test Extents.within(E(X=(2, 3), Y=(3, 4)), E(X=(0, 1), Y=(1, 2))) == false
    @test Extents.within(E(X=(4, 5), Y=(5, 6)), E(X=(1, 3), Y=(3, 4))) == false
    # We just ignore missing dimensions
    @test Extents.within(E(X=(2, 3), Y=(3, 4)), E(X=(2, 4), Y=(1, 4), Z=(1, 2))) == true
    @test Extents.within(E(X=(2, 3), Y=(3, 4), Z=(1, 2)), E(X=(2, 4), Y=(1, 4))) == true
    # Except when `strict` is true
    @test Extents.within(E(X=(2, 3), Y=(3, 4), Z=(1, 2)), E(X=(2, 4), Y=(1, 4)); strict=true) == false
    # When they are present *they can change the result*
    @test Extents.within(E(X=(2, 4), Y=(1, 4), Z=(0, 1)), E(X=(2, 3), Y=(3, 4), Z=(1, 2))) == false
    # Nothing returns false
    @test Extents.within(E(X=(2, 3), Y=(3, 4)), nothing) == false
    @test Extents.within(nothing, E(X=(2, 3), Y=(3, 4))) == false
    @test Extents.within(nothing, nothing) == false
    # Objects that don't have extents also return false
    @test Extents.within(1, :x) == false
    # Objects that have extents can be used
    @test Extents.within(HasExtent(), Extents.extent(HasExtent())) == true
end

@testset "overlaps" begin
    x = E(X=(2, 3), Y=(3, 5))
    # These overlap 
    a = E(X=(0, 3), Y=(1, 4)) 
    @test Extents.overlaps(a, x) == Extents.overlaps(x, a) == true
    # On the line is not overlapping
    b = E(X=(0, 3), Y=(1, 3))
    @test Extents.overlaps(b, x) == Extents.overlaps(x, b) == false
    # Corner touching is also not overlapping
    c = E(X=(0, 2), Y=(1, 3))
    @test Extents.overlaps(c, x) == Extents.overlaps(x, c) == false
    # Same on the high side of the bounds
    d = E(X=(0, 6), Y=(5, 5))
    @test Extents.overlaps(d, x) == Extents.overlaps(x, d) == false
    e = E(X=(0.0, 2.1), Y=(3.0, 5.0)) 
    @test Extents.overlaps(e, x) == Extents.overlaps(x, e) == true
    # We just ignore missing dimensions
    @test Extents.overlaps(E(X=(0, 3), Y=(1, 4)), E(X=(2, 3), Y=(3, 5), Z=(9, 10))) == true
    @test Extents.overlaps(E(X=(0, 3), Y=(1, 4), Z=(9, 10)), E(X=(2, 3), Y=(3, 5))) == true
    # Except when `strict` is true
    @test Extents.overlaps(E(X=(0, 3), Y=(1, 4), Z=(9, 10)), E(X=(2, 3), Y=(3, 5)); strict=true) == false
    # When they are present in both *they can change the result*
    @test Extents.overlaps(E(X=(0, 3), Y=(1, 4), Z=(1, 2)), E(X=(2, 3), Y=(3, 5), Z=(9, 10))) == false
    # Nothing returns false
    @test Extents.overlaps(E(X=(2, 3), Y=(3, 4)), nothing) == false
    @test Extents.overlaps(nothing, E(X=(2, 3), Y=(3, 4))) == false
    @test Extents.overlaps(nothing, nothing) == false
    # Objects that don't have extents also return false
    @test Extents.overlaps(1, :x) == false
    # Objects that have extents can be used
    @test Extents.overlaps(HasExtent(), Extent(X=(0.5, 1.0), Y=(0.5, 1.5))) == true
end

@testset "touches" begin
    # These touch at the corner of (X=2.0, Y=3.0)
    @test Extents.touches(E(X=(0.0, 2.0), Y=(1.0, 3.0)), E(X=(2.0, 3.0), Y=(3.0, 5.0))) == true
    # These touch on the side (X=2.0, Y=3.0) - (X=2.0, Y=5.0) 
    @test Extents.touches(E(X=(0.0, 2.0), Y=(3.0, 5.0)), E(X=(2.0, 3.0), Y=(3.0, 5.0)))
    # Overlapping is not touching
    @test Extents.touches(E(X=(0.0, 3.0), Y=(1.0, 4.0)), E(X=(2.0, 3.0), Y=(3.0, 5.0))) == false
    # We just ignore missing dimensions
    @test Extents.touches(E(X=(0, 3), Y=(1, 3)), E(X=(2, 3), Y=(3, 5), Z=(9, 10))) == true
    @test Extents.touches(E(X=(0, 3), Y=(1, 3), Z=(9, 10)), E(X=(2, 3), Y=(3, 5))) == true
    # Except when `strict` is true
    @test Extents.touches(E(X=(0, 3), Y=(1, 3)), E(X=(2, 3), Y=(3, 5), Z=(9, 10)); strict=true) == false
    @test Extents.touches(E(X=(0, 3), Y=(1, 3), Z=(9, 10)), E(X=(2, 3), Y=(3, 5)); strict=true) == false
    # When they are present in both they can change the result
    @test Extents.touches(E(X=(0, 3), Y=(1, 4), Z=(1, 2)), E(X=(2, 3), Y=(3, 5), Z=(9, 10))) == false
    # Nothing returns false
    @test Extents.touches(E(X=(2, 3), Y=(3, 4)), nothing) == false
    @test Extents.touches(nothing, E(X=(2, 3), Y=(3, 4))) == false
    @test Extents.touches(nothing, nothing) == false
    # Objects that don't have extents also return false
    @test Extents.touches(1, :x) == false
    # Objects that have extents can be used
    @test Extents.touches(HasExtent(), Extent(X=(1.0, 1.0), Y=(1.0, 1.0))) == true
end

@testset "equals" begin
    # Matching bounds are `equal`
    @test Extents.equals(E(X=(2, 3), Y=(1, 4)), E(X=(2, 3), Y=(1, 4))) == true
    # Their order doesn't matter
    @test Extents.equals(E(X=(2, 3), Y=(1, 4)), E(Y=(1, 4), X=(2, 3))) == true
    # We just ignore missing dimensions
    @test Extents.equals(E(X=(2, 3), Y=(1, 4), Z=(1, 2)), E(X=(2, 3), Y=(1, 4))) == true
    # Except when `strict` is true
    @test Extents.equals(E(X=(2, 3), Y=(1, 4), Z=(1, 2)), E(X=(2, 3), Y=(1, 4)); strict=true) == false
    # Adding a dimension can change the result
    @test Extents.equals(E(X=(2, 3), Y=(1, 4), Z=(1, 2)), E(X=(2, 3), Y=(1, 4), Z=(0, 2))) == false
    # Nothing returns false
    @test Extents.equals(E(X=(2, 3), Y=(3, 4)), nothing) == false
    @test Extents.equals(nothing, E(X=(2, 3), Y=(3, 4))) == false
    @test Extents.equals(nothing, nothing) == false
    # Objects that don't have extents also return false
    @test Extents.equals(1, :x) == false
    # Objects that have extents can be used
    @test Extents.equals(HasExtent(), Extents.extent(HasExtent())) == true
    @test Extents.equals(Extents.extent(HasExtent()), HasExtent()) == true
end

@testset "intersection/intersects/disjoint" begin
    a = E(X=(0.1, 0.5), Y=(1.0, 2.0))
    b = E(X=(2.1, 2.5), Y=(3.0, 4.0))
    c = E(X=(0.4, 2.5), Y=(1.5, 4.0), Z=(0.0, 1.0))
    d = E(X=(0.2, 0.45))
    e = E(A=(0.0, 1.0))
    # a and b don't intersect
    @test Extents.intersects(a, b) == false
    @test Extents.intersects(b, a) == false
    # c and d do, despite the extra Z dimension
    @test Extents.intersects(c, a) == true
    @test Extents.intersects(a, c) == true
    # Unless strict is true
    @test Extents.intersects(a, c; strict=true) == false
    @test Extents.intersects(c, a; strict=true) == false
    # a and d do despite missing Y
    @test Extents.intersects(a, d) == true
    @test Extents.intersects(d, a) == true
    # Unless strict is true
    @test Extents.intersects(a, d; strict=true) == false
    @test Extents.intersects(d, a; strict=true) == false
    # If there are no dimensions, intersects is false
    @test Extents.intersects(a, e) == false
    @test Extents.intersects(e, a) == false
    # Nothing returns false
    @test Extents.intersects(E(X=(2, 3), Y=(3, 4)), nothing) == false
    @test Extents.intersects(nothing, E(X=(2, 3), Y=(3, 4))) == false
    @test Extents.intersects(nothing, nothing) == false
    # Objects that don't have extents also return false
    @test Extents.intersects(1, :x) == false
    # Objects that have extents can be used
    @test Extents.intersects(HasExtent(), Extents.extent(HasExtent())) == true
    @test Extents.intersects(Extents.extent(HasExtent()), HasExtent()) == true
    
    # a and b are disjoint
    @test Extents.disjoint(a, b) == true
    @test Extents.disjoint(b, a) == true
    # c and d do, despite the extra Z dimension
    @test Extents.disjoint(c, a) == false
    @test Extents.disjoint(a, c) == false
    # Unless strict is true
    @test Extents.disjoint(a, c; strict=true) == true
    @test Extents.disjoint(c, a; strict=true) == true
    # a and d do despite missing Y
    @test Extents.disjoint(a, d) == false
    @test Extents.disjoint(d, a) == false
    # Unless strict is true
    @test Extents.disjoint(a, d; strict=true) == true
    @test Extents.disjoint(d, a; strict=true) == true
    # If there are no dimensions, disjoint is false
    @test Extents.disjoint(a, e) == true
    @test Extents.disjoint(e, a) == true
    # Nothing returns true
    @test Extents.disjoint(E(X=(2, 3), Y=(3, 4)), nothing) == true
    @test Extents.disjoint(nothing, E(X=(2, 3), Y=(3, 4))) == true
    @test Extents.disjoint(nothing, nothing) == true
    # Objects that don't have extents also return false
    @test Extents.disjoint(1, :x) == true
    # Objects that have extents can be used
    @test Extents.disjoint(HasExtent(), Extent(X=(2, 3), Y=(4, 5))) == true
    @test Extents.disjoint(Extent(X=(2, 3), Y=(4, 5)), HasExtent()) == true

    # a and b do not intersect
    @test Extents.intersection(a, b) === nothing
    @test Extents.intersection(b, a) === nothing
    # a and c do, if we ignore the extra Z dimension
    @test Extents.intersection(a, c) == 
          Extents.intersection(c, a) == 
          Extents.intersection(a, c, a) == 
          Extents.intersection(a, c, a, c) == 
          Extent(X=(0.4, 0.5), Y=(1.5, 2.0))
    @test Extents.intersection(c, a) == Extents.intersection(c, a) == Extent(X=(0.4, 0.5), Y=(1.5, 2.0))
    # Unless strict is true
    @test Extents.intersection(a, c; strict=true) === nothing
    @test Extents.intersection(c, a; strict=true) === nothing

    # a and d intersect on X
    @test Extents.intersection(a, d) == Extents.intersection(d, a) == Extent(X=(0.2, 0.45))
    # If there are no dimensions, there are no intersections
    @test Extents.intersection(a, e) === nothing
    @test Extents.intersection(e, a) === nothing

    # If either argument is nothing, return nothing
    @test Extents.intersection(a, nothing) === nothing
    @test Extents.intersection(nothing, nothing) === nothing
    @test Extents.intersection(nothing, b) === nothing
end

@testset "buffer" begin
    a = Extent(X=(0.1, 0.5), Y=(1.0, 2.0))
    b = Extent(Lat=(0.1, 0.5), Lon=(1.0, 2.0), Elev=(-3, 4))
    c = Extent(X=(0.1, 0.5), Y=(1.0, 2.0), Ti=(DateTime(2000, 1, 1), DateTime(2020, 1, 1)))
    @test Extents.buffer(a,(H=0,)) == a
    @test Extents.buffer(a, (X=1, Y=2)) == Extent(X=(-0.9, 1.5), Y=(-1.0, 4.0))
    @test Extents.buffer(b, (Lat=2, Lon=1)) == Extent(Lat=(-1.9, 2.5), Lon=(0.0, 3.0), Elev=(-3, 4))
    @test Extents.buffer(c, (X=2, Y=2, Ti=Year(1))) == Extent(X=(-1.9, 2.5), Y=(-1.0, 4.0), Ti=(DateTime("1999-01-01T00:00:00"), DateTime("2021-01-01T00:00:00")))
end

@testset "grow" begin
    a = Extent(X=(1.0, 2.0), Y=(4.0, 8.0))
    @test Extents.grow(a, 0.5) == Extent(X=(0.5, 2.5), Y=(2.0, 10.0))
    @test Extents.grow(a, (X=0.1, Y=0.5)) == Extent(X=(0.9, 2.1), Y=(2.0, 10.0))
end

@testset "Single argument predicates" begin
    a = Extent(X=(0.1, 0.5), Y=(1.0, 2.0))
    b = Extent(X=(0.2, 0.4), Y=(2.0, 5.0), Z=(1, 2))
    @test Extents.intersects(b)(a) == Extents.intersects(a, b) == true
    @test Extents.disjoint(b)(a) == Extents.disjoint(a, b) == false
    @test Extents.touches(b)(a) == Extents.touches(a, b) == true
    @test Extents.equals(b)(a) == Extents.equals(a, b) == false
    @test Extents.covers(b)(a) == Extents.covers(a, b) == false
    @test Extents.coveredby(b)(a) == Extents.coveredby(a, b) == false
    @test Extents.contains(b)(a) == Extents.contains(a, b) == false
    @test Extents.within(b)(a) == Extents.within(a, b) == false
    @test Extents.touches(b)(a) == Extents.touches(a, b) == true

    @testset "strict keyword still works" begin
        @test Extents.intersects(b; strict=true)(a) == 
            Extents.intersects(b)(a; strict=true) ==
            Extents.intersects(a, b; strict=true) == false
    end
end