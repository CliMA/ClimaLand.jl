using IndirectArrays, MappedArrays, OrderedCollections
using Test, FixedPointNumbers, Colors

@testset "values::AbstractVector" begin
    colors = [RGB(1,0,0) RGB(0,1,0);
              RGB(0,0,1) RGB(1,0,0)]
    index0 = [1 3;
              2 1]
    for indexT in (Int8, Int16, UInt8, UInt16)
        A = IndirectArray{indexT}(colors)
        @test eltype(A) == RGB{N0f8}
        @test size(A) == (2,2)
        @test ndims(A) == 2
        @test A[1,1] === A[1] === RGB(1,0,0)
        @test A[2,1] === A[2] === RGB(0,0,1)
        @test A[1,2] === A[3] === RGB(0,1,0)
        @test A[2,2] === A[4] === RGB(1,0,0)
        @test isa(eachindex(A), AbstractUnitRange)
        @test A.index == index0
    end
    x = IndirectArray(colors[:])
    @test x == IndirectArray{UInt8}(colors[:])
    xc = copy(x)
    @test x == xc
    xc[2], xc[3] = RGB(0,1,0), RGB(0,1,1)
    @test xc == IndirectArray([RGB(1,0,0), RGB(0,1,0), RGB(0,1,1), RGB(1,0,0)])
    @test append!(x, x) == IndirectArray([colors[:]; colors[:]])
    @test append!(x, IndirectArray([RGB(1,0,0), RGB(1,1,0)])) ==
        IndirectArray([colors[:]; colors[:]; [RGB(1,0,0), RGB(1,1,0)]])
    # Append with non-IndirectArray
    @test append!(IndirectArray(colors[:]), IndirectArray(colors[:])) ==
          append!(IndirectArray(colors[:]), colors[:])

    # Bounds checking upon construction
    index_ob = copy(index0)
    index_ob[1] = 5   # out-of-bounds
    unsafe_ia(idx, vals) = (@inbounds ret = IndirectArray(idx, vals); ret)
      safe_ia(idx, vals) = (ret = IndirectArray(idx, vals); ret)
    @test_throws BoundsError safe_ia(index_ob, colors[1:3])
    # This requires inlining, which means it fails on Travis since we turn
    # off inlining for better coverage stats
    # B = unsafe_ia(index_ob, colors)
    # @test_throws BoundsError B[1]
    # @test B[2] == RGB(0,0,1)

    # Non-Arrays
    a = [0.1  0.4;
        0.33 1.0]
    f(x) = round(Int, 99*x) + 1   # maps 0-1 to 1-100
    m = mappedarray(f, a)
    cmap = colormap("RdBu", 100)
    img = IndirectArray(m, cmap)
    @test img == [cmap[11] cmap[41];
                  cmap[34] cmap[100]]
end

@testset "values::AbstractDict" begin
    # With Dicts
    a = ['a' 'b';
        'q' 'j']
    v = freeze(Dict('a' => "apple", 'b' => "book", 'q' => "quit", 'j' => "jolly"))
    A = IndirectArray(a, v)
    @test A[1,1] == "apple"
    @test A[2,1] == "quit"
    @test A[1,2] == "book"
    @test A[2,2] == "jolly"
    @test_throws BoundsError IndirectArray(a, Dict('a' => "apple", 'b' => "book"))
end
