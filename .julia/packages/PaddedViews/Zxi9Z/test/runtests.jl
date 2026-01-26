using OffsetArrays
using Test
ambs = detect_ambiguities(Base, Core)  # in case these have ambiguities of their own
using PaddedViews
using PaddedViews: filltype, IdentityUnitRange
@testset "ambiguities" begin
    @test isempty(setdiff(detect_ambiguities(PaddedViews, Base, Core), ambs))
end

if VERSION >= v"1.2"
    # array summary changes after v1.2
    using Documenter
    doctest(PaddedViews, manual = false)
end

@testset "PaddedView" begin
    for n = 0:5
        a = @inferred(PaddedView(0, ones(Int,ntuple(d->1,n)), ntuple(x->x+1,n)))
        @test axes(a) == ntuple(x->1:x+1,n)
        @test axes(a, ndims(a)+1) == Base.OneTo(1)
        @test @inferred(a[1]) == 1
        n > 0 && @test @inferred(a[2]) == 0
        @test @inferred(a[ntuple(x->1,n)...]) == 1
        n > 0 && @test @inferred(a[2, ntuple(x->1,n-1)...]) == 0
    end
    a0 = reshape([3])
    a = @inferred(PaddedView(-1, a0, ()))
    @test axes(a) == ()
    @test ndims(a) == 0
    @test a[] == 3

    a = reshape(1:9, 3, 3)
    A = @inferred(PaddedView(0.0, a, (Base.OneTo(4), Base.OneTo(5))))
    @test eltype(A) == Int
    @test ndims(A) == 2
    @test size(A) === (4,5)
    @test @inferred(axes(A)) === (Base.OneTo(4), Base.OneTo(5))
    @test @inferred(axes(A, 3)) === Base.OneTo(1)
    @test A == [1 4 7 0 0;
                2 5 8 0 0;
                3 6 9 0 0;
                0 0 0 0 0]
    A32 = @inferred(PaddedView{Float32}(0.0, a, (Base.OneTo(4), Base.OneTo(5))))
    @test eltype(A32) == Float32
    @test A32 == A
    @test A32[1, 1] === 1.0f0

    A = @inferred(PaddedView(0.0, a, (0:4, -1:5)))
    @test eltype(A) == Int
    @test ndims(A) == 2
    @test size(A) == (5, 7)
    @test @inferred(axes(A)) === map(IdentityUnitRange, (0:4, -1:5))
    @test @inferred(axes(A, 3)) == IdentityUnitRange(1:1)
    @test A == OffsetArray([0 0 0 0 0 0 0;
                            0 0 1 4 7 0 0;
                            0 0 2 5 8 0 0;
                            0 0 3 6 9 0 0;
                            0 0 0 0 0 0 0], 0:4, -1:5)
    A32 = @inferred(PaddedView{Float32}(0.0, a, (0:4, -1:5)))
    @test eltype(A32) == Float32
    @test A32 == A
    @test A32[1, 1] === 1.0f0

    A = @inferred(PaddedView(0.0, a, (Base.OneTo(5), Base.OneTo(5)), (2:4, 2:4)))
    @test A == [0 0 0 0 0;
                0 1 4 7 0;
                0 2 5 8 0;
                0 3 6 9 0;
                0 0 0 0 0]
    A32 = @inferred(PaddedView{Float32}(0.0, a, (Base.OneTo(5), Base.OneTo(5)), (2:4, 2:4)))
    @test eltype(A32) == Float32
    @test A32 == A
    @test A32[2, 2] === 1.0f0

    @test A == @inferred(PaddedView(0.0, a, (5, 5), (2, 2)))
    A32 = @inferred(PaddedView{Float32}(0.0, a, (5, 5), (2, 2)))
    @test eltype(A32) == Float32
    @test A32 == A
    @test A32[2, 2] === 1.0f0

    B = @inferred(PaddedView(0.0, OffsetArray(a, (2:4, 2:4)), (Base.OneTo(5), Base.OneTo(5))))
    @test B == A

    @test PaddedView(0, ones(4), (3,), (1,)) == PaddedView(0, ones(4), (3,))

    let a = reshape(1:6, 2, 3) # [1 3 5; 2 4 6]
        @test @inferred(PaddedView(-1, a, (2, 2))) == [1 3
                                                       2 4]

        @test @inferred(PaddedView(-1, a, (1, 4)))  == [1 3 5 -1 ]

        @test @inferred(PaddedView(-1, a, (2, 4), (1, 1)))  == [1 3 5 -1
                                                                2 4 6 -1]

        @test @inferred(PaddedView(-1, a, (2, 4), (1, 2)))  == [-1 1 3 5
                                                                -1 2 4 6]

        @test @inferred(PaddedView(-1, a, (2, 3), (0, 2)))  == [-1  2  4
                                                                -1 -1 -1]

        @test @inferred(PaddedView(-1, a, (1, 4), (1, -1))) == [5 -1 -1 -1]

        # Create an OffsetArray with axes (0:1, 2:4)
        oa = OffsetArray(a, -1, 1)
        # Make a PaddedView of its elements (1:3, 1:2) (i.e pv[1, 1] == oa[1, 1] == "a[2, 0]" == -1)
        pv = @inferred PaddedView(-1, oa, (3, 2))
        @test axes(pv) === (Base.OneTo(3), Base.OneTo(2))
        @test pv[1:3, 1:2] == [-1 2; -1 -1; -1 -1]
        # Put oa[1, 1] in pv[2, 0] and make a PaddedView of the resulting (1:3, 1:2) elements
            # i.e. same as putting oa[0,2] in pv[1, 1]
        pv = @inferred PaddedView(-1, oa, (3, 2), (2, 0))
        @test pv[1:3, 1:2] == [1 3; 2 4; -1 -1]
    end

    # test for offset axes: the axes should be their own axes
    a = [1,2,3]
    P = PaddedView(-1, a, (0:3,))
    @test first(P) == P[first(LinearIndices(P))]
end

@testset "paddedviews" begin
    @testset "0-d array" begin
        a = reshape([1])
        pa, pb = paddedviews(-1, a, a)
        @test pa === pb
        @test pa == a

        pa = PaddedView{Int}(-1, a, (), ())
        @test pa == a
    end

    a1 = reshape([1,2], 2, 1)
    a2 = [1.0,2.0]'
    a1p, a2p = @inferred(paddedviews(0, a1, a2))
    @test a1p == [1 0; 2 0]
    @test a1p[CartesianIndices(a1)] == a1
    @test a2p == [1.0 2.0; 0.0 0.0]
    @test a2p[CartesianIndices(a2)] == a2
    @test eltype(a1p) == Int
    @test eltype(a2p) == Float64
    @test axes(a1p) === axes(a2p) === (Base.OneTo(2), Base.OneTo(2))
    a1p, a2p, a3p = @inferred(paddedviews(0, a1, a2, a1))
    @test a1p == a3p == [1 0; 2 0]
    @test a2p == [1.0 2.0; 0.0 0.0]
    a1p, a2p = PaddedViews.paddedviews_itr(0, [a1, a2])
    @test a1p == [1 0; 2 0]
    @test a2p == [1.0 2.0; 0.0 0.0]

    a3 = OffsetArray([1.0,2.0]', (0,-1))
    a1p, a3p = @inferred(paddedviews(0, a1, a3))
    @test a1p == OffsetArray([0 1; 0 2], 1:2, 0:1)
    @test a1p[CartesianIndices(a1)] == a1
    @test a3p == OffsetArray([1.0 2.0; 0.0 0.0], 1:2, 0:1)
    @test a3p[CartesianIndices(a3)] == a3
    @test eltype(a1p) == Int
    @test eltype(a3p) == Float64
    @test axes(a1p) === axes(a3p) === map(IdentityUnitRange, (1:2, 0:1))

    @test @inferred(paddedviews(3)) == ()
    @test_throws ErrorException PaddedViews.outerinds()
    # But a zero-dimensional input should not trigger that error
    a = reshape([5])
    @test paddedviews(-1, a) == (a,)


    # paddedviews
    a1 = reshape([1, 2, 3,4,5,6,7,8,9], 3, 3)
    a2 = [4 5;6 7]
    a2b = [4 5 -1; 6 7 -1; -1 -1 -1]
    a2c = [4 5; 6 7; -1 -1]
    a2r = [4 5 -1; 6 7 -1]

    a1f, a2f = paddedviews(-1, a1, a2;dims=(1))
    @test a2f == a2c
    a1f, a2f = paddedviews(-1, a1, a2;dims=(2))
    @test a2f == a2r
    a1f, a2f = paddedviews(-1, a1, a2;dims=(3))
    @test a2f == a2

    a1f, a2f = paddedviews(-1, a1, a2;dims=(1,2))
    @test a2f == a2b
    a1f, a2f = paddedviews(-1, a1, a2;dims=(2,3))
    @test a2f == a2r
    a1f, a2f = paddedviews(-1, a1, a2;dims=(1,2,3))
    @test a2f == a2b

    a1f, a2f = paddedviews(-1, a1, a2;dims=1)
    @test a2f == a2c
    a1f, a2f = paddedviews(-1, a1, a2;dims=2)
    @test a2f == a2r
    a1f, a2f = paddedviews(-1, a1, a2;dims=3)
    @test a2f == a2

    a1f, a2f = paddedviews(-1, a1, a2;dims=1:0)
    @test a2f == a2
    a1f, a2f = paddedviews(-1, a1, a2;dims=1:1)
    @test a2f == a2c
    a1f, a2f = paddedviews(-1, a1, a2;dims=1:2)
    @test a2f == a2b
    a1f, a2f = paddedviews(-1, a1, a2;dims=1:3)
    @test a2f == a2b

    a1f, a2f = paddedviews(-1, a1, a2;dims=2:3)
    @test a2f == a2r
    a1f, a2f = paddedviews(-1, a1, a2;dims=3:4)
    @test a2f == a2
end

@testset "sym_paddedviews" begin
    @testset "0-d array" begin
        a = reshape([1])
        pa, pb = sym_paddedviews(-1, a, a)
        @test pa === pb
        @test pa == a
    end

    # even case
    a1 = reshape([1,2], 2, 1)
    a2 = [1.0,2.0]'
    a1p, a2p = @inferred(sym_paddedviews(0, a1, a2))
    @test a1p == [1 0; 2 0]
    @test a1p[CartesianIndices(a1)] == a1
    @test a2p == [1.0 2.0; 0.0 0.0]
    @test a2p[CartesianIndices(a2)] == a2
    @test eltype(a1p) == Int
    @test eltype(a2p) == Float64
    @test axes(a1p) === axes(a2p) === map(IdentityUnitRange, (1:2, 1:2))

    a1 = reshape([1,2,3], 3, 1)
    a2 = [1.0,2.0,3.0]'
    a1p, a2p = @inferred(sym_paddedviews(0, a1, a2))

    @test @inferred(sym_paddedviews(0, a1)) == (a1,)

    @test a1p == OffsetArray([0 1 0; 0 2 0; 0 3 0], (1:3, 0:2))
    @test a1p[CartesianIndices(a1)] == a1
    @test a2p == OffsetArray([0.0 0.0 0.0; 1.0 2.0 3.0; 0.0 0.0 0.0], (0:2, 1:3))
    @test a2p[CartesianIndices(a2)] == a2
    @test eltype(a1p) == Int
    @test eltype(a2p) == Float64
    @test axes(a1p) === map(IdentityUnitRange, (1:3, 0:2))
    @test axes(a2p) === map(IdentityUnitRange, (0:2, 1:3))

    a1p, a2p, a3p = @inferred(sym_paddedviews(0, a1, a2, a1))
    @test a1p == a3p == OffsetArray([0 1 0; 0 2 0; 0 3 0], (1:3, 0:2))
    @test a2p == OffsetArray([0.0 0.0 0.0; 1.0 2.0 3.0; 0.0 0.0 0.0], (0:2, 1:3))

    a3 = OffsetArray([1.0,2.0,3.0]', (0,-1))
    a1p, a3p = @inferred(sym_paddedviews(0, a1, a3))

    @test a1p == OffsetArray([0 1 0; 0 2 0; 0 3 0], 1:3, 0:2)
    @test a1p[CartesianIndices(a1)] == a1
    @test a3p == OffsetArray([0.0 0.0 0.0; 1.0 2.0 3.0; 0.0 0.0 0.0], 0:2, 0:2)
    @test a3p[CartesianIndices(a3)] == a3
    @test eltype(a1p) == Int
    @test eltype(a3p) == Float64
    @test axes(a1p) == (1:3, 0:2)
    @test axes(a3p) == (0:2, 0:2)

    @test @inferred(sym_paddedviews(3)) == ()

    a1p, a2p = @inferred(PaddedViews.sym_paddedviews_itr(0, Matrix{Int}[a1, a2]))
    @test a1p == OffsetArray([0 1 0; 0 2 0; 0 3 0], (1:3, 0:2))
    @test a2p == OffsetArray([0 0 0; 1 2 3; 0 0 0], (0:2, 1:3))

    # sys_paddedviews cases
    a1 = reshape([1, 2, 3,4,5,6,7,8,9], 3, 3)
    a2 = reshape([5, 6], 2, 1)
    a2b = [-1 5 -1; -1 6 -1; -1 -1 -1]
    a2c = reshape([5, 6, -1], 3, 1)
    a2r =  [-1 5 -1; -1 6 -1]

    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=(1))
    @test a2f == a2c
    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=(2))
    @test collect(a2f) == a2r
    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=(3))
    @test a2f == a2

    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=(1,2))
    @test collect(a2f) == a2b
    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=(2,3))
    @test collect(a2f) == a2r
    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=(1,2,3))
    @test collect(a2f) == a2b

    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=1)
    @test collect(a2f) == a2c
    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=2)
    @test collect(a2f) == a2r
    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=3)
    @test collect(a2f)  == a2

    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=1:0)
    @test collect(a2f)  == a2
    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=1:1)
    @test collect(a2f)  == a2c
    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=1:2)
    @test collect(a2f)  == a2b
    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=1:3)
    @test collect(a2f)  == a2b

    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=2:3)
    @test collect(a2f)  == a2r
    a1f, a2f = sym_paddedviews(-1, a1, a2;dims=3:4)
    @test collect(a2f)  == a2
end

@testset "showarg" begin
    a = collect(reshape(1:9, 3, 3))
    pv = PaddedView(-1, a, (0:4, 1:3))
    io = IOBuffer()
    show(IOContext(io, :displaysize=>(1000,1000)), MIME("text/plain"), pv)
    str = String(take!(io))
    @test endswith(str, "PaddedView(-1, ::$(Array{Int,2}), (0:4, 1:3)) with eltype $Int with indices 0:4×1:3:\n -1  -1  -1\n  1   4   7\n  2   5   8\n  3   6   9\n -1  -1  -1")
end

@testset "similar" begin
    for (T, v) in ((Int, 0),
                   (Int, 0.),
                   (Float64, 0.),
                   (Missing, missing),
                   (Nothing, nothing))
        A = reshape(1:9, 3, 3)
        Ap = @inferred(PaddedView(v, A, (0:4, 0:4)))
        B = similar(Ap)
        @test eltype(B) == eltype(Ap)
        @test size(B) == size(Ap)
        @test axes(B) == axes(Ap)

        B = similar(Ap, (4, 4))
        @test eltype(B) == eltype(Ap)
        @test size(B) == (4, 4)
        @test axes(B) == (Base.OneTo(4), Base.OneTo(4))

        B = similar(Ap, (4, 4))
        @test eltype(B) == eltype(Ap)
        @test size(B) == (4, 4)
        @test axes(B) == (Base.OneTo(4), Base.OneTo(4))

        B = similar(Ap, Int, (4, 4))
        @test eltype(B) == Int
        @test size(B) == (4, 4)
        @test axes(B) == (Base.OneTo(4), Base.OneTo(4))
    end
end

@testset "nothing/missing" begin
    for (FT, T) in ((Missing, Float32),
                    (Nothing, Float32),
                    (Float32, Nothing),
                    (Float32, Missing),
                    (Nothing, Missing),
                    (Missing, Nothing))
        @test @inferred(filltype(FT, T)) === Union{FT, T}
        @test @inferred(filltype(T, Union{FT, T})) === Union{FT, T}
        @test @inferred(filltype(FT, Union{FT, T})) === Union{FT, T}
    end

    for (T, v) in ((Missing, missing),
                 (Nothing, nothing))
        A = reshape(1:9, 3, 3)
        Ap = @inferred(PaddedView(v, A, (0:4, 0:4)))
        @test axes(Ap) == (OffsetArrays.IdOffsetRange(0:4), OffsetArrays.IdOffsetRange(0:4))
        @test eltype(Ap) === Union{T, eltype(A)}
        @test Ap[0, 0] === v
        @test Ap[1, 1] === 1
        @test Ap[axes(A)...] == A

        Ap = @inferred(PaddedView(v, A, (5, 5)))
        @test axes(Ap) === (Base.OneTo(5), Base.OneTo(5))
        @test eltype(Ap) === Union{T, eltype(A)}
        @test Ap[4, 4] === v
        @test Ap[1, 1] === 1
        @test Ap[axes(A)...] == A
    end
end
