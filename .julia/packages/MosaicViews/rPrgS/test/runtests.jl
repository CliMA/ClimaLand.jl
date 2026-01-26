using MosaicViews
using Test
using ImageCore, ColorVectorSpace
using OffsetArrays

# Because of the difference in axes types between paddedviews (Base.OneTo) and
# sym_paddedviews (UnitRange), the return type of `_padded_cat` isn't inferrable to a
# concrete type. But it is inferrable to a Union{A,B} where both A and B are concrete.
# While `@inferred(_padded_cat((A, B), ...))` would therefore fail, this is a close substitute.
function _checkinferred_paddedcat(V, A; kwargs...)
    vd = MosaicViews.valdim(first(A))
    RTs = Base.return_types(MosaicViews._padded_cat, (typeof(A), Bool, eltype(V), typeof(vd)))
    @test length(RTs) == 1
    RT = RTs[1]
    @test isconcretetype(RT) || (isa(RT, Union) && isconcretetype(RT.a) && isconcretetype(RT.b))
    return V
end
function checkinferred_mosaic(As::Tuple; kwargs...)
    V = mosaic(As; kwargs...)
    return _checkinferred_paddedcat(V, As)
end
function checkinferred_mosaic(A...; kwargs...)
    V = mosaic(A...; kwargs...)
    return _checkinferred_paddedcat(V, A)
end
function checkinferred_mosaic(As::AbstractVector{<:AbstractArray}; kwargs...)
    V = mosaic(As; kwargs...)
    return _checkinferred_paddedcat(V, (As...,); kwargs...)
end

@testset "MosaicView" begin
    @test_throws ArgumentError MosaicView(rand(2,2,2,2,2))

    @testset "1D input" begin
        A = [1,2,3]
        mva = @inferred MosaicView(A)
        @test parent(mva) === A
        @test typeof(mva) <: MosaicView
        @test eltype(mva) == eltype(A)
        @test size(mva) == (3, 1) # 1D vector is lifted to 2D matrix
        @test axes(mva) == (Base.OneTo(3), Base.OneTo(1))
        @test @inferred(getindex(mva, 1, 1)) === 1
        @test @inferred(getindex(mva, 2, 1)) === 2
    end

    @testset "2D input" begin
        A = [1 2;3 4]
        mva = @inferred MosaicView(A)
        @test parent(mva) === A
        @test typeof(mva) <: MosaicView
        @test eltype(mva) == eltype(A)
        @test size(mva) == (2, 2)
        @test axes(mva) == (Base.OneTo(2), Base.OneTo(2))
        @test @inferred(getindex(mva, 1, 1)) === 1
        @test @inferred(getindex(mva, 2, 1)) === 3
    end

    @testset "3D input" begin
        A = zeros(Int,2,2,2)
        A[:,:,1] = [1 2; 3 4]
        A[:,:,2] = [5 6; 7 8]
        mva = @inferred MosaicView(A)
        @test parent(mva) === A
        @test typeof(mva) <: MosaicView
        @test eltype(mva) == eltype(A)
        @test size(mva) == (4, 2)
        @test @inferred(getindex(mva,1,1)) === 1
        @test @inferred(getindex(mva,2,1)) === 3
        @test_throws BoundsError mva[0,1]
        @test_throws BoundsError mva[1,0]
        @test_throws BoundsError mva[1,3]
        @test_throws BoundsError mva[5,1]
        @test all(mva .== vcat(A[:,:,1],A[:,:,2]))
        # singleton dimension doesn't change anything
        @test mva == MosaicView(reshape(A,2,2,2,1))
    end

    @testset "4D input" begin
        A = zeros(Int,2,2,1,2)
        A[:,:,1,1] = [1 2; 3 4]
        A[:,:,1,2] = [5 6; 7 8]
        mva = @inferred MosaicView(A)
        @test parent(mva) === A
        @test typeof(mva) <: MosaicView
        @test eltype(mva) == eltype(A)
        @test size(mva) == (2, 4)
        @test @inferred(getindex(mva,1,1)) === 1
        @test @inferred(getindex(mva,2,1)) === 3
        @test_throws BoundsError mva[0,1]
        @test_throws BoundsError mva[1,0]
        @test_throws BoundsError mva[3,1]
        @test_throws BoundsError mva[1,5]
        @test all(mva .== hcat(A[:,:,1,1],A[:,:,1,2]))
        A = zeros(Int,2,2,2,3)
        A[:,:,1,1] = [1 2; 3 4]
        A[:,:,1,2] = [5 6; 7 8]
        A[:,:,1,3] = [9 10; 11 12]
        A[:,:,2,1] = [13 14; 15 16]
        A[:,:,2,2] = [17 18; 19 20]
        A[:,:,2,3] = [21 22; 23 24]
        mva = @inferred MosaicView(A)
        @test parent(mva) === A
        @test typeof(mva) <: MosaicView
        @test eltype(mva) == eltype(A)
        @test size(mva) == (4, 6)
        @test @inferred(getindex(mva,1,1)) === 1
        @test @inferred(getindex(mva,2,1)) === 3
        @test all(mva .== vcat(hcat(A[:,:,1,1],A[:,:,1,2],A[:,:,1,3]), hcat(A[:,:,2,1],A[:,:,2,2],A[:,:,2,3])))
    end
end

@testset "mosaicview" begin
    @testset "1D input" begin
        A1 = [1, 2, 3]
        A2 = [4, 5, 6]
        mva = checkinferred_mosaic(A1, A2)
        @test typeof(mva) <: MosaicView
        @test eltype(mva) == eltype(A1)
        @test size(mva) == (6, 1)
        @test @inferred(getindex(mva, 1, 1)) == 1

        @test checkinferred_mosaic(A1, A2; nrow=1) == [1 4; 2 5; 3 6]
        @test mosaic(A1, A2; nrow=2) == reshape([1, 2, 3, 4, 5, 6], (6, 1))
        @test mosaic(A1, A2) == mosaic([A1, A2]) == mosaic((A1, A2))
    end

    @testset "2D input" begin
        A = [i*ones(Int, 2, 3) for i in 1:4]
        Ao = [i*ones(Int, 0i:1, 0:2) for i in 1:4]

        for B in (A, tuple(A...), Ao)
            @test_throws ArgumentError mosaic(B, nrow=0)
            @test_throws ArgumentError mosaic(B, ncol=0)
            @test_throws ArgumentError mosaic(B, nrow=1, ncol=1)

            mva = checkinferred_mosaic(B)
            @test mosaic(B...) == mva
            @test typeof(mva) <: MosaicView
            @test eltype(mva) == eltype(eltype(B))
            @test size(mva) == (8, 3)
            @test @inferred(getindex(mva,3,1)) === 2
            @test collect(mva) == [
                1  1  1
                1  1  1
                2  2  2
                2  2  2
                3  3  3
                3  3  3
                4  4  4
                4  4  4
            ]

            mva = checkinferred_mosaic(B, nrow=2)
            @test typeof(mva) <: MosaicView
            @test eltype(mva) == eltype(eltype(B))
            @test size(mva) == (4, 6)
            @test collect(mva) == [
             1  1  1  3  3  3
             1  1  1  3  3  3
             2  2  2  4  4  4
             2  2  2  4  4  4
            ]
        end

        @test mosaic(A...) == mosaic(A)
        @test mosaic(A..., nrow=2) == mosaic(A, nrow=2)
        @test mosaic(A..., nrow=2, rowmajor=true) == mosaic(A, nrow=2, rowmajor=true)

        A1 = reshape([1 2 3], (1, 3))
        A2 = reshape([4;5;6], (3, 1))
        @test checkinferred_mosaic([A1, A2]; center=false) == [
         1 2 3;
         0 0 0;
         0 0 0;
         4 0 0;
         5 0 0;
         6 0 0
        ]
        @test mosaic([A1, A2]; center=true) == [
         0 0 0;
         1 2 3;
         0 0 0;
         0 4 0;
         0 5 0;
         0 6 0
        ]
        @test mosaic([A1, A2]) == mosaic([A1, A2]; center=true)

        # same size but different axes
        A1 = fill(1, 1:2, 1:2)
        A2 = fill(2, 2:3, 2:3)
        @test collect(checkinferred_mosaic(A1, A2; center=true)) == [
            1 1 0;
            1 1 0;
            0 0 0;
            2 2 0;
            2 2 0;
            0 0 0;
        ]
        @test collect(mosaic(A1, A2; center=false)) == [
            1 1 0;
            1 1 0;
            0 0 0;
            0 0 0;
            0 2 2;
            0 2 2;
        ]

        @test mosaicview(A1) == A1 # a trivial case
    end

    @testset "3D input" begin
        A = [(k+1)*l-1 for i in 1:2, j in 1:3, k in 1:2, l in 1:2]
        B = reshape(A, 2, 3, :)
        @test_throws ArgumentError mosaicview(B, nrow=0)
        @test_throws ArgumentError mosaicview(B, ncol=0)
        @test_throws ArgumentError mosaicview(B, nrow=1, ncol=1)

        mva = checkinferred_mosaic(B)
        @test typeof(mva) <: MosaicView
        @test eltype(mva) == eltype(B)
        @test size(mva) == (8, 3)
        @test @inferred(getindex(mva,3,1)) === 2
        @test mva == [
            1  1  1
            1  1  1
            2  2  2
            2  2  2
            3  3  3
            3  3  3
            5  5  5
            5  5  5
        ]
        mva = checkinferred_mosaic(B, nrow=2)
        @test mva == MosaicView(A)
        @test typeof(mva) != typeof(MosaicView(A))
        @test parent(parent(mva)).data == B
        @test typeof(mva) <: MosaicView
        @test eltype(mva) == eltype(B)
        @test size(mva) == (4, 6)

        @test mosaic(B, B) == mosaicview(cat(B, B; dims=4))
        @test mosaic(B, B, nrow=2) == mosaicview(cat(B, B; dims=4), nrow=2)
        @test mosaic(B, B, nrow=2, rowmajor=true) == mosaicview(cat(B, B; dims=4), nrow=2, rowmajor=true)
    end

    @testset "4D input" begin
        A = [(k+1)*l-1 for i in 1:2, j in 1:3, k in 1:2, l in 1:2]
        @test_throws ArgumentError mosaicview(A, nrow=0)
        @test_throws ArgumentError mosaicview(A, ncol=0)
        @test_throws ArgumentError mosaicview(A, nrow=1, ncol=1)

        mva = mosaicview(A)
        @test mva == MosaicView(A)
        @test typeof(mva) != typeof(MosaicView(A))
        @test parent(parent(mva)).data == reshape(A, 2, 3, :)
        @test typeof(mva) <: MosaicView
        @test eltype(mva) == eltype(A)
        @test size(mva) == (4, 6)
        mva = mosaicview(A, npad=1)
        @test typeof(mva) <: MosaicView
        @test eltype(mva) == eltype(A)
        @test size(mva) == (5, 7)
        @test mva == mosaicview(A, nrow=2, npad=1)
        @test mva == mosaicview(A, ncol=2, npad=1)
        @test @inferred(getindex(mva,3,1)) === 0
        @test @inferred(getindex(mva,2,5)) === 3
        @test mva == [
            1  1  1  0  3  3  3
            1  1  1  0  3  3  3
            0  0  0  0  0  0  0
            2  2  2  0  5  5  5
            2  2  2  0  5  5  5
        ]
        mva = mosaicview(A, ncol=3, npad=1)
        @test typeof(mva) <: MosaicView
        @test eltype(mva) == eltype(A)
        @test size(mva) == (5, 11)
        @test mva == [
            1  1  1  0  3  3  3  0  0  0  0
            1  1  1  0  3  3  3  0  0  0  0
            0  0  0  0  0  0  0  0  0  0  0
            2  2  2  0  5  5  5  0  0  0  0
            2  2  2  0  5  5  5  0  0  0  0
        ]
        mva = mosaicview(A, rowmajor=true, ncol=3, npad=1)
        @test typeof(mva) <: MosaicView
        @test eltype(mva) == eltype(A)
        @test size(mva) == (5, 11)
        @test @inferred(getindex(mva,3,1)) === 0
        @test @inferred(getindex(mva,2,5)) === 2
        @test mva == [
            1  1  1  0  2  2  2  0  3  3  3
            1  1  1  0  2  2  2  0  3  3  3
            0  0  0  0  0  0  0  0  0  0  0
            5  5  5  0  0  0  0  0  0  0  0
            5  5  5  0  0  0  0  0  0  0  0
        ]
        mva = mosaicview(A, rowmajor=true, ncol=3, npad=2)
        @test typeof(mva) <: MosaicView
        @test eltype(mva) == eltype(A)
        @test size(mva) == (6, 13)
        @test mva == [
            1  1  1  0  0  2  2  2  0  0  3  3  3
            1  1  1  0  0  2  2  2  0  0  3  3  3
            0  0  0  0  0  0  0  0  0  0  0  0  0
            0  0  0  0  0  0  0  0  0  0  0  0  0
            5  5  5  0  0  0  0  0  0  0  0  0  0
            5  5  5  0  0  0  0  0  0  0  0  0  0
        ]
        mva = mosaicview(A, fillvalue=-1.0, rowmajor=true, ncol=3, npad=1)
        @test typeof(mva) <: MosaicView
        @test eltype(mva) == eltype(A)
        @test size(mva) == (5, 11)
        @test mva == [
             1   1   1  -1   2   2   2  -1   3   3   3
             1   1   1  -1   2   2   2  -1   3   3   3
            -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1
             5   5   5  -1  -1  -1  -1  -1  -1  -1  -1
             5   5   5  -1  -1  -1  -1  -1  -1  -1  -1
        ]

        @test mosaic(A, A) == mosaicview(cat(A, A; dims=5))
        @test mosaic(A, A, nrow=2) == mosaicview(cat(A, A; dims=4), nrow=2)
        @test mosaic(A, A, nrow=2, rowmajor=true) == mosaicview(cat(A, A; dims=4), nrow=2, rowmajor=true)
    end

    @testset "Colorant Array" begin
        A = rand(RGB{Float32}, 2, 3, 2, 2)
        mvaa = mosaicview(A)
        @test eltype(mvaa) == eltype(A)
        @test mvaa == @inferred(MosaicView(A))
        mvaa = mosaicview(A, rowmajor=true, ncol=3)
        @test eltype(mvaa) == eltype(A)
        @test @inferred(getindex(mvaa, 3, 4)) == RGB(0,0,0)
        mvaa = mosaicview(A, fillvalue=colorant"white", rowmajor=true, ncol=3)
        @test eltype(mvaa) == eltype(A)
        @test @inferred(getindex(mvaa, 3, 4)) == RGB(1,1,1)

        # this should work regardless they're of different size and color
        @test_nowarn mosaic(rand(RGB{Float32}, 4, 4),
                            rand(Gray{N0f8}, 5, 5))
    end

    # all arrays should have the same dimension
    @test_throws ArgumentError mosaic(ones(2), ones(1, 2))
    @test_throws ArgumentError mosaic((ones(2), ones(1, 2)))
    @test_throws ArgumentError mosaic([ones(2), ones(1, 2)])

    @testset "filltype" begin
        # always a concrete type
        A = checkinferred_mosaic(rand(N0f8, 4, 4), rand(Float64, 4, 4), rand(Float32, 4, 4))
        @test eltype(A) == Float64

        A = mosaic(Any[1 2 3; 4 5 6], rand(Float32, 4, 4))
        @test eltype(A) == Float32
        A = mosaic(rand(Float32, 4, 4), Any[1 2 3; 4 5 6])
        @test eltype(A) == Float32

        # FIXME:
        # A = checkinferred_mosaic(rand(Float64, 4, 4), Union{Missing, Float32}[1 2 3; 4 5 6])
        A = mosaic(rand(Float64, 4, 4), Union{Missing, Float32}[1 2 3; 4 5 6])
        @test eltype(A) == Union{Missing, Float64}
    end
end


@testset "deprecations" begin
    @info "deprecations are expected"
    # mosaicview -> mosaic deprecations
    A = [(k+1)*l-1 for i in 1:2, j in 1:3, k in 1:2, l in 1:2]
    mva_old = mosaicview(A, A, rowmajor=true, ncol=3, npad=1)
    mva_new = mosaic(A, A, rowmajor=true, ncol=3, npad=1)
    @test mva_old == mva_new
    mva_old = mosaicview(A, A, A, rowmajor=true, ncol=3, npad=1)
    mva_new = mosaic(A, A, A, rowmajor=true, ncol=3, npad=1)
    @test mva_old == mva_new
    mva_old = mosaicview([A, A], rowmajor=true, ncol=3, npad=1)
    mva_new = mosaic([A, A], rowmajor=true, ncol=3, npad=1)
    @test mva_old == mva_new
    mva_old = mosaicview((A, A), rowmajor=true, ncol=3, npad=1)
    mva_new = mosaic((A, A), rowmajor=true, ncol=3, npad=1)
    @test mva_old == mva_new
    # center keyword for `mosaicview` (still applies for `mosaic`)
    A = [(k+1)*l-1 for i in 1:2, j in 1:3, k in 1:2, l in 1:2]
    B = reshape(A, 2, 3, :)
    @test mosaicview(B, center=2) == mosaicview(B) # no op
    @test mosaicview(A, center=2) == mosaicview(A) # no op
end
