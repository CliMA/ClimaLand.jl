using StaticArrayInterface
using LinearAlgebra
using StaticArrays
using Static
using Test

x = @SVector [1,2,3]
@test @inferred(StaticArrayInterface.device(typeof(x))) === StaticArrayInterface.CPUTuple()

x = @MVector [1,2,3]
@test @inferred(StaticArrayInterface.device(typeof(x))) === StaticArrayInterface.CPUPointer()

@test isone(StaticArrayInterface.known_first(typeof(StaticArrays.SOneTo(7))))
@test StaticArrayInterface.known_last(typeof(StaticArrays.SOneTo(7))) == 7
@test StaticArrayInterface.known_length(typeof(StaticArrays.SOneTo(7))) == 7

@test StaticArrayInterface.parent_type(SizedVector{1, Int, Vector{Int}}) <: Vector{Int}
@test StaticArrayInterface.known_length(@inferred(StaticArrayInterface.indices(SOneTo(7)))) == 7

x = view(SArray{Tuple{3,3,3}}(ones(3,3,3)), :, SOneTo(2), 2)
@test @inferred(StaticArrayInterface.known_length(x)) == 6
@test @inferred(StaticArrayInterface.known_length(x')) == 6

v = @SVector rand(8);
A = @MMatrix rand(7, 6);
T = SizedArray{Tuple{5,4,3}}(zeros(5,4,3));
@test @inferred(StaticArrayInterface.static_length(v)) === StaticInt(8)
@test @inferred(StaticArrayInterface.static_length(A)) === StaticInt(42)
@test @inferred(StaticArrayInterface.static_length(T)) === StaticInt(60)

Am = @MMatrix rand(2,10);
@test @inferred(StaticArrayInterface.static_strides(view(Am,1,:))) === (StaticInt(2),)

@test @inferred(StaticArrayInterface.contiguous_axis(@SArray(zeros(2,2,2)))) === Static.StaticInt(1)
@test @inferred(StaticArrayInterface.contiguous_axis_indicator(@SArray(zeros(2,2,2)))) == (true,false,false)
@test @inferred(StaticArrayInterface.contiguous_batch_size(@SArray(zeros(2,2,2)))) === Static.StaticInt(0)
@test @inferred(StaticArrayInterface.stride_rank(@SArray(zeros(2,2,2)))) == (1, 2, 3)
@test @inferred(StaticArrayInterface.is_column_major(@SArray(zeros(2,2,2)))) === True()
@test @inferred(StaticArrayInterface.dense_dims(@SArray(zeros(2,2,2)))) == (true,true,true)

