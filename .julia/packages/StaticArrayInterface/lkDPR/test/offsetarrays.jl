using StaticArrayInterface
using OffsetArrays
using Static
using Test

A = zeros(3, 4, 5);
O = OffsetArray(A, 3, 7, 10);
Op = PermutedDimsArray(O,(3,1,2));
@test @inferred(StaticArrayInterface.offsets(O)) === (4, 8, 11)
@test @inferred(StaticArrayInterface.offsets(Op)) === (11, 4, 8)

@test @inferred(StaticArrayInterface.static_to_indices(O, (:, :, :))) == (4:6, 8:11, 11:15)
@test @inferred(StaticArrayInterface.static_to_indices(Op, (:, :, :))) == (11:15, 4:6, 8:11)

@test @inferred(StaticArrayInterface.offsets((1,2,3))) === (StaticInt(1),)
o = OffsetArray(vec(A), 8);
@test @inferred(StaticArrayInterface.offset1(o)) === 9

@test @inferred(StaticArrayInterface.device(OffsetArray(view(PermutedDimsArray(A, (3,1,2)), 1, :, 2:4)', 3, -173))) === StaticArrayInterface.CPUPointer()
@test @inferred(StaticArrayInterface.device(view(OffsetArray(A,2,3,-12), 4, :, -11:-9))) === StaticArrayInterface.CPUPointer()
@test @inferred(StaticArrayInterface.device(view(OffsetArray(A,2,3,-12), 3, :, [-11,-10,-9])')) === StaticArrayInterface.CPUIndex()

@test @inferred(StaticArrayInterface.indices(OffsetArray(view(PermutedDimsArray(A, (3,1,2)), 1, :, 2:4)', 3, -173),1)) === Base.Slice(Static.OptionallyStaticUnitRange(4,6))
@test @inferred(StaticArrayInterface.indices(OffsetArray(view(PermutedDimsArray(A, (3,1,2)), 1, :, 2:4)', 3, -173),2)) === Base.Slice(Static.OptionallyStaticUnitRange(-172,-170))

@test @inferred(StaticArrayInterface.device(OffsetArray(1:10))) === StaticArrayInterface.CPUIndex()
@test @inferred(StaticArrayInterface.device(OffsetArray(@view(reshape(1:8, 2,2,2)[1,1:2,:]),-3,4))) === StaticArrayInterface.CPUIndex()
@test @inferred(StaticArrayInterface.device(OffsetArray(zeros(2,2,2),8,-2,-5))) === StaticArrayInterface.CPUPointer()

offset_view = @view OffsetArrays.centered(zeros(eltype(A), 5, 5))[:, begin]; # SubArray of OffsetArray
@test @inferred(StaticArrayInterface.offsets(offset_view)) == (-2,)

B = OffsetArray(PermutedDimsArray(rand(2,3,4), (2,3,1)));
@test @inferred(StaticArrayInterface.StrideIndex(B)) === StaticArrayInterface.StrideIndex{3, (2, 3, 1), 3}((2, 6, static(1)), (1, 1, 1))
