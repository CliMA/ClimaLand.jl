A = zeros(3, 4, 5);
A[:] = 1:60
Ap = @view(PermutedDimsArray(A,(3,1,2))[:,1:2,1])';

ap_index = StaticArrayInterface.StrideIndex(Ap)
for x_i in axes(Ap, 1)
    for y_i in axes(Ap, 2)
        @test ap_index[x_i, y_i] == ap_index[x_i, y_i]
    end
end
@test @inferred(StaticArrayInterface.known_offsets(ap_index)) === StaticArrayInterface.known_offsets(Ap)
@test @inferred(StaticArrayInterface.known_offset1(ap_index)) === StaticArrayInterface.known_offset1(Ap)
@test @inferred(StaticArrayInterface.offsets(ap_index, 1)) === StaticArrayInterface.offset1(Ap)
@test @inferred(StaticArrayInterface.offsets(ap_index, static(1))) === StaticArrayInterface.offset1(Ap)
@test @inferred(StaticArrayInterface.known_strides(ap_index)) === StaticArrayInterface.known_strides(Ap)
@test @inferred(StaticArrayInterface.contiguous_axis(ap_index)) == 1
@test @inferred(StaticArrayInterface.contiguous_axis(StaticArrayInterface.StrideIndex{2,(1,2),nothing,NTuple{2,Int},NTuple{2,Int}})) === nothing
@test @inferred(StaticArrayInterface.stride_rank(ap_index)) == (1, 3)

let v = Float64.(1:10)', v2 = transpose(parent(v))
  sv = @view(v[1:5])'
  sv2 = @view(v2[1:5])'
  @test @inferred(StaticArrayInterface.StrideIndex(sv)) === @inferred(StaticArrayInterface.StrideIndex(sv2)) === StaticArrayInterface.StrideIndex{2, (2, 1), 2}((StaticInt(1), StaticInt(1)), (StaticInt(1), StaticInt(1)))
  @test @inferred(StaticArrayInterface.stride_rank(parent(sv))) === @inferred(StaticArrayInterface.stride_rank(parent(sv2))) === (StaticInt(1),)
end


@testset "IndicesInfo" begin

  struct SplatFirst end

  import StaticArrayInterface: IndicesInfo

  StaticArrayInterface.is_splat_index(::Type{SplatFirst}) = true

  @test @inferred(IndicesInfo(SubArray{Float64, 2, Vector{Float64}, Tuple{Base.ReshapedArray{Int64, 2, UnitRange{Int64}, Tuple{}}}, true})) isa
      IndicesInfo{1,(1,),((1,2),)}

  @test @inferred(IndicesInfo{1}((Tuple{Vector{Int}}))) isa IndicesInfo{1, (1,), (1,)}

  @test @inferred(IndicesInfo{2}(Tuple{Vector{Int}})) isa IndicesInfo{2, (:,), (1,)}

  @test @inferred(IndicesInfo{1}(Tuple{SplatFirst})) isa IndicesInfo{1, (1,), (1,)}

  @test @inferred(IndicesInfo{2}(Tuple{SplatFirst})) isa IndicesInfo{2, ((1,2),), ((1, 2),)}

  @test @inferred(IndicesInfo{5}(typeof((:,[CartesianIndex(1,1),CartesianIndex(1,1)], 1, ones(Int, 2, 2), :, 1)))) isa
      IndicesInfo{5, (1, (2, 3), 4, 5, 0, 0), (1, 2, 0, (3, 4), 5, 0)}

  @test @inferred(IndicesInfo{10}(Tuple{Vararg{Int,10}})) isa
      IndicesInfo{10, (1, 2, 3, 4, 5, 6, 7, 8, 9, 10), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)}

  @test @inferred(IndicesInfo{10}(typeof((1, CartesianIndex(2, 1), 2, CartesianIndex(1, 2), 1, CartesianIndex(2, 1), 2)))) isa
      IndicesInfo{10, (1, (2, 3), 4, (5, 6), 7, (8, 9), 10), (0, 0, 0, 0, 0, 0, 0)}

  @test @inferred(IndicesInfo{10}(typeof((fill(true, 4, 4), 2, fill(true, 4, 4), 2, 1, fill(true, 4, 4), 1)))) isa
      IndicesInfo{10, ((1, 2), 3, (4, 5), 6, 7, (8, 9), 10), (1, 0, 2, 0, 0, 3, 0)}

  @test @inferred(IndicesInfo{10}(typeof((1, SplatFirst(), 2, SplatFirst(), CartesianIndex(1, 1))))) isa
      IndicesInfo{10, (1, (2, 3, 4, 5, 6), 7, 8, (9, 10)), (0, (1, 2, 3, 4, 5), 0, 6, 0)}
end
