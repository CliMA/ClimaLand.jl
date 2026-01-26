using LayoutPointers, StaticArrayInterface, Aqua, Test
struct SizedWrapper{M,N,T,AT<:AbstractMatrix{T}} <: AbstractMatrix{T}
  A::AT
end
StaticArrayInterface.is_forwarding_wrapper(::Type{<:SizedWrapper}) = true
SizedWrapper{M,N}(A::AT) where {M,N,T,AT<:AbstractMatrix{T}} = SizedWrapper{M,N,T,AT}(A)
Base.size(::SizedWrapper{M,N}) where {M,N} = (M, N)
Base.getindex(A::SizedWrapper, i...) = getindex(parent(A), i...)
Base.parent(dw::SizedWrapper) = dw.A
StaticArrayInterface.parent_type(::Type{SizedWrapper{M,N,T,AT}}) where {M,N,T,AT} = AT
LayoutPointers.memory_reference(dw::SizedWrapper) =
  LayoutPointers.memory_reference(parent(dw))
StaticArrayInterface.contiguous_axis(dw::SizedWrapper) =
  LayoutPointers.contiguous_axis(parent(dw))
StaticArrayInterface.contiguous_batch_size(dw::SizedWrapper) =
  LayoutPointers.contiguous_batch_size(parent(dw))
LayoutPointers.val_stride_rank(dw::SizedWrapper) =
  LayoutPointers.val_stride_rank(parent(dw))
function StaticArrayInterface.static_strides(dw::SizedWrapper{M,N,T}) where {M,N,T}
  if LayoutPointers.val_stride_rank(dw) === Val((1, 2))
    return LayoutPointers.One(), LayoutPointers.StaticInt{M}()
  else#if LayoutPointers.val_stride_rank(dw) === Val((2,1))
    return LayoutPointers.StaticInt{N}(), LayoutPointers.One()
  end
end
StaticArrayInterface.offsets(dw::SizedWrapper) = LayoutPointers.offsets(parent(dw))
LayoutPointers.val_dense_dims(dw::SizedWrapper{T,N}) where {T,N} =
  LayoutPointers.val_dense_dims(parent(dw))

@testset "LayoutPointers.jl" begin
  Aqua.test_all(LayoutPointers)

  println("Grouped Strided Pointers")
  @time @testset "Grouped Strided Pointers" begin
    M, K, N = 4, 5, 6
    A = Matrix{Float64}(undef, M, K)
    B = Matrix{Float64}(undef, K, N)
    C = Matrix{Float64}(undef, M, N)

    GC.@preserve A B C begin
      fs = (false, true)#[identity, adjoint]
      for ai ∈ fs, bi ∈ fs, ci ∈ fs
        At = ai ? A : (similar(A')')
        Bt = bi ? B : (similar(B')')
        Ct = ci ? C : (similar(C')')
        spdw = LayoutPointers.DensePointerWrapper{(true, true)}(
          LayoutPointers.stridedpointer(At),
        )
        gsp, pres = @inferred(
          LayoutPointers.grouped_strided_pointer(
            (spdw, Bt, Ct),
            Val{(((1, 1), (3, 1)), ((1, 2), (2, 1)), ((2, 2), (3, 2)))}(),
          )
        )
        if ai === ci
          @test sizeof(gsp.strides) == 2sizeof(Int)
        end
        # Test to confirm that redundant strides are not stored in the grouped strided pointer
        @test sizeof(gsp) == sizeof(Int) * (6 - (ai & ci) - ((!ai) & bi) - ((!bi) & (!ci)))
        @test sizeof(gsp.offsets) == 0
        pA, pB, pC = @inferred(LayoutPointers.stridedpointers(gsp))
        @test pA === stridedpointer(At)
        @test pB === stridedpointer(Bt)
        @test pC === stridedpointer(Ct)
        Btsw = SizedWrapper{K,N}(Bt)
        gsp2, pres2 = @inferred(
          LayoutPointers.grouped_strided_pointer(
            (At, Btsw, Ct),
            Val{(((1, 1), (3, 1)), ((1, 2), (2, 1)), ((2, 2), (3, 2)))}(),
          )
        )
        @test sizeof(gsp2) == sizeof(Int) * (5 - (ai & ci) - ((!ai) & bi) - ((!bi) & (!ci)))

        pA2, pB2, pC2 = @inferred(LayoutPointers.stridedpointers(gsp2))
        @test pointer(pA2) == pointer(At)
        @test pointer(pB2) == pointer(Bt)
        @test pointer(pC2) == pointer(Ct)
        @test strides(pA2) == strides(pA)
        @test strides(pB2) == strides(pB)
        @test strides(pC2) == strides(pC)
      end
    end

    data_in_large = Array{Float64}(undef, 4, 4, 4, 4, 1)
    data_in = view(data_in_large, :, 1, :, :, 1)
    tmp1 = Array{Float64}(undef, 4, 4, 4)
    sp_data_in, sp_tmp1 = LayoutPointers.stridedpointers(
      LayoutPointers.grouped_strided_pointer((data_in, tmp1), Val((((1, 1), (2, 1)),)))[1],
    )
    @test sp_data_in === stridedpointer(data_in)
    @test sp_tmp1 === stridedpointer(tmp1)
  end

end
