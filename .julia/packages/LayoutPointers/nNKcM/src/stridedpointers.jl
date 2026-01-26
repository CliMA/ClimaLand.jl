@inline function mulsizeof(::Type{T}, x::Number) where {T}
  (Base.allocatedinline(T) ? sizeof(T) : sizeof(Int)) * x
end
@generated function mulsizeof(::Type{T}, ::StaticInt{N}) where {T,N}
  st = Base.allocatedinline(T) ? sizeof(T) : sizeof(Int)
  Expr(:block, Expr(:meta, :inline), Expr(:call, Expr(:curly, :StaticInt, N * st)))
end
@inline mulsizeof(::Type{T}, ::Tuple{}) where {T} = ()
@inline mulsizeof(::Type{T}, x::Tuple{X}) where {T,X} =
  (mulsizeof(T, getfield(x, 1, false)),)
@inline mulsizeof(::Type{T}, x::Tuple{X1,X2,Vararg}) where {T,X1,X2} =
  (mulsizeof(T, getfield(x, 1, false)), mulsizeof(T, Base.tail(x))...)

@inline bytestrides(A::AbstractArray{T}) where {T} =
  mulsizeof(T, StaticArrayInterface.static_strides(A))

@inline memory_reference(A::NTuple) = memory_reference(StaticArrayInterface.device(A), A)
@inline memory_reference(A::AbstractArray) =
  memory_reference(StaticArrayInterface.device(A), A)
@inline memory_reference(A::BitArray) =
  Base.unsafe_convert(Ptr{Bit}, pointer(A.chunks)), A.chunks
@inline memory_reference(::CPUPointer, A) = pointer(A), preserve_buffer(A)
@inline memory_reference(
  ::CPUPointer,
  A::Union{
    LinearAlgebra.Adjoint,
    Base.ReshapedArray,
    Base.PermutedDimsArray,
    LinearAlgebra.Transpose,
  },
) = memory_reference(CPUPointer(), parent(A))
@inline function memory_reference(::CPUPointer, A::Base.ReinterpretArray{T}) where {T}
  p, m = memory_reference(CPUPointer(), parent(A))
  reinterpret(Ptr{T}, p), m
end
@inline ind_diff(::Base.Slice, ::Any) = Zero()
@inline ind_diff(x::AbstractRange, o) = static_first(x) - o
@inline ind_diff(x::IntegerTypes, o) = x - o
@inline memory_reference(::CPUPointer, A::SubArray) =
  memory_reference_subarray(CPUPointer(), A)
@inline memory_reference(::CPUTuple, A::SubArray) = memory_reference_subarray(CPUTuple(), A)
@inline function memory_reference_subarray(::PT, A::SubArray) where {PT}
  p, m = memory_reference(PT(), parent(A))
  pA = parent(A)
  offset = StaticArrayInterface.reduce_tup(
    +,
    _map(*, _map(ind_diff, A.indices, offsets(pA)), static_strides(pA)),
  )
  if offset isa Zero
    p, m
  else
    p + sizeof(eltype(A)) * offset, m
  end
end
@inline function memory_reference(::CPUTuple, A)
  r = Ref(A)
  Base.unsafe_convert(Ptr{eltype(A)}, r), r
end
@inline function memory_reference(::StaticArrayInterface.CheckParent, A)
  P = parent(A)
  if P === A
    memory_reference(StaticArrayInterface.CPUIndex(), A)
  else
    memory_reference(StaticArrayInterface.device(P), P)
  end
end
@inline memory_reference(::StaticArrayInterface.CPUIndex, A) =
  throw("Memory access for $(typeof(A)) not implemented yet.")

@inline StaticArrayInterface.contiguous_axis(
  ::Type{A},
) where {T,N,C,A<:AbstractStridedPointer{T,N,C}} = StaticInt{C}()
@inline StaticArrayInterface.contiguous_batch_size(
  ::Type{A},
) where {T,N,C,B,A<:AbstractStridedPointer{T,N,C,B}} = StaticInt{B}()
@inline StaticArrayInterface.stride_rank(
  ::Type{A},
) where {T,N,C,B,R,A<:AbstractStridedPointer{T,N,C,B,R}} = _map(StaticInt, R)
@inline memory_reference(A::AbstractStridedPointer) = pointer(A), nothing

@inline Base.eltype(::AbstractStridedPointer{T}) where {T} = T

struct StridedPointer{T,N,C,B,R,X,O} <: AbstractStridedPointer{T,N,C,B,R,X,O}
  p::Ptr{T}
  si::StrideIndex{N,R,C,X,O}
end
@inline stridedpointer(
  p::Ptr{T},
  si::StrideIndex{N,R,C,X,O},
  ::StaticInt{B},
) where {T,N,R,C,B,X,O} = StridedPointer{T,N,C,B,R,X,O}(p, si)
@inline stridedpointer(p::Ptr{T}, si::StrideIndex{N,R,C,X,O}) where {T,N,R,C,X,O} =
  StridedPointer{T,N,C,0,R,X,O}(p, si)

@inline bytestrideindex(A::AbstractArray{T}) where {T} = mulstrides(T, StrideIndex(A))

@inline function mulstrides(::Type{T}, si::StrideIndex{N,R,C,X,O}) where {N,R,C,X,O,T}
  StrideIndex{N,R,C}(mulsizeof(T, si.strides), si.offsets)
end
@inline function stridedpointer(A::AbstractArray)
  p, r = memory_reference(A)
  stridedpointer(p, bytestrideindex(A), StaticArrayInterface.contiguous_batch_size(A))
end
@inline function stridedpointer_preserve(A::AbstractArray)
  p, r = memory_reference(A)
  stridedpointer(p, bytestrideindex(A), StaticArrayInterface.contiguous_batch_size(A)), r
end
@inline function stridedpointer_preserve(t::NTuple)
  p, r = memory_reference(t)
  stridedpointer(
    p,
    StaticArrayInterface.StrideIndex{1,(1,),1}((static(sizeof(eltype(t))),), (static(1),)),
    static(0),
  ),
  r
end
@inline val_stride_rank(::AbstractStridedPointer{T,N,C,B,R}) where {T,N,C,B,R} = Val{R}()
@generated val_dense_dims(::AbstractStridedPointer{T,N}) where {T,N} =
  Val{ntuple(==(0), Val(N))}()
@inline val_stride_rank(A) = Val(known(stride_rank(A)))
@inline val_dense_dims(A) = Val(known(StaticArrayInterface.dense_dims(A)))

function zerotupleexpr(N::Int)
  t = Expr(:tuple)
  for n = 1:N
    push!(t.args, :(Zero()))
  end
  Expr(:block, Expr(:meta, :inline), t)
end
@generated zerotuple(::Val{N}) where {N} = zerotupleexpr(N)
# @inline center(p::AbstractStridedPointer{T,N}) where {T,N} = gesp(p, zerotuple(Val(N)))
@inline zero_offsets(si::StrideIndex{N,R,C}) where {N,R,C} =
  StrideIndex{N,R,C}(getfield(si, :strides), zerotuple(Val(N)))
@inline zero_offsets(sptr::AbstractStridedPointer) = stridedpointer(
  pointer(sptr),
  zero_offsets(StrideIndex(sptr)),
  contiguous_batch_size(sptr),
)
@inline zstridedpointer(A) = zero_offsets(stridedpointer(A))
@inline function zstridedpointer_preserve(A::AbstractArray{T,N}) where {T,N}
  strd = mulsizeof(T, StaticArrayInterface.static_strides(A))
  si = StrideIndex{N,known(stride_rank(A)),Int(contiguous_axis(A))}(
    strd,
    zerotuple(Val(N)),
    Zero(),
  )
  p, r = memory_reference(A)
  stridedpointer(p, si, contiguous_batch_size(A)), r
end

@inline Base.similar(sptr::StridedPointer, ptr::Ptr) =
  stridedpointer(ptr, StrideIndex(sptr), contiguous_batch_size(sptr))

Base.unsafe_convert(::Type{Ptr{T}}, ptr::AbstractStridedPointer{T}) where {T} = pointer(ptr)
# Shouldn't need to special case Array
@generated function nopromote_axis_indicator(::AbstractStridedPointer{<:Any,N}) where {N}
  t = Expr(:tuple)
  foreach(n -> push!(t.args, True()), 1:N)
  Expr(:block, Expr(:meta, :inline), t)
end

@inline dynamic_offsets(si::StrideIndex{N,R,C}) where {N,R,C} =
  StrideIndex{N,R,C}(static_strides(si), _map(Int, offsets(si)))
struct StridedBitPointer{N,C,B,R,X,O} <: AbstractStridedPointer{Bit,N,C,B,R,X,O}
  p::Ptr{Bit}
  si::StrideIndex{N,R,C,X,O}
end
@inline stridedpointer(
  p::Ptr{Bit},
  si::StrideIndex{N,R,C,X,O},
  ::StaticInt{B},
) where {N,R,C,B,X,O} = StridedBitPointer{N,C,B,R,X,NTuple{N,Int}}(p, dynamic_offsets(si))
@inline stridedpointer(p::Ptr{Bit}, si::StrideIndex{N,R,C,X,O}) where {N,R,C,X,O} =
  StridedBitPointer{N,C,0,R,X,NTuple{N,Int}}(p, dynamic_offsets(si))

@inline Base.pointer(p::Union{StridedPointer,StridedBitPointer}) = getfield(p, :p)
@inline StaticArrayInterface.StrideIndex(sptr::Union{StridedPointer,StridedBitPointer}) =
  getfield(sptr, :si)
@inline bytestrideindex(sptr::AbstractStridedPointer) = StrideIndex(sptr)

@inline bytestrides(si::StrideIndex) =
  _map(Base.Fix2(*, StaticInt{8}()), getfield(si, :strides))
@inline bytestrides(ptr::AbstractStridedPointer) = getfield(StrideIndex(ptr), :strides)
@inline Base.strides(ptr::AbstractStridedPointer) = getfield(StrideIndex(ptr), :strides)
@inline StaticArrayInterface.static_strides(ptr::AbstractStridedPointer) =
  getfield(StrideIndex(ptr), :strides)
@inline StaticArrayInterface.offsets(ptr::AbstractStridedPointer) =
  offsets(StrideIndex(ptr))
@inline StaticArrayInterface.contiguous_axis_indicator(
  ptr::AbstractStridedPointer{T,N,C},
) where {T,N,C} = contiguous_axis_indicator(StaticInt{C}(), Val{N}())


@inline function similar_with_offset(
  sptr::AbstractStridedPointer{T,N,C,B,R},
  ptr::Ptr,
  offset::Tuple,
) where {T,N,C,B,R}
  si = StrideIndex{N,R,C}(static_strides(sptr), offset)
  stridedpointer(ptr, si, contiguous_batch_size(sptr))
end
@inline function similar_no_offset(
  sptr::AbstractStridedPointer{T,N,C,B,R},
  ptr::Ptr,
) where {T,N,C,B,R}
  si = StrideIndex{N,R,C}(static_strides(sptr), zerotuple(Val(N)))
  stridedpointer(ptr, si, contiguous_batch_size(sptr))
end

# There is probably a smarter way to do indexing adjustment here.
# The reasoning for the current approach of geping for Zero() on extracted inds
# and for offsets otherwise is best demonstrated witht his motivational example:
#
# A = OffsetArray(rand(10,11), 6:15, 5:15);
# for i in 6:15
# s += A[i,i]
# end
# first access is at zero-based index
# (first(6:16) - StaticArrayInterface.offsets(a)[1]) * StaticArrayInterface.static_strides(A)[1] + (first(6:16) - StaticArrayInterface.offsets(a)[2]) * StaticArrayInterface.static_strides(A)[2]
# equal to
#  (6 - 6)*1 + (6 - 5)*10 = 10
# i.e., the 1-based index 11.
# So now we want to adjust the offsets and pointer's value
# so that indexing it with a single `i` (after summing strides) is correct.
# Moving to 0-based indexing is easiest for combining strides. So we gep by 0 on these inds.
# E.g., gep(stridedpointer(A), (0,0))
# ptr += (0 - 6)*1 + (0 - 5)*10 = -56
# new stride = 1 + 10 = 11
# new_offset = 0
# now if we access the 6th element
# (6 - new_offse) * new_stride
# ptr += (6 - 0) * 11 = 66
# cumulative:
# ptr = pointer(A) + 66 - 56  = pointer(A) + 10
# so initial load is of pointer(A) + 10 -> the 11th element w/ 1-based indexing

@inline stridedpointer(ptr::AbstractStridedPointer) = ptr
@inline stridedpointer_preserve(ptr::AbstractStridedPointer) = (ptr, nothing)

struct FastRange{T,F,S,O}# <: AbstractStridedPointer{T,1,1,0,(1,),Tuple{S},Tuple{O}}# <: AbstractRange{T}
  f::F
  s::S
  o::O
end
FastRange{T}(f::F, s::S) where {T<:IntegerTypes,F,S} = FastRange{T,Zero,S,F}(Zero(), s, f)
FastRange{T}(f::F, s::S, o::O) where {T,F,S,O} = FastRange{T,F,S,O}(f, s, o)

FastRange{T}(f::F, s::S) where {T<:FloatingTypes,F,S} = FastRange{T,F,S,Int}(f, s, 0)
# FastRange{T}(f,s) where {T<:FloatingTypes} = FastRange{T}(f,s,fast_int64_to_double())
# FastRange{T}(f::F,s::S,::True) where {T<:FloatingTypes,F,S} = FastRange{T,F,S,Int}(f,s,0)
# FastRange{T}(f::F,s::S,::False) where {T<:FloatingTypes,F,S} = FastRange{T,F,S,Int32}(f,s,zero(Int32))

@inline function memory_reference(r::AbstractRange{T}) where {T}
  s = StaticArrayInterface.static_step(r)
  FastRange{T}(StaticArrayInterface.static_first(r) - s, s), nothing
end
@inline memory_reference(r::FastRange) = (r, nothing)
@inline bytestrides(::FastRange{T}) where {T} = (StaticInt(sizeof(T)),)
@inline StaticArrayInterface.offsets(::FastRange) = (One(),)
@inline val_stride_rank(::FastRange) = Val{(1,)}()
@inline val_dense_dims(::FastRange) = Val{(true,)}()
@inline StaticArrayInterface.contiguous_axis(::FastRange) = One()
@inline StaticArrayInterface.contiguous_batch_size(::FastRange) = Zero()

@inline stridedpointer(fr::FastRange, ::StrideIndex, ::StaticInt{0}) = fr
struct NoStrides end
@inline bytestrideindex(::FastRange) = NoStrides()
@inline StaticArrayInterface.offsets(::NoStrides) = NoStrides()
@inline reconstruct_ptr(r::FastRange{T}, o) where {T} =
  FastRange{T}(getfield(r, :f), getfield(r, :s), o)


@inline function zero_offsets(fr::FastRange{T,Zero}) where {T<:IntegerTypes}
  s = getfield(fr, :s)
  FastRange{T}(Zero(), s, getfield(fr, :o) + s)
end
@inline function zero_offsets(fr::FastRange{T}) where {T<:FloatingTypes}
  FastRange{T}(getfield(fr, :f), getfield(fr, :s), getfield(fr, :o) + 1)
end

# `FastRange{<:FloatingTypes}` must use an integer offset because `ptrforcomparison` needs to be exact/integral.

# discard unnueeded align/reg size info
@inline Base.eltype(::FastRange{T}) where {T} = T

# @inline Base.pointer(p::AbstractStridedPointer, i::Tuple) = pointer(p) + StrideIndex(p)[i...]
function reinterpret_strideindex_quote(
  sz_new::Int,
  sz_old::Int,
  C::Int,
  N::Int,
  curly::Expr,
)
  if sz_old == sz_new
    # size_expr = :size_A
    bx_expr = :bx
  else
    @assert 1 ≤ C ≤ N
    # size_expr = Expr(:tuple)
    bx_expr = Expr(:tuple)
    for n ∈ 1:N
      # sz_n = Expr(:call, GlobalRef(Core,:getfield), :size_A, n, false)
      bx_n = Expr(:call, GlobalRef(Core, :getfield), :bx, n, false)
      if n ≠ C
        # push!(size_expr.args, sz_n)
        push!(bx_expr.args, bx_n)
      elseif sz_old > sz_new
        si = :(StaticInt{$(sz_old ÷ sz_new)}())
        # push!(size_expr.args, Expr(:call, :*, sz_n, si))
        push!(bx_expr.args, Expr(:call, :÷, bx_n, si))
      else
        si = :(StaticInt{$(sz_new ÷ sz_old)}())
        # push!(size_expr.args, Expr(:call, :÷, sz_n, si))
        push!(bx_expr.args, Expr(:call, :*, bx_n, si))
      end
    end
  end
  quote
    $(Expr(:meta, :inline))
    bx = si.strides
    $curly($bx_expr, si.offsets)
  end
end
@generated function reinterpret_strideindex(
  ::Type{Tnew},
  ::Type{Told},
  si::StrideIndex{N,R,C},
) where {Tnew,Told,N,R,C}
  # 1+2
  reinterpret_strideindex_quote(
    offsetsize(Tnew),
    offsetsize(Told),
    C,
    N,
    :(StrideIndex{$N,$R,$C}),
  )
end
@generated function Base.reinterpret(
  ::Type{Tnew},
  sp::AbstractStridedPointer{Told,N,C,B,R},
) where {Tnew,Told,N,C,B,R}
  # 1+2
  sz_old = offsetsize(Told)
  sz_new = offsetsize(Tnew)
  Bnew = ((B > 0) & (sz_new ≠ sz_old)) ? ((B * sz_old) ÷ sz_new) : B
  quote
    $(Expr(:meta, :inline))
    si = reinterpret_strideindex(Tnew, Told, StrideIndex(sp))
    stridedpointer(reinterpret(Ptr{$Tnew}, pointer(sp)), si, StaticInt{$Bnew}())
  end
end

@inline vload(::StaticInt{N}, args...) where {N} = StaticInt{N}()
@inline stridedpointer(::StaticInt{N}) where {N} = StaticInt{N}()
@inline zero_offsets(::StaticInt{N}) where {N} = StaticInt{N}()
