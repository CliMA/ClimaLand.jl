@inline undef_memory_buffer(::Type{T}, ::StaticInt{L}) where {T,L} =
  MemoryBuffer{L,T}(undef)
@inline undef_memory_buffer(::Type{T}, L) where {T} = Vector{T}(undef, L)
@inline undef_memory_buffer(::Type{Bit}, ::StaticInt{L}) where {L} =
  MemoryBuffer{(L + 7) >>> 3,UInt8}(undef)
@inline undef_memory_buffer(::Type{Bit}, L) =
  Vector{UInt8}(undef, (L + 7) >>> 3)
# @inline undef_memory_buffer(::Type{Bit}, L) = BitVector(undef, L)

struct AbstractStrideArrayImpl{T,N,R,S,X,O,P,A} <:
       AbstractStrideArray{T,N,R,S,X,O}
  ptr::AbstractPtrArray{T,N,R,S,X,O,P}
  data::A
end
const StrideArray{T,N,R,S,X,O,A} = AbstractStrideArrayImpl{T,N,R,S,X,O,T,A}
const StrideArray0{T,N,R,S,X,A} =
  AbstractStrideArrayImpl{T,N,R,S,X,NTuple{N,Zero},T,A}
const StrideArray1{T,N,R,S,X,A} =
  AbstractStrideArrayImpl{T,N,R,S,X,NTuple{N,One},T,A}
const BitStrideArray{N,R,S,X,O,A} =
  AbstractStrideArrayImpl{Bool,N,R,S,X,O,Bit,A}
const BitStrideArray0{N,R,S,X,A} =
  AbstractStrideArrayImpl{Bool,N,R,S,X,NTuple{N,Zero},Bit,A}
const BitStrideArray1{N,R,S,X,A} =
  AbstractStrideArrayImpl{Bool,N,R,S,X,NTuple{N,One},Bit,A}

const StrideVector{T,R,S,X,O,A} = StrideArray{T,1,R,S,X,O,A}
const StrideMatrix{T,R,S,X,O,A} = StrideArray{T,2,R,S,X,O,A}

@inline function StrideArray(
  B::AbstractPtrArray{T,N,R,S,X,O,P},
  data::A
) where {T,N,R,S,X,O,A,P}
  AbstractStrideArrayImpl{T,N,R,S,X,O,P,A}(B, data)
end
function StrideArray(::AbstractPtrArray, ::Integer)
  throw("StrideArray(::AbstractPtrArray, ::Integer) doesn't make any sense.")
end
function StrideArray(::AbstractPtrArray, ::Type{T}) where {T}
  throw("StrideArray(::AbstractPtrArray, ::Type{$T}) doesn't make any sense.")
end
function StrideArray(::AbstractPtrArray, ::Tuple{Vararg{Integer}})
  throw(
    "StrideArray(::AbstractPtrArray, ::Tuple{Vararg{Integer}}) doesn't make any sense."
  )
end

@inline Base.pointer(A::AbstractStrideArrayImpl) = pointer(getfield(A, :ptr))
@inline PtrArray(A::AbstractStrideArrayImpl) = getfield(A, :ptr)

@inline StrideArray(A::AbstractArray) = StrideArray(PtrArray(A), A)

@inline function StrideArray{T}(
  ::UndefInitializer,
  s::Tuple{Vararg{Union{Integer,StaticInt},N}}
) where {N,T}
  L = Static.reduce_tup(*, s)
  b = undef_memory_buffer(T, L)
  # For now, just trust Julia's alignment heuristics are doing the right thing
  # to save us from having to over-allocate
  # ptr = LayoutPointers.align(pointer(b))
  A = PtrArray(
    pointer(b),
    s,
    ntuple(Returns(nothing), Val(N)),
    ntuple(Returns(One()), Val(N)),
    Val(ntuple(identity, Val(N)))
  )
  StrideArray(A, b)
end
@inline function StrideArray(
  ptr::Ptr{T},
  s::S,
  x::X,
  b,
  ::Val{D}
) where {S,X,T,D}
  StrideArray(PtrArray(ptr, s, x, Val{D}()), b)
end
@inline StrideArray(f, s::Vararg{Union{Integer,StaticInt},N}) where {N} =
  StrideArray{Float64}(f, s)
@inline StrideArray(f, s::Tuple{Vararg{Union{Integer,StaticInt},N}}) where {N} =
  StrideArray{Float64}(f, s)
@inline StrideArray(
  f,
  ::Type{T},
  s::Vararg{Union{Integer,StaticInt},N}
) where {T,N} = StrideArray{T}(f, s)
@inline StrideArray{T}(f, s::Vararg{Union{Integer,StaticInt},N}) where {T,N} =
  StrideArray{T}(f, s)
@inline function StrideArray(
  A::PtrArray{T,N},
  s::Tuple{Vararg{Union{Integer,StaticInt},N}}
) where {T,N}
  PtrArray(stridedpointer(A), s, val_dense_dims(A))
end
@inline function StrideArray(
  A::AbstractArray{T,N},
  sz::Tuple{Vararg{Union{Integer,StaticInt},N}}
) where {T,N}
  # what is the point of this method? Why not `view(PtrArray(A), ...)`?
  p = pointer(A)
  sx = _sparse_strides(dense_dims(A), static_strides(A))
  R = map(Int, stride_rank(A))
  B = PtrArray(p, sz, sx, offsets(A), Val(R))
  StrideArray(B, preserve_buffer(A))
end
@inline function StrideArray{T}(
  f,
  s::Tuple{Vararg{Union{Integer,StaticInt},N}}
) where {T,N}
  A = StrideArray{T}(undef, s)
  @inbounds for i ∈ eachindex(A)
    A[i] = f(T)
  end
  A
end
@inline function StrideArray{T}(
  ::typeof(zero),
  s::Tuple{Vararg{Union{Integer,StaticInt},N}}
) where {T,N}
  ptr = Ptr{T}(Libc.calloc(prod(map(Int, s)), sizeof(T)))
  A = unsafe_wrap(Array{T}, ptr, s; own = true)
  StrideArray(A)
end
mutable struct StaticStrideArray{T,N,R,S,X,O,L} <:
               AbstractStrideArray{T,N,R,S,X,O}
  data::NTuple{L,T}
  @inline StaticStrideArray{T,N,R,S,X,O}(
    ::UndefInitializer
  ) where {T,N,R,S,X,O} =
    new{T,N,R,S,X,O,Int(ArrayInterface.reduce_tup(*, to_static_tuple(Val(S))))}()
  @inline StaticStrideArray{T,N,R,S,X,O,L}(
    ::UndefInitializer
  ) where {T,N,R,S,X,O,L} = new{T,N,R,S,X,O,L}()
end

@generated function to_static_tuple(::Val{S}) where {S}
  t = Expr(:tuple)
  for s ∈ S.parameters
    push!(t.args, Expr(:new, s))
  end
  t
end
@inline StrideArray(A::StaticStrideArray) = A

@inline ArrayInterface.static_size(
  ::StaticStrideArray{<:Any,<:Any,<:Any,S}
) where {S} = to_static_tuple(Val(S))
@inline function ArrayInterface.static_strides(
  A::StaticStrideArray{<:Any,N,R}
) where {N,R}
  _strides_entry(
    static_size(A),
    ntuple(Returns(nothing), Val{N}()),
    Val{R}(),
    Val(false)
  )
end
@inline ArrayInterface.offsets(
  ::StaticStrideArray{T,N,R,S,X,O}
) where {T,N,R,S,X,O} = to_static_tuple(Val(O))
@inline Base.unsafe_convert(::Type{Ptr{T}}, A::StaticStrideArray{T}) where {T} =
  Base.unsafe_convert(Ptr{T}, pointer_from_objref(A))
@inline Base.pointer(A::StaticStrideArray{T}) where {T} =
  Base.unsafe_convert(Ptr{T}, pointer_from_objref(A))

@inline StrideArray{T}(
  ::UndefInitializer,
  s::Tuple{Vararg{StaticInt,N}}
) where {N,T} = StaticStrideArray{T}(undef, s)
@inline StrideArray{T}(f, s::Tuple{Vararg{StaticInt,N}}) where {N,T} =
  StaticStrideArray{T}(f, s)
@inline StrideArray{T}(
  ::typeof(zero),
  s::Tuple{Vararg{StaticInt,N}}
) where {N,T} = StaticStrideArray{T}(zero, s) # Eager when static; assumed small

@inline function StaticStrideArray{T}(
  ::UndefInitializer,
  s::S
) where {N,T,S<:Tuple{Vararg{StaticInt,N}}}
  StaticStrideArray{
    T,
    N,
    ntuple(identity, Val(N)),
    S,
    NTuple{N,Nothing},
    NTuple{N,One}
  }(
    undef
  )
end

@inline function StaticStrideArray{T}(
  f,
  s::Tuple{Vararg{StaticInt,N}}
) where {T,N}
  A = StaticStrideArray{T}(undef, s)
  @inbounds for i ∈ eachindex(A)
    A[i] = f(T)
  end
  A
end

@inline LayoutPointers.preserve_buffer(A::MemoryBuffer) = A
@inline LayoutPointers.preserve_buffer(A::StrideArray) =
  preserve_buffer(getfield(A, :data))

@inline PtrArray(A::StrideArray) = getfield(A, :ptr)

@inline maybe_ptr_array(A) = A
@inline maybe_ptr_array(A::AbstractArray) =
  maybe_ptr_array(ArrayInterface.device(A), A)
@inline maybe_ptr_array(::ArrayInterface.CPUPointer, A::AbstractArray) =
  PtrArray(A)
@inline maybe_ptr_array(_, A::AbstractArray) = A

@inline ArrayInterface.static_size(A::AbstractStrideArrayImpl) =
  getfield(getfield(A, :ptr), :sizes)

@inline ArrayInterface.static_strides(A::AbstractStrideArrayImpl) =
  static_strides(getfield(A, :ptr))
@inline ArrayInterface.offsets(A::AbstractStrideArrayImpl) =
  offsets(getfield(A, :ptr))

@inline zeroindex(r::ArrayInterface.OptionallyStaticUnitRange{One}) =
  CloseOpen(Zero(), last(r))
@inline zeroindex(r::Base.OneTo) = CloseOpen(Zero(), last(r))
@inline zeroindex(r::AbstractUnitRange) = Zero():(last(r)-first(r))

@inline zeroindex(r::CloseOpen{Zero}) = r
@inline zeroindex(r::ArrayInterface.OptionallyStaticUnitRange{Zero}) = r

@inline LayoutPointers.zero_offsets(A::StrideArray) =
  StrideArray(LayoutPointers.zero_offsets(PtrArray(A)), preserve_buffer(A))
@inline LayoutPointers.zero_offsets(A::StaticStrideArray) =
  StrideArray(LayoutPointers.zero_offsets(PtrArray(A)), A)

@generated rank_to_sortperm_val(::Val{R}) where {R} =
  :(Val{$(rank_to_sortperm(R))}())
@inline function similar_layout(A::AbstractStrideArray{T,N,R}) where {T,N,R}
  permutedims(similar(permutedims(A, rank_to_sortperm_val(Val{R}()))), Val{R}())
end
@inline function similar_layout(A::AbstractArray)
  b = preserve_buffer(A)
  GC.@preserve b begin
    similar_layout(PtrArray(A))
  end
end
@inline function Base.similar(A::AbstractStrideArray{T}) where {T}
  StrideArray{T}(undef, static_size(A))
end
@inline Base.similar(A::BitPtrArray) = StrideArray{Bit}(undef, static_size(A))
@inline function Base.similar(A::AbstractStrideArray, ::Type{T}) where {T}
  StrideArray{T}(undef, static_size(A))
end

@inline function Base.view(
  A::AbstractStrideArray,
  i::Vararg{Union{Integer,AbstractRange,Colon},K}
) where {K}
  StrideArray(view(PtrArray(A), i...), preserve_buffer(A))
end
@inline function Base.view(
  A::AbstractStrideArray{<:Any,K},
  ::Vararg{Colon,K}
) where {K}
  A
end

@inline function zview(
  A::AbstractStrideArray,
  i::Vararg{Union{Integer,AbstractRange,Colon},K}
) where {K}
  StrideArray(zview(PtrArray(A), i...), preserve_buffer(A))
end
@inline function Base.permutedims(A::AbstractStrideArray, ::Val{P}) where {P}
  StrideArray(permutedims(PtrArray(A), Val{P}()), preserve_buffer(A))
end
@inline Base.adjoint(a::AbstractStrideVector{<:Real}) =
  StrideArray(transpose(PtrArray(a)), preserve_buffer(a))
@inline Base.transpose(a::AbstractStrideVector) =
  StrideArray(transpose(PtrArray(a)), preserve_buffer(a))

function gc_preserve_call(ex::Expr, skip::Int = 0, escape::Bool = true)
  q = Expr(:block)
  call = Expr(:call, escape ? esc(ex.args[1]) : ex.args[1])
  gcp = Expr(:gc_preserve, call)
  if length(ex.args) >= 2 && Meta.isexpr(ex.args[2], :parameters)
    params = ex.args[2]
    deleteat!(ex.args, 2)
    for i in eachindex(params.args)
      param = params.args[i]
      if Meta.isexpr(param, :kw, 2)
        push!(ex.args, param)
      else
        push!(ex.args, Expr(:kw, param, param))
      end
    end
  end
  for i ∈ 2:length(ex.args)
    arg = ex.args[i]
    if i + 1 ≤ skip
      push!(call.args, arg)
      continue
    end
    A = gensym(:A)
    buffer = gensym(:buffer)
    if arg isa Expr && arg.head === :kw
      push!(call.args, Expr(:kw, arg.args[1], Expr(:call, maybe_ptr_array, A)))
      arg = arg.args[2]
    else
      push!(call.args, Expr(:call, maybe_ptr_array, A))
    end
    if escape
      arg = esc(arg)
    end
    push!(q.args, :($A = $(arg)))
    push!(q.args, Expr(:(=), buffer, Expr(:call, preserve_buffer, A)))
    push!(gcp.args, buffer)
  end
  push!(q.args, gcp)
  q
end
"""
@gc_preserve foo(A, B, C)

Apply to a single, non-nested, function call. It will `GC.@preserve` all the arguments, and substitute suitable arrays with `PtrArray`s.
This has the benefit of potentially allowing statically sized mutable arrays to be both stack allocated, and passed through a non-inlined function boundary.
"""
macro gc_preserve(ex)
  @assert ex.head === :call
  gc_preserve_call(ex)
end

@inline LayoutPointers.zero_offsets(A::AbstractArray) =
  StrideArray(LayoutPointers.zero_offsets(PtrArray(A)), A)

@generated function Base.IndexStyle(
  ::Type{
    <:Union{
      <:BitPtrArray{N,R,S,NTuple{N,Nothing}},
      <:StrideArray{Bit,N,R,S,NTuple{N,Nothing}}
    }
  }
) where {N,R,S}
  # if is column major || is a transposed contiguous vector
  if (
    (R === ntuple(identity, Val(N))) ||
    (R === (2, 1) && S <: Tuple{One,Integer})
  )
    if N > 1
      ks1 = known(S)[1]
      if ks1 === nothing || ((ks1 & 7) != 0)
        return :(IndexCartesian())
      end
    end
    :(IndexLinear())
  else
    :(IndexCartesian())
  end
end

@inline Base.reinterpret(::Type{T}, A::AbstractStrideArrayImpl) where {T} =
  StrideArray(reinterpret(T, PtrArray(A)), preserve_buffer(A))
@inline Base.reinterpret(
  ::typeof(reshape),
  ::Type{T},
  A::AbstractStrideArrayImpl
) where {T} =
  StrideArray(reinterpret(reshape, T, PtrArray(A)), preserve_buffer(A))

@inline function Base.reinterpret(
  ::typeof(reshape),
  ::Type{U},
  B::AbstractStrideArrayImpl{T,N,R,S,X,O,P,A}
) where {U,T,N,R,S,X,O,P,A<:StaticStrideArray{U}}
  C = B.data
  if typeof(reinterpret(reshape, T, C)) === typeof(B)
    C
  else
    StrideArray(reinterpret(reshape, T, PtrArray(B)), preserve_buffer(B))
  end
end
@inline function Base.reinterpret(
  ::Type{U},
  B::AbstractStrideArrayImpl{T,N,R,S,X,O,P,A}
) where {U,T,N,R,S,X,O,P,A<:StaticStrideArray{U}}
  C = B.data
  if typeof(reinterpret(T, C)) === typeof(B)
    C
  else
    StrideArray(reinterpret(T, PtrArray(B)), preserve_buffer(B))
  end
end

@inline Base.reinterpret(::Type{T}, A::StaticStrideArray) where {T} =
  StrideArray(reinterpret(T, PtrArray(A)), A)
@inline Base.reinterpret(
  ::typeof(reshape),
  ::Type{T},
  A::StaticStrideArray
) where {T} = StrideArray(reinterpret(reshape, T, PtrArray(A)), A)
