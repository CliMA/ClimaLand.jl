"""
P are the pointers
I contains indexes into `strides` (dynamic) and `X` (static)
X static strides
"""
struct GroupedStridedPointers{P,C,B,R,I,X,O}
  ptrs::P
  strides::X
  offsets::O
end

@inline function GroupedStridedPointers{P,C,B,R,I}(
  ptrs::P,
  strides::X,
  offsets::O,
) where {P,C,B,R,I,X,O}
  GroupedStridedPointers{P,C,B,R,I,X,O}(ptrs, strides, offsets)
end

@inline function map_mem_ref(A::Tuple{P}) where {P}
  p, r = memory_reference(first(A))
  (p,), (r,)
end
@inline function map_mem_ref(A::Tuple{P,T,Vararg}) where {P,T}
  p, r = memory_reference(first(A))
  pt, rt = map_mem_ref(Base.tail(A))
  (p, pt...), (r, rt...)
end

struct DensePointerWrapper{D,T,N,C,B,R,X,O,P<:AbstractStridedPointer{T,N,C,B,R,X,O}} <:
       AbstractStridedPointer{T,N,C,B,R,X,O}
  p::P
end

@inline val_dense_dims(::DensePointerWrapper{D}) where {D} = Val{D}()

@inline Base.pointer(A::DensePointerWrapper) = pointer(getfield(A, :p))
@inline StaticArrayInterface.StrideIndex(sptr::DensePointerWrapper) =
  StrideIndex(getfield(sptr, :p))


@inline DensePointerWrapper{D}(
  sp::P,
) where {D,T,N,C,B,R,X,O,P<:AbstractStridedPointer{T,N,C,B,R,X,O}} =
  DensePointerWrapper{D,T,N,C,B,R,X,O,P}(sp)

@inline _gp_strides(x::StrideIndex) = getfield(x, :strides)
@inline _gp_strides(x) = static_strides(x)
@inline _gp_strides(::NoStrides) = NoStrides()
grouped_strided_pointer(::Tuple{}, ::Val{()}) = ((), ())

"""
G is a tuple(tuple((A_ind,A's dim),(A_ind,A's dim)), ())
it gives the groups.
"""
@inline function grouped_strided_pointer(
  A::Tuple{Vararg{Union{AbstractArray,AbstractStridedPointer,FastRange},N}},
  ::Val{G},
) where {N,G}
  m, r = map_mem_ref(A)
  sis = _map(bytestrideindex, A)
  grouped_strided_pointer(
    m,
    _map(contiguous_axis, A),
    _map(contiguous_batch_size, A),
    _map(val_stride_rank, A),
    _map(_gp_strides, sis),
    _map(offsets, sis),
    _map(val_dense_dims, A),
    Val{G}(),
  ),
  r
end

@generated function grouped_strided_pointer(
  ptrs::P,
  contig_axis::C,
  batch_sz::B,
  r::R,
  x::X,
  o::O,
  d::D,
  ::Val{()},
) where {P,C,B,R,X,O,D}
  # no groups
  N = length(P.parameters)
  q = Expr(:block, Expr(:meta, :inline))
  i = 0
  Ct = Expr(:tuple)
  Bt = Expr(:tuple)
  Rt = Expr(:tuple)
  Xt = Expr(:tuple)
  Ot = Expr(:tuple)
  It = Expr(:tuple)
  for n = 1:N
    push!(Ct.args, C.parameters[n].parameters[1])
    push!(Bt.args, B.parameters[n].parameters[1])
    push!(Rt.args, R.parameters[n].parameters[1])
    Xₙ = X.parameters[n]# Oₙ = O.parameters[n];
    Itt = Expr(:tuple)
    if !iszero(length(Xₙ.parameters))
      xₙ = Symbol(:x_, n)
      oₙ = Symbol(:o_, n)
      push!(q.args, Expr(:(=), xₙ, Expr(:call, Core.getfield, :x, n, false)))
      push!(q.args, Expr(:(=), oₙ, Expr(:call, Core.getfield, :o, n, false)))
      for j ∈ 1:length(Xₙ.parameters)
        push!(Xt.args, Expr(:call, Core.getfield, xₙ, j, false))
        push!(Ot.args, Expr(:call, Core.getfield, oₙ, j, false))
        push!(Itt.args, (i += 1))
      end
    end
    push!(It.args, Itt)
  end
  push!(q.args, :(GroupedStridedPointers{$P,$Ct,$Bt,$Rt,$It}(ptrs, $Xt, $Ot)))
  q
end

function ordered_rank_and_sort(R::NTuple{N,Int}) where {N}
  sp = ntuple(zero, Val{N}())
  r = ntuple(n -> sum(R[n] .≥ R), Val{N}())
  @inbounds for n = 1:N
    sp = Base.setindex(sp, n, r[n])
  end
  (r, sp)::NTuple{2,NTuple{N,Int}}
end

function matching_values(Xₙ, j, Xₚ, k)::Bool
  xₙⱼ = Xₙ[j]
  xₙⱼ === Int && return false
  xₚₖ = Xₚ[k]
  xₚₖ === Int && return false
  xₙⱼ === xₚₖ
end

function check_match(pm, Rₙ, Sₙ, Xₙ, Dₙ, j, Rₚ, Sₚ, Xₚ, Dₚ, k)::Bool
  # if strides match, we're done
  # @show Xₙ, j, Xₚ, k
  matching_values(Xₙ, j, Xₚ, k) && return true
  Rₙⱼ = Rₙ[j]
  Rₚₖ = Rₚ[k]
  if Rₚₖ > 1 && Rₙⱼ > 1
    # Otherwise, we check preceding axis
    nind = Sₙ[Rₙⱼ-1]
    pind = Sₚ[Rₚₖ-1]
    # matching_values(Xₙ, nind, Xₚ, pind) && return true
    # for them being of equal size,
    # @show pm[nind, pind], Dₙ[j], Dₚ[k], nind, pind
    if pm[nind, pind] && Dₙ[j] && Dₚ[k]
      return check_match(pm, Rₙ, Sₙ, Xₙ, Dₙ, nind, Rₚ, Sₚ, Xₚ, Dₚ, pind)
    end
  end
  false
end
# _typesize(@nospecialize(_::Type{Ptr{T}})) where {T} = (Base.allocatedinline(T) ? sizeof(T) : sizeof(Int))::Int
# _typesize(@nospecialize(_)) = -1
@generated function grouped_strided_pointer(
  ptrs::P,
  contig_axis::C,
  batch_sz::B,
  r::R,
  x::X,
  o::O,
  d::D,
  ::Val{G},
) where {P,C,B,R,X,O,D,G}
  # 1+1
  N = length(P.parameters)
  # We search for stride matches here
  # first we go over the groups, looking for matches
  m = Matrix{Matrix{Bool}}(undef, N, N)
  for g ∈ G
    # each `g` is a tuple of (array ind, paired dim)
    ng = length(g)
    for i ∈ 1:ng
      ai, di = g[i]
      # @show ai, di
      idim = length(X.parameters[ai].parameters)
      for j ∈ i+1:ng
        aj, dj = g[j]
        # @show aj, dj
        jdim = length(X.parameters[aj].parameters)
        pm = if isassigned(m, LinearIndices(m)[ai, aj])
          m[ai, aj]
        elseif ai < aj
          m[ai, aj] = m[aj, ai] = fill(false, jdim, idim)
          # m[aj,ai] = fill(false, idim, jdim)
        elseif ai > aj
          m[ai, aj] = m[aj, ai] = fill(false, idim, jdim)
          # m[ai,aj] = fill(false, jdim, idim)
        end
        if ai < aj
          if checkbounds(Bool, pm, dj, di)
            @inbounds pm[dj, di] = true
          end
        elseif ai > aj
          if checkbounds(Bool, pm, di, dj)
            @inbounds pm[di, dj] = true
          end
        end
      end
    end
  end
  # eltypesizes = Vector{Int}(undef, N)
  # for n in 1:N
  #   eltypesizes[n] = _typesize(P.parameters[n])
  # end
  # @show m G
  # Here we check for static size equivalence info
  # m contains information on matching sizes
  q = Expr(:block, Expr(:meta, :inline))
  i = 0
  staticxs = Bool[]
  staticos = Bool[]
  Ct = Expr(:tuple)
  Bt = Expr(:tuple)
  Rt = Expr(:tuple)
  Xt = Expr(:tuple)
  Ot = Expr(:tuple)
  It = Expr(:tuple)
  for n = 1:N
    push!(Ct.args, C.parameters[n].parameters[1])
    push!(Bt.args, B.parameters[n].parameters[1])
    Rₙ = R.parameters[n].parameters[1]
    push!(Rt.args, Rₙ)
    Xₙ = X.parameters[n].parameters
    Itt = Expr(:tuple)
    ndim = length(Xₙ)
    if ndim > 0
      Rₙ, Sₙ = ordered_rank_and_sort(Rₙ)
      Dₙ = D.parameters[n].parameters[1]
      xₙ = Symbol(:x_, n)
      oₙ = Symbol(:o_, n)
      xₙ_not_extracted = true
      oₙ_not_extracted = true
      for j = 1:ndim
        # Here, we now check if we actually need to add info, or if we can find it
        match = false
        Rₙⱼ = (Rₙ[j])::Int
        # nprev_max = isone(Rₙⱼ) ? 0 : n - 1
        for nprev = 1:n-1#ifelse(eltypesizes[n] == -1, 0, n - 1)
          # eltypematch = (eltypesizes[nprev] == eltypesizes[n])
          # @show n, nprev
          isassigned(m, LinearIndices(m)[n, nprev]) || continue
          pm = m[n, nprev] # just accessing lower triangle, but matrix is symmetric
          # pm is ndim_n x ndim_nprev
          # we need to search back through stride ranks...
          for k ∈ axes(pm, 2)
            Rₚ, Sₚ = ordered_rank_and_sort(R.parameters[nprev].parameters[1])
            Dₚ = D.parameters[nprev].parameters[1]
            Xₚ = X.parameters[nprev].parameters
            match |= check_match(pm, Rₙ, Sₙ, Xₙ, Dₙ, j, Rₚ, Sₚ, Xₚ, Dₚ, k)
            if match
              # must push nprev's ind
              storedindex = It.args[nprev].args[k]
              push!(Itt.args, storedindex)
              if (!staticxs[storedindex]) && (Xₙ[j] <: StaticInt)
                Xt.args[storedindex] = Expr(:call, Xₙ[j])
                staticxs[storedindex] = true
              end
              if (!staticos[storedindex]) && (O.parameters[n].parameters[j] <: StaticInt)
                Ot.args[storedindex] = Expr(:call, O.parameters[n].parameters[j])
                staticos[storedindex] = true
              end
              break
            end
          end
          match && break
        end
        if !match
          # @show Xₙ[j]
          xparam = Xₙ[j]
          xstatic = xparam <: StaticInt
          push!(staticxs, xstatic)
          if xstatic
            push!(Xt.args, Expr(:call, xparam))
          else
            if xₙ_not_extracted
              push!(q.args, Expr(:(=), xₙ, Expr(:call, Core.getfield, :x, n, false)))
              xₙ_not_extracted = false
            end
            push!(Xt.args, Expr(:call, Core.getfield, xₙ, j, false))
          end
          oparam = O.parameters[n].parameters[j]
          ostatic = oparam <: StaticInt
          push!(staticos, ostatic)
          if ostatic
            push!(Ot.args, Expr(:call, oparam))
          else
            if oₙ_not_extracted
              push!(q.args, Expr(:(=), oₙ, Expr(:call, Core.getfield, :o, n, false)))
              oₙ_not_extracted = false
            end
            push!(Ot.args, Expr(:call, Core.getfield, oₙ, j, false))
          end
          push!(Itt.args, (i += 1))
        end
      end
    end
    push!(It.args, Itt)
  end
  push!(q.args, :(GroupedStridedPointers{$P,$Ct,$Bt,$Rt,$It}(ptrs, $Xt, $Ot)))
  q
end

@generated function stridedpointers(
  gsp::GroupedStridedPointers{P,C,B,R,I,X,O},
) where {P,C,B,R,N,X<:Tuple{Vararg{Any,N}},O<:Tuple{Vararg{Any,N}},I}
  t = Expr(:tuple)
  gf = Core.getfield
  q = Expr(
    :block,
    Expr(:meta, :inline),
    :(ptrs = $gf(gsp, :ptrs)),
    :(strds = $gf(gsp, :strides)),
    :(offs = $gf(gsp, :offsets)),
  )
  for n ∈ 1:N
    push!(
      q.args,
      Expr(:(=), Symbol(:x_, n), Expr(:call, gf, :strds, n, false)),
      Expr(:(=), Symbol(:o_, n), Expr(:call, gf, :offs, n, false)),
    )
  end
  for i ∈ eachindex(I)
    p = Expr(:call, gf, :ptrs, i, false)
    if P.parameters[i] <: FastRange
      push!(t.args, p)
      continue
    end
    Iᵢ = I[i]
    Nᵢ = length(Iᵢ)
    x = Expr(:tuple)
    o = Expr(:tuple)
    for j ∈ Iᵢ
      push!(x.args, Symbol(:x_, j))
      push!(o.args, Symbol(:o_, j))
    end
    si = Expr(:call, Expr(:curly, :StrideIndex, Nᵢ, R[i], C[i]), x, o)
    push!(t.args, Expr(:call, :stridedpointer, p, si, :(StaticInt{$(B[i])}())))
  end
  push!(q.args, t)
  return q
end
