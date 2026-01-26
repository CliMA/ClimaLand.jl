module ManualMemory

mutable struct MemoryBuffer{N,T}
  data::NTuple{N,T}
  @inline function MemoryBuffer{N,T}(::UndefInitializer) where {N,T}
    @assert Base.allocatedinline(T)
    new{N,T}()
  end
  @inline function MemoryBuffer(data::NTuple{N,T}) where {N,T}
    @assert Base.allocatedinline(T)
    new{N,T}(data)
  end
end
@inline Base.unsafe_convert(::Type{Ptr{T}}, m::MemoryBuffer) where {T} = Ptr{T}(pointer_from_objref(m))
@inline Base.pointer(m::MemoryBuffer{N,T}) where {N,T} = Ptr{T}(pointer_from_objref(m))
@inline Base.:(==)(::MemoryBuffer, ::MemoryBuffer) = false
@inline function Base.:(==)(a::MemoryBuffer{N,A}, b::MemoryBuffer{N,B}) where {N,A,B}
    GC.@preserve a b begin
        pa = pointer(a)
        pb = pointer(b)
        for n in 0:N-1
            load(pa + n*offsetsize(A)) == load(pb + n*offsetsize(B)) || return false
        end
        return true
    end
end

"""
    PseudoPtr(data, position=firstindex(data))

Provides a convenient wrapper that functions like `pointer(data)` for instances where `data`
cannot produce a viable pointer.
"""
struct PseudoPtr{T,D} <: Ref{T}
  data::D
  position::Int

  PseudoPtr(data::D, position) where {D} = new{eltype(D),D}(data, position)
  PseudoPtr(data) = PseudoPtr(data, firstindex(data))
end
Base.:(+)(x::PseudoPtr, y::Int) = PseudoPtr(getfield(x, :data), getfield(x, :position) + y)
Base.:(+)(x::Int, y::PseudoPtr) = y + x

@inline load(x::PseudoPtr) = @inbounds(getindex(getfield(x, :data), getfield(x, :position)))
@generated function load(p::Ptr{T}) where {T}
  if Base.allocatedinline(T)
    Expr(:block, Expr(:meta,:inline), :(unsafe_load(p)))
  else
    Expr(:block, Expr(:meta,:inline), :(ccall(:jl_value_ptr, Ref{$T}, (Ptr{Cvoid},), unsafe_load(Base.unsafe_convert(Ptr{Ptr{Cvoid}}, p)))))
  end
end
@inline load(p::Ptr{UInt}, ::Type{T}) where {T} = load(p, T, 0)[2]

@inline function store!(x::PseudoPtr, val)
    @inbounds(setindex!(getfield(x, :data), val, getfield(x, :position)))
end
@generated function store!(p::Ptr{T}, v::T) where {T}
  if Base.allocatedinline(T)
    Expr(:block, Expr(:meta,:inline), :(unsafe_store!(p, v); return nothing))
  else
    Expr(:block, Expr(:meta,:inline), :(unsafe_store!(Base.unsafe_convert(Ptr{Ptr{Cvoid}}, p), Base.pointer_from_objref(v)); return nothing))
  end
end
@generated offsetsize(::Type{T}) where {T} = Base.allocatedinline(T) ? sizeof(T) : sizeof(Int)

@inline store!(p::Ptr{T}, v) where {T} = store!(p, convert(T, v))

mutable struct Reference{T}
  data::T
  Reference{T}() where {T} = new{T}()
  Reference(x) = new{typeof(x)}(x)
end
@inline Base.pointer(r::Reference{T}) where {T} = Ptr{T}(pointer_from_objref(r))
@inline load(p::Ptr{Reference{T}}) where {T} = getfield(ccall(:jl_value_ptr, Ref{Reference{T}}, (Ptr{Cvoid},), unsafe_load(Base.unsafe_convert(Ptr{Ptr{Cvoid}}, p))), :data)
@inline dereference(r::Reference) = getfield(r, :data)
@inline dereference(x) = x

"""
    LazyPreserve(x, ptrcall=nothing)

Used by [`preserve`](@ref) to identify arguments that will be unwrapped with
[`preserve_buffer`](@ref), which is in turn converted to a pointer. If `ptrcall` is
specified conversion to a pointer occurs through a call equivalent to
`ptrcall(preserve_buffer(x), x)`. Otherwise, a call equivalent to
`pointer(preserve_buffer(x))` occurs.
"""
struct LazyPreserve{A,P}
  arg::A
  ptrcall::P
end
LazyPreserve(arg) = LazyPreserve(arg, nothing)
(p::LazyPreserve)(x) = p.ptrcall(x, p.arg)
(p::LazyPreserve{A,Nothing})(x) where {A} = pointer(x)


"""
    preserve_buffer(x)

For structs wrapping arrays, using `GC.@preserve` can trigger heap allocations.
`preserve_buffer` attempts to extract the heap-allocated part. Isolating it by itself
will often allow the heap allocations to be elided. For example:

```julia
julia> using StaticArrays, BenchmarkTools
julia> # Needed until a release is made featuring https://github.com/JuliaArrays/StaticArrays.jl/commit/a0179213b741c0feebd2fc6a1101a7358a90caed
       Base.elsize(::Type{<:MArray{S,T}}) where {S,T} = sizeof(T)
julia> @noinline foo(A) = unsafe_load(A,1)
foo (generic function with 1 method)
julia> function alloc_test_1()
           A = view(MMatrix{8,8,Float64}(undef), 2:5, 3:7)
           A[begin] = 4
           GC.@preserve A foo(pointer(A))
       end
alloc_test_1 (generic function with 1 method)
julia> function alloc_test_2()
           A = view(MMatrix{8,8,Float64}(undef), 2:5, 3:7)
           A[begin] = 4
           pb = parent(A) # or `LoopVectorization.preserve_buffer(A)`; `perserve_buffer(::SubArray)` calls `parent`
           GC.@preserve pb foo(pointer(A))
       end
alloc_test_2 (generic function with 1 method)
julia> @benchmark alloc_test_1()
BenchmarkTools.Trial:
  memory estimate:  544 bytes
  allocs estimate:  1
  --------------
  minimum time:     17.227 ns (0.00% GC)
  median time:      21.352 ns (0.00% GC)
  mean time:        26.151 ns (13.33% GC)
  maximum time:     571.130 ns (78.53% GC)
  --------------
  samples:          10000
  evals/sample:     998
julia> @benchmark alloc_test_2()
BenchmarkTools.Trial:
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     3.275 ns (0.00% GC)
  median time:      3.493 ns (0.00% GC)
  mean time:        3.491 ns (0.00% GC)
  maximum time:     4.998 ns (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1000
```
"""
@inline preserve_buffer(x::LazyPreserve) = preserve_buffer(x.arg)
@inline preserve_buffer(x) = x
@inline preserve_buffer(A::AbstractArray) = _preserve_buffer(A, parent(A))
@inline _preserve_buffer(a::A, p::P) where {A,P<:AbstractArray} = _preserve_buffer(p, parent(p))
@inline _preserve_buffer(a::A, p::A) where {A<:AbstractArray} = a
@inline _preserve_buffer(a::A, p::P) where {A,P} = p

function load_aggregate(::Type{T}, offset::Int) where {T}
  numfields = fieldcount(T)
  call = (T <: Tuple) ? Expr(:tuple) : Expr(:new, T)
  for f ∈ 1:numfields
    TF = fieldtype(T, f)
    if Base.issingletontype(TF)
      push!(call.args, TF.instance)
    elseif (fieldcount(TF) ≡ 0) || (TF <: Reference)
      ptr = :(p + (offset + $offset))
      ptr = TF === UInt ? ptr : :(reinterpret(Ptr{$TF}, $ptr))
      push!(call.args, :(load($ptr)))
      offset += offsetsize(TF)
    else
      arg, offset = load_aggregate(TF, offset)
      push!(call.args, arg)
    end
  end
  return call, offset
end
@generated function load(p::Ptr{UInt}, ::Type{T}, offset::Int) where {T}
  if Base.issingletontype(T)
    call = Expr(:tuple, :offset, T.instance)
  elseif (fieldcount(T) ≡ 0) || (T <: Reference)
    ptr = :(p + offset)
    ptr = T === UInt ? ptr : :(reinterpret(Ptr{$T}, $ptr))
    call = :(((offset + $(offsetsize(T)), load($ptr))))
  else
    call, off = load_aggregate(T, 0)
    call = Expr(:tuple, :(offset + $off), call)
  end
  Expr(:block, Expr(:meta,:inline), call)
end

function store_aggregate!(q::Expr, sym, ::Type{T}, offset::Int) where {T}
  gf = GlobalRef(Core,:getfield)
  for f ∈ 1:fieldcount(T)
    TF = fieldtype(T, f)
    Base.issingletontype(TF) && continue
    gfcall = Expr(:call, gf, sym, f)
    if (fieldcount(TF) ≡ 0) || (TF <: Reference)
      ptr = :(p + (offset + $offset))
      ptr = TF === UInt ? ptr : :(reinterpret(Ptr{$TF}, $ptr))
      push!(q.args, :(store!($ptr, $gfcall)))
      offset += offsetsize(TF)
    else
      newsym = gensym(sym)
      push!(q.args, Expr(:(=), newsym, gfcall))
      offset = store_aggregate!(q, newsym, TF, offset)
    end
  end
  return offset
end
@generated function store!(p::Ptr{UInt}, x::T, offset::Int) where {T}
  Base.issingletontype(T) && return :offset
  body = Expr(:block, Expr(:meta,:inline))
  if (fieldcount(T) ≡ 0) || (T <: Reference)
    ptr = :(p + offset)
    ptr = T === UInt ? ptr : :(reinterpret(Ptr{$T}, $ptr))
    push!(body.args, :(store!($ptr, x)))
    off = offsetsize(T)
  else
    off = store_aggregate!(body, :x, T, 0)
  end
  push!(body.args, Expr(:call, +, :offset, off))
  return body
end

"""
    preserve(op, args...; kwargs...)

Searches through `args` and `kwargs` for instances of [`LazyPreserve`](@ref), which are
unwrapped using [`preserve_buffer`](@ref) and preserved from garbage collection
(`GC.@preserve`). the resulting buffers are converted to pointers and passed in order to `op`.

# Examples

```julia
julia> using ManualMemory: store!, preserve, LazyPreserve

julia> x = [0 0; 0 0];

julia> preserve(store!, LazyPreserve(x), 1)

julia> x[1]
1

```
"""
preserve(op, args...; kwargs...) = _preserve(op, args, kwargs.data)
@generated function _preserve(op, args::A, kwargs::NamedTuple{syms,K}) where {A,syms,K}
  _preserve_expr(A, syms, K)
end
function _preserve_expr(::Type{A}, syms::Tuple{Vararg{Symbol}}, ::Type{K}) where {A,K}
  body = Expr(:block, Expr(:meta,:inline))
  call = Expr(:call, :op)
  pres = :(GC.@preserve)
  @inbounds for i in 1:length(A.parameters)
    arg_i = _unwrap_preserve(body, pres, :(getfield(args, $i)), A.parameters[i])
    push!(call.args, arg_i)
  end
  if length(syms) > 0
    kwargs = Expr(:parameters)
    @inbounds for i in 1:length(syms)
      arg_i = _unwrap_preserve(body, pres, :(getfield(kwargs, $i)), K.parameters[i])
      push!(call.args, Expr(:kw, syms[i], arg_i))
    end
    push!(call.args, kwargs)
  end
  push!(pres.args, call)
  push!(body.args, pres)
  return body
end
function _unwrap_preserve(body::Expr, pres::Expr, argexpr::Expr, argtype::Type)
  if argtype <: LazyPreserve
    bufsym = gensym()
    push!(body.args, Expr(:(=), bufsym, Expr(:call, :preserve_buffer, argexpr)))
    push!(pres.args, bufsym)
    return :($argexpr($bufsym))
  else
    return argexpr
  end
end

end
