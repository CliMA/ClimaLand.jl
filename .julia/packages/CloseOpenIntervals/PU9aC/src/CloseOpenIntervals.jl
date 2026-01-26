module CloseOpenIntervals
if isdefined(Base, :Experimental) &&
   isdefined(Base.Experimental, Symbol("@max_methods"))
    @eval Base.Experimental.@max_methods 1
end

using Static: StaticInt, Zero, One
export CloseOpen, SafeCloseOpen

const IntegerType = Union{Integer,StaticInt}

abstract type AbstractCloseOpen{L <: IntegerType, U <: IntegerType} <: AbstractUnitRange{Int} end
for T ∈ (:CloseOpen,:SafeCloseOpen)
  @eval begin
    struct $T{L <: IntegerType, U <: IntegerType} <: AbstractCloseOpen{L,U}
      start::L
      upper::U
      @inline $T{L,U}(l::L,u::U) where {L <: IntegerType, U <: IntegerType} = new{L,U}(l,u)
    end
    @inline $T(s::S, u::U) where {S,U} = $T{S,U}(s, u)
    @inline $T(len::T) where {T<:IntegerType} = $T{Zero,T}(Zero(), len)
  end
end

"""
    CloseOpen([start=0], U) <: AbstractUnitRange{Int}
    SafeCloseOpen([start=0], U) <: AbstractUnitRange{Int}

Define close-open unit ranges, i.e. `CloseOpen(0,10)` iterates from from `0:9`.
Close-open ranges can be more convenient in some circumstances, e.g. when
partitioning a larger array.

```julia

function foo(x)
  nt = Threads.nthreads()
  d, r = divrem(length(x), nt)
  i = firstindex(x)
  Threads.@sync for j in 1:nt
    stop = i + d + (r >= j)
    Threads.@spawn bar!(\$(@view(x[CloseOpen(i, stop)])))
    i = stop
  end
end
```
This saves a few `-1`s on the ends of the ranges if using a normal unit range.

`SafeCloseOpen` will not iterate if `U <= start`, while `CloseOpen` doesn't check
for this, and thus the behavior is undefined in that case.
"""
AbstractCloseOpen
"""
    CloseOpen([start=0], U) <: AbstractUnitRange

Iterates over the range start:U-1. Guaranteed to iterate at least once, skipping initial empty check.
See `?AbstractCloseOpen` for more information.
"""
CloseOpen
"""
    SafeCloseOpen([start=0], U) <: AbstractUnitRange

Iterates over the range start:U-1. Will not iterate if it `isempty`, like a typical range, making it
generally the preferred choice over `CloseOpen`.
See `?AbstractCloseOpen` for more information.
"""
SafeCloseOpen

@inline Base.first(r::AbstractCloseOpen) = getfield(r,:start)
@inline Base.first(r::AbstractCloseOpen{StaticInt{F}}) where {F} = F
@inline Base.step(::AbstractCloseOpen) = 1
@inline Base.last(r::AbstractCloseOpen{<:IntegerType,I}) where {I} = getfield(r,:upper) - 1
@inline Base.last(r::AbstractCloseOpen{<:IntegerType,StaticInt{L}}) where {L} = L - 1
@inline Base.length(r::AbstractCloseOpen) = Int(getfield(r,:upper) - getfield(r,:start))
@inline Base.length(r::AbstractCloseOpen{Zero}) = Int(getfield(r,:upper))

@inline Base.iterate(r::CloseOpen) = (i = Int(first(r)); (i, i))
@inline Base.iterate(r::SafeCloseOpen) = (i = Int(first(r)); i ≥ getfield(r, :upper) ? nothing : (i, i))
@inline Base.iterate(r::AbstractCloseOpen, i::IntegerType) = (i += one(i)) ≥ getfield(r, :upper) ? nothing : (i, i)

import StaticArrayInterface
StaticArrayInterface.known_first(::Type{<:AbstractCloseOpen{StaticInt{F}}}) where {F} = F
StaticArrayInterface.known_step(::Type{<:AbstractCloseOpen}) = 1
StaticArrayInterface.known_last(::Type{<:AbstractCloseOpen{<:Any,StaticInt{L}}}) where {L} = L - 1
StaticArrayInterface.known_length(::Type{<:AbstractCloseOpen{StaticInt{F},StaticInt{L}}}) where {F,L} = L - F

Base.IteratorSize(::Type{<:AbstractCloseOpen}) = Base.HasShape{1}()
Base.IteratorEltype(::Type{<:AbstractCloseOpen}) = Base.HasEltype()
@inline Base.size(r::AbstractCloseOpen) = (length(r),)
@inline Base.eltype(r::AbstractCloseOpen) = Int
@inline Base.eachindex(r::AbstractCloseOpen) = StaticInt(1):StaticArrayInterface.static_length(r)

@static if isdefined(Base.IteratorsMD, :OrdinalRangeInt)
  @inline function Base.IteratorsMD.__inc(state::Tuple{Int,Int,Vararg{Int}}, indices::Tuple{AbstractCloseOpen,Vararg{Base.IteratorsMD.OrdinalRangeInt}})
    rng = indices[1]
    I1 = state[1] + step(rng)
    if Base.IteratorsMD.__is_valid_range(I1, rng) && state[1] != last(rng)
      return true, (I1, Base.tail(state)...)
    end
    valid, I = Base.IteratorsMD.__inc(Base.tail(state), Base.tail(indices))
    return valid, (convert(typeof(I1), first(rng)), I...)
  end
else
  @inline function Base.IteratorsMD.__inc(state::Tuple{Int,Int,Vararg{Int}}, indices::Tuple{AbstractCloseOpen,Vararg})
    rng = indices[1]
    I1 = state[1] + step(rng)
    if Base.IteratorsMD.__is_valid_range(I1, rng) && state[1] != last(rng)
      return true, (I1, Base.tail(state)...)
    end
    valid, I = Base.IteratorsMD.__inc(Base.tail(state), Base.tail(indices))
    return valid, (convert(typeof(I1), first(rng)), I...)
  end
end
end
