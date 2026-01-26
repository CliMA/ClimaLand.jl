import CPUSummary

function worker_bits()
  wts = nextpow2(CPUSummary.sys_threads()) # Typically sys_threads (i.e. Sys.CPU_THREADS) does not change between runs, thus it will precompile well.
  ws = static(8sizeof(UInt))               # For testing purposes it can be overridden by JULIA_CPU_THREADS,
  ifelse(Static.lt(wts, ws), ws, wts)
end
function worker_mask_count()
  bits = worker_bits()
  (bits + StaticInt{63}()) ÷ StaticInt{64}() # cld not defined on `StaticInt`
end

worker_pointer() = Base.unsafe_convert(Ptr{UInt}, pointer_from_objref(WORKERS))

function free_threads!(freed_threads::U) where {U<:Unsigned}
  _atomic_or!(worker_pointer(), freed_threads)
  nothing
end
function free_threads!(freed_threads_tuple::NTuple{1,U}) where {U<:Unsigned}
  _atomic_or!(worker_pointer(), freed_threads_tuple[1])
  nothing
end
function free_threads!(
  freed_threads_tuple::Tuple{U,Vararg{U,N}}
) where {N,U<:Unsigned}
  wp = worker_pointer()
  for freed_threads in freed_threads_tuple
    _atomic_or!(wp, freed_threads)
    wp += sizeof(UInt)
  end
  nothing
end

@inline _remaining(x::Tuple) = Base.tail(x)
@inline _remaining(@nospecialize(x)) = nothing
@inline _first(::Tuple{}) = nothing
@inline _first(x::Tuple{X,Vararg}) where {X<:Unsigned} = getfield(x, 1)
@inline _first(x::Union{Unsigned,Nothing}) = x
@inline function _request_threads(
  num_requested::UInt32,
  wp::Ptr,
  ::StaticInt{N},
  threadmask
) where {N}
  ui, ft, num_requested, wp =
    __request_threads(num_requested, wp, _first(threadmask))
  uit, ftt = _request_threads(
    num_requested,
    wp,
    StaticInt{N}() - StaticInt{1}(),
    _remaining(threadmask)
  )
  (ui, uit...), (ft, ftt...)
end
@inline function _request_threads(
  num_requested::UInt32,
  wp::Ptr,
  ::StaticInt{1},
  threadmask
)
  ui, ft, num_requested, wp =
    __request_threads(num_requested, wp, _first(threadmask))
  (ui,), (ft,)
end

@inline function _exchange_mask!(wp, ::Nothing)
  all_threads = _atomic_xchg!(wp, zero(UInt))
  all_threads, all_threads
end
@inline function _exchange_mask!(wp, threadmask::Unsigned)
  all_threads = _atomic_xchg!(wp, zero(UInt))
  tm = threadmask % UInt
  saved = all_threads & (~tm)
  _atomic_store!(wp, saved)
  all_threads | saved, all_threads & tm
end
@inline function __request_threads(num_requested::UInt32, wp::Ptr, threadmask)
  no_threads = zero(UInt)
  if (num_requested ≢ StaticInt{-1}()) && (num_requested % Int32 ≤ zero(Int32))
    return UnsignedIteratorEarlyStop(zero(UInt), 0x00000000),
    no_threads,
    0x00000000,
    wp
  end
  # to get more, we xchng, setting all to `0`
  # then see which we need, and free those we aren't using.
  wpret = wp + 8 # (UInt === UInt64) | (worker_mask_count() === StaticInt(1)) #, so adding 8 is fine.
  # _all_threads = all_threads = _apply_mask(_atomic_xchg!(wp, no_threads), threadmask)
  _all_threads, all_threads = _exchange_mask!(wp, threadmask)
  additional_threads = count_ones(all_threads) % UInt32
  # num_requested === StaticInt{-1}() && return reserved_threads, all_threads
  if num_requested === StaticInt{-1}()
    return UnsignedIteratorEarlyStop(all_threads),
    all_threads,
    num_requested,
    wpret
  end
  nexcess = num_requested - additional_threads
  if signed(nexcess) ≥ 0
    return UnsignedIteratorEarlyStop(all_threads), all_threads, nexcess, wpret
  end
  # we need to return the `excess` to the pool.
  lz = leading_zeros(all_threads) % UInt32
  while true
    # start by trying to trim off excess from lz
    lz += (-nexcess) % UInt32
    m = (one(UInt) << (UInt32(8sizeof(UInt)) - lz)) - one(UInt)
    masked = (all_threads & m) ⊻ all_threads
    nexcess += count_ones(masked) % UInt32
    all_threads &= (~masked)
    nexcess == zero(nexcess) && break
  end
  _atomic_store!(wp, _all_threads & (~all_threads))
  return UnsignedIteratorEarlyStop(all_threads, num_requested),
  all_threads,
  0x00000000,
  wpret
end

@inline function request_threads(num_requested, threadmask)
  _request_threads(
    num_requested % UInt32,
    worker_pointer(),
    worker_mask_count(),
    threadmask
  )
end
@inline request_threads(num_requested) = request_threads(num_requested, nothing)

"""
`request_threads` is called to request a set of `ThreadingUtilities` worker tasks to run on.

```
julia> using PolyesterWeave

julia> PolyesterWeave.request_threads(12)
((Thread (12) Iterator: U[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],), (0x0000000000000fff,))
```

This returns an iterator over the available threads 1-12
`U[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]`, as well as a mask that must be used for freeing
these threads when done.

Note that the thread doing the requesting should also always do part of the work. Thus if
you want to run on 13 threads, you'd request 12. If you want to run on 8 threads, you'd
request 7. Etc.

If you request more threads than available on the system, you get only what is available:

```
julia> Sys.CPU_THREADS, Threads.nthreads()
(16, 8)

julia> r1 = PolyesterWeave.request_threads(30) # Note you get only 7, not 8
((Thread (7) Iterator: U[1, 2, 3, 4, 5, 6, 7],), (0x000000000000007f,))

julia> r2 = PolyesterWeave.request_threads(30)
((Thread (0) Iterator: UInt64[],), (0x0000000000000000,))
```

The following are unsupported features that might be yanked out at any time:

During request you can also mask out which threads you can get. For instance, this call:

```
julia> Sys.CPU_THREADS, Threads.nthreads()
(96, 96)

julia> r1 = PolyesterWeave.request_threads(96, 0xa)
((Thread (2) Iterator: U[2, 4], Thread (31) Iterator: U[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]), (0x000000000000000a, 0x000000007fffffff))
```

Is saying you want 96 threads, but you have the mask `0x0a` which turns off all threads from
the first set of 64 except for 2 and 4:

```
julia> bitstring(0xa)
"00001010"
```

Hence you got the iterator `U[2,4]`, plus an iterator over all the available works from the
next batch of 64.
"""
function request_threads end
