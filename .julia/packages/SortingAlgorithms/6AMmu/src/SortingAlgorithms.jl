__precompile__()

module SortingAlgorithms

using DataStructures
using Base.Sort
using Base.Order

import Base.Sort: sort!
import DataStructures: heapify!, percolate_down!

export HeapSort, TimSort, RadixSort, CombSort, PagedMergeSort

struct HeapSortAlg  <: Algorithm end
struct TimSortAlg   <: Algorithm end
struct RadixSortAlg <: Algorithm end
struct CombSortAlg  <: Algorithm end
struct PagedMergeSortAlg  <: Algorithm end

function maybe_optimize(x::Algorithm)
    isdefined(Base.Sort, :InitialOptimizations) ? Base.Sort.InitialOptimizations(x) : x
end
const HeapSort  = maybe_optimize(HeapSortAlg())
const TimSort   = maybe_optimize(TimSortAlg())
# Whenever InitialOptimizations is defined, RadixSort falls
# back to Base.DEFAULT_STABLE which already includes them.
const RadixSort = RadixSortAlg()

"""
    CombSort

Indicates that a sorting function should use the comb sort
algorithm. Comb sort traverses the collection multiple times
ordering pairs of elements with a given interval between them.
The interval decreases exponentially until it becomes 1, then
it switches to insertion sort on the whole input.

Characteristics:
 - *not stable* does not preserve the ordering of elements which
   compare equal (e.g. "a" and "A" in a sort of letters which
   ignores case).
 - *in-place* in memory.
 - *parallelizable* suitable for vectorization with SIMD instructions
   because it performs many independent comparisons.
 - *pathological inputs* such as `repeat(1:5.0, 4^8)` can make this algorithm perform very poorly.
 - *`n log n` average runtime* measured for random inputs of length up to 100 million, but theoretical runtime of `Θ(n^2)` for extremely long inputs.

## References
- Dobosiewicz, Wlodzimierz, (1980). "An efficient variation of bubble sort", Information Processing Letters, 11(1), pp. 5-6, https://doi.org/10.1016/0020-0190(80)90022-8.
 - Werneck, N. L., (2020). "ChipSort: a SIMD and cache-aware sorting module. JuliaCon Proceedings, 1(1), 12, https://doi.org/10.21105/jcon.00012
 - H. Inoue, T. Moriyama, H. Komatsu and T. Nakatani, "AA-Sort: A New Parallel Sorting Algorithm for Multi-Core SIMD Processors," 16th International Conference on Parallel Architecture and Compilation Techniques (PACT 2007), 2007, pp. 189-198, doi: 10.1109/PACT.2007.4336211.
"""
const CombSort  = maybe_optimize(CombSortAlg())

"""
    PagedMergeSort

Indicates that a sorting function should use the paged merge sort
algorithm. Paged merge sort uses is a merge sort, that uses different
merge routines to achieve stable sorting with a scratch space of size O(√n).
The merge routine for merging large subarrays merges
pages of size O(√n) almost in place, before reordering them using a page table.
At deeper recursion levels, where the scratch space is big enough,
normal merging is used, where one input is copied into the scratch space.
When the scratch space is large enough to hold the complete subarray,
the input is merged interleaved from both sides, which increases performance
for random data.

Characteristics:
 - *stable*: does preserve the ordering of elements which
   compare equal (e.g. "a" and "A" in a sort of letters which
   ignores case).
 - *`O(√n)`* auxilary memory usage.
 - *`O(n log n)`* garuanteed runtime.

## References
 - Dvořák, S., Ďurian, B. (1986). Towards an efficient merging. In: Gruska, J., Rovan, B., Wiedermann,
   J. (eds) Mathematical Foundations of Computer Science 1986. MFCS 1986. Lecture Notes in Computer Science, vol 233.
   Springer, Berlin, Heidelberg. https://doi.org/10.1007/BFb0016253
 - https://max-arbuzov.blogspot.com/2021/10/merge-sort-with-osqrtn-auxiliary-memory.html
"""
const PagedMergeSort  = maybe_optimize(PagedMergeSortAlg())

## Heap sort

function sort!(v::AbstractVector, lo::Int, hi::Int, a::HeapSortAlg, o::Ordering)
    if lo > 1 || hi < length(v)
        return sort!(view(v, lo:hi), 1, length(v), a, o)
    end
    r = ReverseOrdering(o)
    heapify!(v, r)
    for i = length(v):-1:2
        # Swap the root with i, the last unsorted position
        x = v[i]
        v[i] = v[1]
        # The heap portion now ends at position i-1, but needs fixing up
        # starting with the root
        percolate_down!(v,1,x,r,i-1)
    end
    v
end


# Implementation of TimSort based on the algorithm description at:
#
#   http://svn.python.org/projects/python/trunk/Objects/listsort.txt
#   http://en.wikipedia.org/wiki/Timsort
#
# Original author: @kmsquire

@static if v"1.9.0-alpha" <= VERSION <= v"1.9.1"
    function Base.getindex(v::Base.Sort.WithoutMissingVector, i::UnitRange)
        out = Vector{eltype(v)}(undef, length(i))
        out .= v.data[i]
        out
    end

    # skip MissingOptimization due to JuliaLang/julia#50171
    const _FIVE_ARG_SAFE_DEFAULT_STABLE = Base.DEFAULT_STABLE.next

    # Explicitly define conversion from _sort!(v, alg, order, kw) to sort!(v, lo, hi, alg, order)
    # To avoid excessively strict dispatch loop detection
    function Base.Sort._sort!(v::AbstractVector, a::Union{HeapSortAlg, TimSortAlg, RadixSortAlg, CombSortAlg}, o::Base.Order.Ordering, kw)
        Base.Sort.@getkw lo hi scratch
        sort!(v, lo, hi, a, o)
        scratch
    end
else
    const _FIVE_ARG_SAFE_DEFAULT_STABLE = Base.DEFAULT_STABLE
end

const Run = UnitRange{Int}

const MIN_GALLOP = 7

mutable struct MergeState
    runs::Vector{Run}
    min_gallop::Int
end
MergeState() = MergeState(Run[], MIN_GALLOP)

# Determine a good minimum run size for efficient merging
# For details, see "Computing minrun" in
# http://svn.python.org/projects/python/trunk/Objects/listsort.txt
function merge_compute_minrun(N::Int, bits::Int)
    r = 0
    max_val = 2^bits
    while N >= max_val
        r |= (N & 1)
        N >>= 1
    end
    N + r
end
merge_compute_minrun(N::Int) = merge_compute_minrun(N, 6)

# Galloping binary search starting at left
# Finds the last value in v <= x
function gallop_last(o::Ordering, v::AbstractVector, x, lo::Int, hi::Int)
    i = lo
    inc = 1
    lo = lo-1
    hi = hi+1
    while i < hi && !lt(o, x, v[i])
        lo = i
        i += inc
        inc <<= 1
    end
    hi = min(i+1, hi)
    # Binary search
    while lo < hi-1
        i = (lo+hi)>>>1
        if lt(o, x, v[i])
            hi = i
        else
            lo = i
        end
    end
    lo
end

# Galloping binary search starting at right
# Finds the last value in v <= x
function rgallop_last(o::Ordering, v::AbstractVector, x, lo::Int, hi::Int)
    i = hi
    dec = 1
    lo = lo-1
    hi = hi+1
    while i > lo && lt(o, x, v[i])
        hi = i
        i -= dec
        dec <<= 1
    end
    lo = max(lo, i-1)
    # Binary search
    while lo < hi-1
        i = (lo+hi)>>>1
        if lt(o, x, v[i])
            hi = i
        else
            lo = i
        end
    end
    lo
end

# Galloping binary search starting at left
# Finds the first value in v >= x
function gallop_first(o::Ordering, v::AbstractVector, x, lo::Int, hi::Int)
    i = lo
    inc = 1
    lo = lo-1
    hi = hi+1
    while i < hi && lt(o, v[i], x)
        lo = i
        i += inc
        inc <<= 1
    end
    hi = min(i+1, hi)
    # Binary search
    while lo < hi-1
        i = (lo+hi)>>>1
        if lt(o, v[i], x)
            lo = i
        else
            hi = i
        end
    end
    hi
end

# Galloping binary search starting at right
# Finds the first value in v >= x
function rgallop_first(o::Ordering, v::AbstractVector, x, lo::Int, hi::Int)
    i = hi
    dec = 1
    lo = lo-1
    hi = hi+1
    while i > lo && !lt(o, v[i], x)
        hi = i
        i -= dec
        dec <<= 1
    end
    lo = max(lo, i-1)
    # Binary search
    while lo < hi-1
        i = (lo+hi)>>>1
        if lt(o, v[i], x)
            lo = i
        else
            hi = i
        end
    end
    hi
end

# Get the next run
# Returns the v range a:b, or b:-1:a for a reversed sequence
function next_run(o::Ordering, v::AbstractVector, lo::Int, hi::Int)
    lo == hi && return lo:hi
    if !lt(o, v[lo+1], v[lo])
        for i = lo+2:hi
            if lt(o, v[i], v[i-1])
                return lo:i-1
            end
        end
        return lo:hi
    else
        for i = lo+2:hi
            if !lt(o, v[i], v[i-1])
                return i-1:-1:lo
            end
        end
        return hi:-1:lo
    end
end

function merge_at(o::Ordering, v::AbstractVector, state::MergeState, n::Integer)
    a = state.runs[n]
    b = state.runs[n+1]
    merge(o,v,a,b,state)
    state.runs[n] = first(a):last(b)
    deleteat!(state.runs, n+1)
    nothing
end

# Merge consecutive runs
# For A,B,C,D = last four lengths, merge_collapse!()
# maintains 3 invariants:
#
#  A > B + C
#  B > C + D
#  C > D
#
# If any of these are violated, a merge occurs to
# correct it
function merge_collapse(o::Ordering, v::AbstractVector, state::MergeState)
    while true
        n = length(state.runs)
        n <= 1 && break

        # Check invariants 1 and 2
        if (n >= 3 && length(state.runs[end-2]) <= length(state.runs[end-1]) + length(state.runs[end])) ||
            (n >= 4 && length(state.runs[end-3]) <= length(state.runs[end-2]) + length(state.runs[end-1]))

            if length(state.runs[end-2]) < length(state.runs[end])
                merge_at(o,v,state,n-2)
            else
                merge_at(o,v,state,n-1)
            end

        # Check invariant 3
        elseif length(state.runs[end-1]) <= length(state.runs[end])
            merge_at(o,v,state,n-1)

        else # Invariant is satisfied
            break
        end
    end
end

# Merge runs a and b in vector v
function merge(o::Ordering, v::AbstractVector, a::Run, b::Run, state::MergeState)

    # First elements in a <= b[1] are already in place
    a = gallop_last(o, v, v[first(b)], first(a), last(a))+1: last(a)

    if length(a) == 0  return  end

    # Last elements in b >= a[end] are already in place
    b = first(b) : rgallop_first(o, v, v[last(a)], first(b), last(b))-1

    # Choose merge_lo or merge_hi based on the amount
    # of temporary memory needed (smaller of a and b)
    if length(a) < length(b)
        merge_lo(o, v, a, b, state)
    else
        merge_hi(o, v, a, b, state)
    end
end

# Merge runs a and b in vector v (a is smaller)
function merge_lo(o::Ordering, v::AbstractVector, a::Run, b::Run, state::MergeState)

    # Copy a
    v_a = v[a]

    # Pointer into (sub)arrays
    i = first(a)
    from_a = 1
    from_b = first(b)

    mode = :normal
    while true
        if mode == :normal
            # Compare and copy element by element
            count_a = count_b = 0
            while from_a <= length(a) && from_b <= last(b)
                if lt(o, v[from_b], v_a[from_a])
                    v[i] = v[from_b]
                    from_b += 1
                    count_a = 0
                    count_b += 1
                else
                    v[i] = v_a[from_a]
                    from_a += 1
                    count_a += 1
                    count_b = 0
                end
                i += 1

                # Switch to galloping mode if a string of elements
                # has come from the same set
                if count_b >= state.min_gallop || count_a >= state.min_gallop
                    mode = :galloping
                    break
                end
            end
            # Finalize if we've exited the loop normally
            if mode == :normal
                mode = :finalize
            end
        end

        if mode == :galloping
            # Use binary search to find range to copy
            while from_a <= length(a) && from_b <= last(b)
                # Copy the next run from b
                b_run = from_b : gallop_first(o, v, v_a[from_a], from_b, last(b)) - 1
                i_end = i + length(b_run) - 1
                v[i:i_end] = v[b_run]
                i = i_end + 1
                from_b = last(b_run) + 1

                # ... then copy the first element from a
                v[i] = v_a[from_a]
                i += 1
                from_a += 1

                if from_a > length(a) || from_b > last(b) break end

                # Copy the next run from a
                a_run = from_a : gallop_last(o, v_a, v[from_b], from_a, length(a))
                i_end = i + length(a_run) - 1
                v[i:i_end] = v_a[a_run]
                i = i_end + 1
                from_a = last(a_run) + 1

                # ... then copy the first element from b
                v[i] = v[from_b]
                i += 1
                from_b += 1

                # Return to normal mode if we haven't galloped...
                if length(a_run) < MIN_GALLOP && length(b_run) < MIN_GALLOP
                    mode = :normal
                    break
                end
                # Reduce min_gallop if this gallop was successful
                state.min_gallop -= 1
            end
            if mode == :galloping
                mode = :finalize
            end
            state.min_gallop = max(state.min_gallop,0) + 2  # penalty for leaving gallop mode
        end

        if mode == :finalize
            # copy end of a
            i_end = i + (length(a) - from_a)
            v[i:i_end] = v_a[from_a:end]
            break
        end
    end
end

# Merge runs a and b in vector v (b is smaller)
function merge_hi(o::Ordering, v::AbstractVector, a::Run, b::Run, state::MergeState)

    # Copy b
    v_b = v[b]

    # Pointer into (sub)arrays
    i = last(b)
    from_a = last(a)
    from_b = length(b)

    mode = :normal
    while true
        if mode == :normal
            # Compare and copy element by element
            count_a = count_b = 0
            while from_a >= first(a) && from_b >= 1
                if !lt(o, v_b[from_b], v[from_a])
                    v[i] = v_b[from_b]
                    from_b -= 1
                    count_a = 0
                    count_b += 1
                else
                    v[i] = v[from_a]
                    from_a -= 1
                    count_a += 1
                    count_b = 0
                end
                i -= 1

                # Switch to galloping mode if a string of elements
                # has come from the same set
                if count_b >= state.min_gallop || count_a >= state.min_gallop
                   mode = :galloping
                   break
                end
            end
            # Finalize if we've exited the loop normally
            if mode == :normal
                mode = :finalize
            end
        end

        if mode == :galloping
            # Use binary search to find range to copy
            while from_a >= first(a) && from_b >= 1
                # Copy the next run from b
                b_run = rgallop_first(o, v_b, v[from_a], 1, from_b) : from_b
                i_start = i - length(b_run) + 1
                v[i_start:i] = v_b[b_run]
                i = i_start - 1
                from_b = first(b_run) - 1

                # ... then copy the first element from a
                v[i] = v[from_a]
                i -= 1
                from_a -= 1

                if from_a < first(a) || from_b < 1 break end

                # Copy the next run from a
                a_run = rgallop_last(o, v, v_b[from_b], first(a), from_a) + 1: from_a
                i_start = i - length(a_run) + 1
                v[i_start:i] = v[a_run]
                i = i_start - 1
                from_a = first(a_run) - 1

                # ... then copy the first element from b
                v[i] = v_b[from_b]
                i -= 1
                from_b -= 1

                # Return to normal mode if we haven't galloped...
                if length(a_run) < MIN_GALLOP && length(b_run) < MIN_GALLOP
                    mode = :normal
                    break
                end
                # Reduce min_gallop if this gallop was successful
                state.min_gallop -= 1
            end
            if mode == :galloping
                mode = :finalize
            end
            state.min_gallop = max(state.min_gallop, 0) + 2  # penalty for leaving gallop mode
        end

        if mode == :finalize
            # copy start of b
            i_start = i - from_b + 1
            v[i_start:i] = v_b[1:from_b]
            break
        end
    end
end

# TimSort main method
function sort!(v::AbstractVector, lo::Int, hi::Int, ::TimSortAlg, o::Ordering)
    minrun = merge_compute_minrun(hi-lo+1)
    state = MergeState()
    i = lo
    while i <= hi
        run_range = next_run(o, v, i, hi)
        count = length(run_range)
        if count < minrun
            # Make a run of length minrun
            count = min(minrun, hi-i+1)
            run_range = i:i+count-1
            sort!(v, i, i+count-1, _FIVE_ARG_SAFE_DEFAULT_STABLE, o)
        else
            if !issorted(run_range)
                run_range = last(run_range):first(run_range)
                reverse!(view(v, run_range))
            end
        end
        # Push this run onto the queue and merge if needed
        push!(state.runs, run_range)
        i = i+count
        merge_collapse(o, v, state)
    end
    # Force merge at the end
    while true
        n = length(state.runs)
        n <= 1 && break
        merge_at(o, v, state, n-1)
    end
    return v
end


function sort!(v::AbstractVector, lo::Int, hi::Int, ::CombSortAlg, o::Ordering)
    interval = (3 * (hi-lo+1)) >> 2

    while interval > 1
        @inbounds for j in lo:hi-interval
            a, b = v[j], v[j+interval]
            v[j], v[j+interval] = lt(o, b, a) ? (b, a) : (a, b)
        end
        interval = (3 * interval) >> 2
    end

    return sort!(v, lo, hi, InsertionSort, o)
end


## Radix sort
@static if VERSION >= v"1.9.0-DEV.482" # Base introduced radixsort in 1.9
    function sort!(vs::AbstractVector{T}, lo::Int, hi::Int, ::RadixSortAlg, o::Ordering, ts::Union{Nothing, AbstractVector{T}}=nothing) where T
        sort!(vs, lo, hi, Base.DEFAULT_STABLE, o)
    end
    function sort!(vs::AbstractVector{T}, lo::Int, hi::Int, ::RadixSortAlg, o::Ordering, ts::Vector{T}) where T # disambiguate
        sort!(vs, lo, hi, Base.DEFAULT_STABLE, o)
    end
else

    # Map a bits-type to an unsigned int, maintaining sort order
    uint_mapping(::ForwardOrdering, x::Unsigned) = x
    for (signedty, unsignedty) in ((Int8, UInt8), (Int16, UInt16), (Int32, UInt32), (Int64, UInt64), (Int128, UInt128))
        # In Julia 0.4 we can just use unsigned() here
        @eval uint_mapping(::ForwardOrdering, x::$signedty) = reinterpret($unsignedty, xor(x, typemin(typeof(x))))
    end
    uint_mapping(::ForwardOrdering, x::Float32)  = (y = reinterpret(Int32, x); reinterpret(UInt32, ifelse(y < 0, ~y, xor(y, typemin(Int32)))))
    uint_mapping(::ForwardOrdering, x::Float64)  = (y = reinterpret(Int64, x); reinterpret(UInt64, ifelse(y < 0, ~y, xor(y, typemin(Int64)))))

    uint_mapping(::Sort.Float.Left, x::Float16)  = ~reinterpret(Int16, x)
    uint_mapping(::Sort.Float.Right, x::Float16)  = reinterpret(Int16, x)
    uint_mapping(::Sort.Float.Left, x::Float32)  = ~reinterpret(Int32, x)
    uint_mapping(::Sort.Float.Right, x::Float32)  = reinterpret(Int32, x)
    uint_mapping(::Sort.Float.Left, x::Float64)  = ~reinterpret(Int64, x)
    uint_mapping(::Sort.Float.Right, x::Float64)  = reinterpret(Int64, x)

    uint_mapping(rev::ReverseOrdering, x) = ~uint_mapping(rev.fwd, x)
    uint_mapping(::ReverseOrdering{ForwardOrdering}, x::Real) = ~uint_mapping(Forward, x) # maybe unnecessary; needs benchmark

    uint_mapping(o::By,   x     ) = uint_mapping(Forward, o.by(x))
    uint_mapping(o::Perm, i::Int) = uint_mapping(o.order, o.data[i])
    uint_mapping(o::Lt,   x     ) = error("uint_mapping does not work with general Lt Orderings")

    const RADIX_SIZE = 11
    const RADIX_MASK = 0x7FF

    function sort!(vs::AbstractVector{T}, lo::Int, hi::Int, ::RadixSortAlg, o::Ordering, ts::AbstractVector{T}=similar(vs)) where T
        # Input checking
        if lo >= hi;  return vs;  end

        # Make sure we're sorting a bits type
        OT = Base.Order.ordtype(o, vs)
        if !isbitstype(OT)
            error("Radix sort only sorts bits types (got $OT)")
        end

        # Init
        iters = ceil(Integer, sizeof(OT)*8/RADIX_SIZE)
        bin = zeros(UInt32, 2^RADIX_SIZE, iters)
        if lo > 1;  bin[1,:] .= lo-1;  end

        # Histogram for each element, radix
        for i = lo:hi
            v = uint_mapping(o, vs[i])
            for j = 1:iters
                idx = Int((v >> ((j-1)*RADIX_SIZE)) & RADIX_MASK) + 1
                @inbounds bin[idx,j] += 1
            end
        end

        # Sort!
        swaps = 0
        len = hi-lo+1
        for j = 1:iters
            # Unroll first data iteration, check for degenerate case
            v = uint_mapping(o, vs[hi])
            idx = Int((v >> ((j-1)*RADIX_SIZE)) & RADIX_MASK) + 1

            # are all values the same at this radix?
            if bin[idx,j] == len;  continue;  end

            cbin = cumsum(bin[:,j])
            ci = cbin[idx]
            ts[ci] = vs[hi]
            cbin[idx] -= 1

            # Finish the loop...
            @inbounds for i in hi-1:-1:lo
                v = uint_mapping(o, vs[i])
                idx = Int((v >> ((j-1)*RADIX_SIZE)) & RADIX_MASK) + 1
                ci = cbin[idx]
                ts[ci] = vs[i]
                cbin[idx] -= 1
            end
            vs,ts = ts,vs
            swaps += 1
        end

        if isodd(swaps)
            vs,ts = ts,vs
            for i = lo:hi
                vs[i] = ts[i]
            end
        end
        vs
    end

    function Base.Sort.Float.fpsort!(v::AbstractVector, ::RadixSortAlg, o::Ordering)
        @static if VERSION >= v"1.7.0-DEV"
            lo, hi = Base.Sort.Float.specials2end!(v, RadixSort, o)
        else
            lo, hi = Base.Sort.Float.nans2end!(v, o)
        end
        sort!(v, lo, hi, RadixSort, o)
    end
end

###
# PagedMergeSort
###

# unsafe version of copyto!
# as workaround for https://github.com/JuliaLang/julia/issues/50900
function _unsafe_copyto!(dest, doffs, src, soffs, n)
    @inbounds for i in 0:n-1
        dest[doffs + i] = src[soffs + i]
    end
    dest
end

function _unsafe_copyto!(dest::Array, doffs, src::Array, soffs, n)
    unsafe_copyto!(dest, doffs, src, soffs, n)
end

# merge v[lo:m] and v[m+1:hi] ([A;B]) using scratch[1:1+hi-lo]
# This is faster than merge! but requires twice as much auxiliary memory.
function twoended_merge!(v::AbstractVector{T}, lo::Integer, m::Integer, hi::Integer, o::Ordering, scratch::AbstractVector{T}) where T
    @assert lo ≤ m ≤ hi
    @assert abs((m-lo) - (hi-(m+1))) ≤ 1 "twoended_merge! only supports balanced merges"
    len = 1 + hi - lo
    # input array indices
    a_lo = lo
    a_hi = m
    b_lo = m + 1
    b_hi = hi
    # output array indices
    k_lo = 1
    k_hi = len
    @inbounds begin
        # two ended merge
        while k_lo <= len ÷ 2
            if lt(o, v[b_lo], v[a_lo])
                scratch[k_lo] = v[b_lo]
                b_lo += 1
            else
                scratch[k_lo] = v[a_lo]
                a_lo += 1
            end
            k_lo +=1
            if !lt(o, v[b_hi], v[a_hi])
                scratch[k_hi] = v[b_hi]
                b_hi -= 1
            else
                scratch[k_hi] = v[a_hi]
                a_hi -= 1
            end
            k_hi -=1
        end
        # if the input length is odd,
        # one item remains
        if a_lo <= a_hi
            scratch[k_lo] = v[a_lo]
        elseif b_lo <= b_hi
            scratch[k_lo] = v[b_lo]
        end
        # copy back from t to v
        offset = lo-1
        for i = 1:len
            v[offset+i] = scratch[i]
        end
    end
end

# core merging loop used throughout PagedMergeSort
Base.@propagate_inbounds function merge!(f::Function,
    target::AbstractVector{T}, source_a::AbstractVector{T}, source_b::AbstractVector{T},
    o::Ordering, a::Integer, b::Integer, k::Integer) where T
    @inbounds while f(a,b,k)
        if lt(o, source_b[b], source_a[a])
            target[k] = source_b[b]
            b += 1
        else
            target[k] = source_a[a]
            a += 1
        end
        k += 1
    end
    a,b,k
end

# merge v[lo:m] and v[m+1:hi] using scratch[1:1+m-lo]
# based on Base.Sort MergeSort
function merge!(v::AbstractVector{T}, lo::Integer, m::Integer, hi::Integer, o::Ordering, scratch::AbstractVector{T}) where {T}
    _unsafe_copyto!(scratch, 1, v, lo, m - lo + 1)
    f(_, b, k) = k < b <= hi
    a, b, k = merge!(f, v, scratch, v, o, 1, m + 1, lo)
    _unsafe_copyto!(v, k, scratch, a, b - k)
end

struct Pages
    current::Int       # current page being merged into
    currentNumber::Int # number of current page (=index in pageLocations)
    nextA::Int         # next possible page in A
    nextB::Int         # next possible page in B
end

next_page_A(pages::Pages) = Pages(pages.nextA, pages.currentNumber + 1, pages.nextA + 1, pages.nextB)
next_page_B(pages::Pages) = Pages(pages.nextB, pages.currentNumber + 1, pages.nextA, pages.nextB + 1)

Base.@propagate_inbounds function next_page!(pageLocations, pages, pagesize, lo, a)
    if a > pages.nextA * pagesize + lo
        pages = next_page_A(pages)
    else
        pages = next_page_B(pages)
    end
    pageLocations[pages.currentNumber] = pages.current
    pages
end

Base.@propagate_inbounds function permute_pages!(f, v, pageLocations, page_offset, pagesize, page)
    while f(page)
        plc = pageLocations[page-3] # plc has data belonging to page
        pageLocations[page-3] = page
        _unsafe_copyto!(v, page_offset(page) + 1, v, page_offset(plc) + 1, pagesize)
        page = plc
    end
    page
end

# merge v[lo:m] (A) and v[m+1:hi] (B) using scratch[] in O(sqrt(n)) space
function paged_merge!(v::AbstractVector{T}, lo::Integer, m::Integer, hi::Integer, o::Ordering, scratch::AbstractVector{T}, pageLocations::AbstractVector{<:Integer}) where {T}
    @assert lo < m < hi
    lenA = 1 + m - lo
    lenB = hi - m

    # this function only supports merges with length(A) <= length(B),
    # which is guaranteed by pagedmergesort!
    @assert lenA <= lenB

    # regular merge if scratch is big enough
    lenA <= length(scratch) && return merge!(v, lo, m, hi, o, scratch)

    len = lenA + lenB
    pagesize = isqrt(len)
    nPages = len ÷ pagesize # a partial page at the end does not count
    @assert length(scratch) >= 3pagesize
    @assert length(pageLocations) >= nPages - 3

    @inline page_offset(page) = (page - 1) * pagesize + lo - 1

    @inbounds begin
        ##################
        # merge
        ##################
        # merge the first 3 pages into scratch
        a, b, _ = merge!((_, _, k) -> k <= 3pagesize, scratch, v, v, o, lo, m + 1, 1)
        # initialize variables for merging into pages
        pages = Pages(-17, 0, 1, (m - lo) ÷ pagesize + 2) # first argument is unused
        # more efficient loop while more than pagesize elements of A and B are remaining
        while_condition1(offset) = (_, _, k) -> k <= offset + pagesize
        while a < m - pagesize && b < hi - pagesize
            pages = next_page!(pageLocations, pages, pagesize, lo, a)
            offset = page_offset(pages.current)
            a, b, _ = merge!(while_condition1(offset), v, v, v, o, a, b, offset + 1)
        end
        # merge until either A or B is empty or the last page is reached
        k, offset = nothing, nothing
        while_condition2(offset) = (a, b, k) -> k <= offset + pagesize && a <= m && b <= hi
        while a <= m && b <= hi && pages.currentNumber + 3 < nPages
            pages = next_page!(pageLocations, pages, pagesize, lo, a)
            offset = page_offset(pages.current)
            a, b, k = merge!(while_condition2(offset), v, v, v, o, a, b, offset + 1)
        end
        # if the last page is reached, merge the remaining elements into the final partial page
        if pages.currentNumber + 3 == nPages && a <= m && b <= hi
            a, b, k = merge!((a, b, _) -> a <= m && b <= hi, v, v, v, o, a, b, nPages * pagesize + lo)
            _unsafe_copyto!(v, k, v, a <= m ? a : b, hi - k + 1)
        else
            use_a = a <= m
            # copy the incomplete page
            partial_page_size = offset + pagesize - k + 1
            _unsafe_copyto!(v, k, v, use_a ? a : b, partial_page_size)
            use_a && (a += partial_page_size)
            use_a || (b += partial_page_size)
            # copy the remaining full pages
            while use_a ? a <= m - pagesize + 1 : b <= hi - pagesize + 1
                pages = next_page!(pageLocations, pages, pagesize, lo, a)
                offset = page_offset(pages.current)
                _unsafe_copyto!(v, offset + 1, v, use_a ? a : b, pagesize)
                use_a && (a += pagesize)
                use_a || (b += pagesize)
            end
            # copy the final partial page only if sourcing from A.
            # If sourcing from B, it is already in place.
            use_a && _unsafe_copyto!(v, hi - m + a, v, a, m - a + 1)
        end

        ##################
        # rearrange pages
        ##################
        # copy pages belonging to the 3 permutation chains ending with a page in the scratch space
        nextA, nextB = pages.nextA, pages.nextB

        for _ in 1:3
            page = (nextB > nPages ? (nextA += 1) : (nextB += 1)) - 1
            page = permute_pages!(>(3), v, pageLocations, page_offset, pagesize, page)
            _unsafe_copyto!(v, page_offset(page) + 1, scratch, (page - 1) * pagesize + 1, pagesize)
        end

        # copy remaining permutation cycles
        for donePageIndex = 5:nPages
            # linear scan through pageLocations to make sure no cycle is missed
            page = pageLocations[donePageIndex-3]
            page == donePageIndex && continue

            # copy the data belonging to donePageIndex into scratch
            _unsafe_copyto!(scratch, 1, v, page_offset(page) + 1, pagesize)

            # follow the cycle starting with the newly freed page
            permute_pages!(!=(donePageIndex), v, pageLocations, page_offset, pagesize, page)
            _unsafe_copyto!(v, page_offset(donePageIndex) + 1, scratch, 1, pagesize)
        end
    end
end

# midpoint was added to Base.sort in version 1.4 and later moved to Base
# -> redefine for compatibility with earlier versions
midpoint(lo::Integer, hi::Integer) = lo + ((hi - lo) >>> 0x01)

function pagedmergesort!(v::AbstractVector{T}, lo::Integer, hi::Integer, o::Ordering, scratch::AbstractVector{T}, pageLocations) where {T}
    len = hi + 1 - lo
    if len <= Base.SMALL_THRESHOLD
        return Base.Sort.sort!(v, lo, hi, Base.Sort.InsertionSortAlg(), o)
    end
    m = midpoint(lo, hi - 1) # hi-1: ensure midpoint is rounded down. OK, because lo < hi is satisfied here
    pagedmergesort!(v, lo, m, o, scratch, pageLocations)
    pagedmergesort!(v, m + 1, hi, o, scratch, pageLocations)
    if len <= length(scratch)
        twoended_merge!(v, lo, m, hi, o, scratch)
    else
        paged_merge!(v, lo, m, hi, o, scratch, pageLocations)
    end
    return v
end

function sort!(v::AbstractVector, lo::Integer, hi::Integer, ::PagedMergeSortAlg, o::Ordering)
    lo >= hi && return v
    n = hi + 1 - lo
    pagesize = isqrt(n)
    scratch = Vector{eltype(v)}(undef, 3pagesize)
    nPages = n ÷ pagesize
    pageLocations = Vector{Int}(undef, max(0, nPages - 3))
    pagedmergesort!(v, lo, hi, o, scratch, pageLocations)
    return v
end
end # module
