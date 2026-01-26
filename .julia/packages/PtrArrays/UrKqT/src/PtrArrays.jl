module PtrArrays

export malloc, free, PtrArray

"""
    PtrArray(ptr::Ptr{T}, dims::Int...; check_dims=true) <: DenseArray{T}

Wrap a pointer in an `DenseArray` interface conformant `PtrArray` using the standard
Julia memory order.

Validates that `dims` are non-negative and don't overflow when multiplied if `check_dims` is
true. Weird things might happen if you set `check_dims=false` and use negative or
overflowing `dims`.

!!! note
    The Julia garbage collector is not able to track `Ptr`s, so the user is responsible for
    ensuring that the memory pointed to by `ptr` through `ptr + prod(dims) - 1` remains
    allocated throughout the time that the `PtrArray` is used and that mutable objects
    stored in a `PtrArray` are prevented from garbage collection.

see also [`malloc`](@ref), [`free`](@ref)
"""
struct PtrArray{T, N} <: DenseArray{T, N}
    ptr::Ptr{T}
    size::NTuple{N, Int}
    function PtrArray(ptr::Ptr{T}, dims::Vararg{Int, N}; check_dims=true) where {T, N}
        check_dims && checked_dims(sizeof(T), dims...; message=:PtrArray)
        new{T, N}(ptr, dims)
    end
end

# Because Core.checked_dims is buggy 😢
checked_dims(elsize::Int; message) = elsize
function checked_dims(elsize::Int, d0::Int, d::Int...; message)
    overflow = false
    neg = (d0+1) < 1
    zero = false # of d0==0 we won't have overflow since we go left to right
    len = d0
    for di in d
        len, o = Base.mul_with_overflow(len, di)
        zero |= di === 0
        overflow |= o
        neg |= (di+1) < 1
    end
    len, o = Base.mul_with_overflow(len, elsize)
    err = o | neg | overflow & !zero
    err && throw(ArgumentError("invalid $message dimensions"))
    len
end

"""
    malloc(T::Type, dims::Int...) -> PtrArray{T, N} <: AbstractArray{T, N}

Allocate a new array of type `T` and dimensions `dims` using the C stdlib's `malloc`.

`T` must be an `isbitstype`.

This array is not tracked by Julia's garbage collector, so it is the user's responsibility
to call [`free`](@ref) on it when it is no longer needed.
"""
function malloc(::Type{T}, dims::Int...) where T
    isbitstype(T) || throw(ArgumentError("malloc only supports isbits types"))
    ptr = Libc.malloc(checked_dims(sizeof(T), dims...; message=:malloc))
    ptr === C_NULL && throw(OutOfMemoryError())
    PtrArray(Ptr{T}(ptr), dims..., check_dims=false)
end

"""
    free(p::PtrArray)

Free the memory allocated by a [`PtrArray`](@ref) allocated by [`malloc`](@ref).

It is only safe to call this function on `PtrArray`s returned by `malloc`, and it is unsafe
to perform any operation on a `PtrArray` after calling `free`.
"""
free(p::PtrArray) = Libc.free(p.ptr)

Base.size(p::PtrArray) = p.size
Base.IndexStyle(::Type{<:PtrArray}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(p::PtrArray, i::Int)
    @boundscheck checkbounds(p, i)
    unsafe_load(p.ptr, i)
end
Base.@propagate_inbounds function Base.setindex!(p::PtrArray, v, i::Int)
    @boundscheck checkbounds(p, i)
    unsafe_store!(p.ptr, v, i)
    p
end

# Strided array interface https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-strided-arrays
Base.unsafe_convert(::Type{Ptr{T}}, p::PtrArray{T}) where T = p.ptr
Base.elsize(::Type{P}) where P<:PtrArray = sizeof(eltype(P))

end
