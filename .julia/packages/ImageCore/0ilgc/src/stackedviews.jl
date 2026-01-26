struct StackedView{T<:Number,N,A<:Tuple{Vararg{AbstractArray{T}}}} <: AbstractArray{T,N}
    parents::A

    function StackedView{T,N,A}(parents::A) where {T<:Number,N,A<:Tuple{Vararg{AbstractArray{T}}}}
        inds = axes(parents[1])
        length(inds) == N-1 || throw(DimensionMismatch("component arrays must be of dimension \$(N-1), got \$(length(inds))"))
        for i = 2:length(parents)
            axes(parents[i]) == inds || throw(DimensionMismatch("all arrays must have the same indices, got \$inds and \$(axes(parents[i]))"))
        end
        new(parents)
    end
end

"""
    StackedView(B, C, ...) -> A

Present arrays `B`, `C`, etc, as if they are separate channels along
the first dimension of `A`. In particular,

    B == A[1,:,:...]
    C == A[2,:,:...]

and so on. Combined with `colorview`, this allows one to combine two
or more grayscale images into a single color image.

See also: [`colorview`](@ref).
"""
StackedView

@inline Base.size(V::StackedView) = (length(V.parents), size(V.parents[1])...)
@inline Base.axes(V::StackedView) = (Base.OneTo(length(V.parents)), axes(V.parents[1])...)

@inline function Base.getindex(V::StackedView{T,N}, I::Vararg{Int,N}) where {T,N}
    i1, itail = I[1], tail(I)
    P = V.parents
    @boundscheck (1 <= i1) & (i1 <= length(P)) & checkbounds(Bool, P[1], itail...) ||
        Base.throw_boundserror(V, I)
    _unsafe_getindex(i1, itail, P...)
end

@inline function Base.setindex!(V::StackedView{T,N}, val, I::Vararg{Int,N}) where {T,N}
    i1, itail = I[1], tail(I)
    P = V.parents
    @boundscheck (1 <= i1) & (i1 <= length(P)) & checkbounds(Bool, P[1], itail...) ||
        Base.throw_boundserror(V, I)
    _unsafe_setindex!(i1, itail, convert(T, val), P...)
    val
end

# This does the same thing that V.parents[i1][itail...] would do,
# except in a type-stable way. The cost is the introduction of
# branches.
@inline function _unsafe_getindex(idx, I, A::AbstractArray, As...)
    if idx == 1
        @inbounds ret = A[I...]
        return ret
    end
    _unsafe_getindex(idx-1, I, As...)
end
_unsafe_getindex(idx, I) = error("ran out of arrays; this shouldn't happen")

# For getting all of the channels (e.g., for ColorView)
@inline function _unsafe_getindex_all(I, A, As...)
    @inbounds ret = A[I...]
    (ret, _unsafe_getindex_all(I, As...)...)
end
_unsafe_getindex_all(I) = ()

@inline function _unsafe_setindex!(idx, I, val, A::AbstractArray, As...)
    if idx == 1
        @inbounds A[I...] = val
        return val
    end
    _unsafe_setindex!(idx-1, I, val, As...)
end
_unsafe_setindex!(idx, I, val) = error("ran out of arrays; this shouldn't happen")

# For setting all of the channels (e.g., for ColorView), one `val` per parent
@inline function _unsafe_setindex_all!(I, vals::NTuple{N}, As::NTuple{N,AbstractArray}) where N
    val1, valrest = vals[1], tail(vals)
    A1, Arest = As[1], tail(As)
    @inbounds A1[I...] = val1
    _unsafe_setindex_all!(I, valrest, Arest)
end
_unsafe_setindex_all!(I, ::Tuple{}, ::Tuple{}) = nothing

# Performance optimizations
@inline getchannels(P::StackedView, ::Type{C}, I) where {C<:Color2} = _unsafe_getindex_all(I, P.parents...)
@inline getchannels(P::StackedView, ::Type{C}, I) where {C<:Color3} = _unsafe_getindex_all(I, P.parents...)
@inline getchannels(P::StackedView, ::Type{C}, I) where {C<:Color4} = _unsafe_getindex_all(I, P.parents...)

@inline setchannels!(P::StackedView, val::Color2, I) = _unsafe_setindex_all!(I, (comp1(val),alpha(val)), P.parents)
@inline setchannels!(P::StackedView, val::Color3, I) = _unsafe_setindex_all!(I, (comp1(val),comp2(val),comp3(val)), P.parents)
@inline setchannels!(P::StackedView, val::Color4, I) = _unsafe_setindex_all!(I, (comp1(val),comp2(val),comp3(val),alpha(val)), P.parents)


# When overlaying 2 images, you often might want one color channel to be black
struct ZeroArrayPromise{T} end
const zeroarray = ZeroArrayPromise{Union{}}()

struct ZeroArray{T,N,R<:AbstractUnitRange} <: AbstractArray{T,N}
    inds::NTuple{N,R}
end

ZeroArrayPromise{T}(inds::Tuple{R,Vararg{R,N}}) where {T,N,R<:AbstractUnitRange} = ZeroArray{T,N+1,R}(inds)
ZeroArrayPromise{T}(inds::NTuple{N,AbstractUnitRange}) where {T,N} = ZeroArrayPromise{T}(promote(inds...))
Base.eltype(::Type{ZeroArrayPromise{T}}) where {T} = T

Base.axes(A::ZeroArray) = A.inds
Base.size(A::ZeroArray) = length.(A.inds)
Base.getindex(A::ZeroArray{T,N}, I::Vararg{Int,N}) where {T,N} = zero(T)


@inline function StackedView(arrays::Union{AbstractArray,ZeroArrayPromise}...)
    T = promote_eleltype_all(arrays...)
    stackedview(T, arrays...)
end
@inline function StackedView{T}(arrays::Union{AbstractArray,ZeroArrayPromise}...) where T<:Number
    stackedview(T, arrays...)
end

# Now, it seems we should be able to do this with uppercase-typed
# calls, but there seems to be some inference bug, and using a
# function with a different name works around it.
@inline function stackedview(::Type{T}, arrays::Union{AbstractArray,ZeroArrayPromise}...) where T<:Number
    inds = firstinds(arrays...)
    arrays_take = take_zeros(T, inds, arrays...)
    arrays_T = map(A->of_eltype(zero(T), A), arrays_take)
    # To compute N = length(inds)+1 in a type-stable fashion, we have
    # to use tuple tricks (i.e., make a tuple of length(inds)+1)
    _stackedview(T, (length(arrays), inds...), arrays_T)
end
@inline _stackedview(::Type{T}, ::Tuple{Vararg{Any,N}}, arrays) where {T,N} = StackedView{T,N,typeof(arrays)}(arrays)


@inline firstinds(A::AbstractArray, Bs...) = axes(A)
@inline firstinds(::ZeroArrayPromise, Bs...) = firstinds(Bs...)
firstinds() = error("not all arrays can be zeroarray")

@inline take_zeros(::Type{T}, inds, ::ZeroArrayPromise, Bs...) where T = (ZeroArrayPromise{T}(inds), take_zeros(T, inds, Bs...)...)
@inline take_zeros(::Type{T}, inds, A::AbstractArray, Bs...) where T = (A, take_zeros(T, inds, Bs...)...)
take_zeros(T, inds) = ()

# Extensions of PaddedViews
function PaddedViews.paddedviews(fillvalue, As::Union{AbstractArray,ZeroArrayPromise}...)
    inds = PaddedViews.outerinds(As...)
    map(A->PaddedView(fillvalue, A, inds), As)
end

PaddedViews.PaddedView(fillvalue, zap::ZeroArrayPromise, indices) = zap

@inline PaddedViews.outerinds(A::ZeroArrayPromise, Bs...) = PaddedViews.outerinds(Bs...)
@inline PaddedViews._outerinds(inds, A::ZeroArrayPromise, Bs...) =
    PaddedViews._outerinds(inds, Bs...)
