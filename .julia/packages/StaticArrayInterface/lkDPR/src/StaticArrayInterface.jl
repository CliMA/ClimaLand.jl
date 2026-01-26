module StaticArrayInterface

@static if isdefined(Base, Symbol("@assume_effects"))
    using Base: @assume_effects
else
    macro assume_effects(args...)
        n = nfields(args)
        call = getfield(args, n)
        if n === 2 && getfield(args, 1) === QuoteNode(:total)
            return esc(:(Base.@pure $(call)))
        else
            return esc(call)
        end
    end
end

using PrecompileTools

@recompile_invalidations begin
    using ArrayInterface
    import ArrayInterface: allowed_getindex, allowed_setindex!, aos_to_soa, buffer,
                            parent_type, fast_matrix_colors, findstructralnz,
                            has_sparsestruct,
                            issingular, isstructured, matrix_colors, restructure,
                            lu_instance,
                            safevec, zeromatrix, undefmatrix, ColoringAlgorithm,
                            fast_scalar_indexing, parameterless_type,
                            is_forwarding_wrapper,
                                map_tuple_type, flatten_tuples, GetIndex, SetIndex!,
                            defines_strides, ndims_index, ndims_shape,
                            stride_preserving_index

    # ArrayIndex subtypes and methods
    import ArrayInterface: ArrayIndex, MatrixIndex, VectorIndex, BidiagonalIndex,
                            TridiagonalIndex
    # managing immutables
    import ArrayInterface: ismutable, can_change_size, can_setindex
    # constants
    import ArrayInterface: MatAdjTrans, VecAdjTrans, UpTri, LoTri
    # device pieces
    import ArrayInterface: AbstractDevice, AbstractCPU, CPUPointer, CPUTuple, CheckParent,
                            CPUIndex, GPU, can_avx, device

    using Static
    using Static: Zero, One, nstatic, eq, ne, gt, ge, lt, le, eachop, eachop_tuple,
                permute, invariant_permutation, field_type, reduce_tup, find_first_eq,
                OptionallyStaticUnitRange, OptionallyStaticStepRange, OptionallyStaticRange,
                IntType,
                SOneTo, SUnitRange

    using IfElse

    using Base.Cartesian
    using Base: @propagate_inbounds, tail, OneTo, LogicalIndex, Slice, ReinterpretArray,
                ReshapedArray, AbstractCartesianIndex

    using Base.Iterators: Pairs
    using LinearAlgebra

    import Compat
end

"""
    StrideIndex(x)
Subtype of `ArrayIndex` that transforms and index using stride layout information
derived from `x`.
"""
struct StrideIndex{N,R,C,S,O} <: ArrayIndex{N}
    strides::S
    offsets::O
    @inline function StrideIndex{N,R,C}(s::S, o::O) where {N,R,C,S,O}
        return new{N,R::NTuple{N,Int},C,S,O}(s, o)
    end
end

"""
    LazyAxis{N}(parent::AbstractArray)
A lazy representation of `axes(parent, N)`.
"""
struct LazyAxis{N,P} <: AbstractUnitRange{Int}
    parent::P

    function LazyAxis{N}(parent::P) where {N,P}
        N > 0 && return new{N::Int,P}(parent)
        throw_dim_error(parent, N)
    end
    @inline LazyAxis{:}(parent::P) where {P} = new{ifelse(ndims(P) === 1, 1, :),P}(parent)
end

function throw_dim_error(@nospecialize(x), @nospecialize(dim))
    throw(DimensionMismatch("$x does not have dimension corresponding to $dim"))
end

abstract type AbstractArray2{T, N} <: AbstractArray{T, N} end

"""
    BroadcastAxis
An abstract trait that is used to determine how axes are combined when calling `broadcast_axis`.
"""
abstract type BroadcastAxis end

@assume_effects :total function _find_first_true(isi::Tuple{Vararg{Union{Bool, Static.StaticBool}, N}}) where {N}
    for i in 1:N
        x = getfield(isi, i) 
        if (x isa Bool && x === true) || x isa Static.True
            return i
        end
    end
    return nothing
end

"""
    IndicesInfo{N}(inds::Tuple) -> IndicesInfo{N}(typeof(inds))
    IndicesInfo{N}(T::Type{<:Tuple}) -> IndicesInfo{N,pdims,cdims}()
    IndicesInfo(inds::Tuple) -> IndicesInfo(typeof(inds))
    IndicesInfo(T::Type{<:Tuple}) -> IndicesInfo{maximum(pdims),pdims,cdims}()


Maps a tuple of indices to `N` dimensions. The resulting `pdims` is a tuple where each
field in `inds` (or field type in `T`) corresponds to the parent dimensions accessed.
`cdims` similarly maps indices to the resulting child array produced after indexing with
`inds`. If `N` is not provided then it is assumed that all indices are represented by parent
dimensions and there are no trailing dimensions accessed. These may be accessed by through
`parentdims(info::IndicesInfo)` and `childdims(info::IndicesInfo)`. If `N` is not provided,
it is assumed that no indices are accessing trailing dimensions (which are represented as
`0` in `parentdims(info)[index_position]`).

The the fields and types of `IndicesInfo` should not be accessed directly.
Instead [`parentdims`](@ref), [`childdims`](@ref), [`ndims_index`](@ref), and
[`ndims_shape`](@ref) should be used to extract relevant information.

# Examples

```julia
julia> using StaticArrayInterface: IndicesInfo, parentdims, childdims, ndims_index, ndims_shape

julia> info = IndicesInfo{5}(typeof((:,[CartesianIndex(1,1),CartesianIndex(1,1)], 1, ones(Int, 2, 2), :, 1)));

julia> parentdims(info)  # the last two indices access trailing dimensions
(1, (2, 3), 4, 5, 0, 0)

julia> childdims(info)
(1, 2, 0, (3, 4), 5, 0)

julia> childdims(info)[3]  # index 3 accesses a parent dimension but is dropped in the child array
0

julia> ndims_index(info)
5

julia> ndims_shape(info)
5

julia> info = IndicesInfo(typeof((:,[CartesianIndex(1,1),CartesianIndex(1,1)], 1, ones(Int, 2, 2), :, 1)));

julia> parentdims(info)  # assumed no trailing dimensions
(1, (2, 3), 4, 5, 6, 7)

julia> ndims_index(info)  # assumed no trailing dimensions
7

```
"""
struct IndicesInfo{Np, pdims, cdims, Nc}
    function IndicesInfo{N}(@nospecialize(T::Type{<:Tuple})) where {N}
        SI = _find_first_true(map_tuple_type(is_splat_index, T))
        NI = map_tuple_type(ndims_index, T)
        NS = map_tuple_type(ndims_shape, T)
        if SI === nothing
            ndi = NI
            nds = NS
        else
            nsplat = N - sum(NI)
            if nsplat === 0
                ndi = NI
                nds = NS
            else
                splatmul = max(0, nsplat + 1)
                ndi = _map_splats(splatmul, SI, NI)
                nds = _map_splats(splatmul, SI, NS)
            end
        end
        if ndi === (1,) && N !== 1
            ns1 = getfield(nds, 1)
            new{N, (:,), (ns1 > 1 ? ntuple(identity, ns1) : ns1,), ns1}()
        else
            nds_cumsum = cumsum(nds)
            if sum(ndi) > N
                init_pdims = _accum_dims(cumsum(ndi), ndi)
                pdims = ntuple(nfields(init_pdims)) do i
                    dim_i = getfield(init_pdims, i)
                    if dim_i isa Tuple
                        ntuple(length(dim_i)) do j
                            dim_i_j = getfield(dim_i, j)
                            dim_i_j > N ? 0 : dim_i_j
                        end
                    else
                        dim_i > N ? 0 : dim_i
                    end
                end
                new{N, pdims, _accum_dims(nds_cumsum, nds), last(nds_cumsum)}()
            else
                new{N, _accum_dims(cumsum(ndi), ndi), _accum_dims(nds_cumsum, nds),
                    last(nds_cumsum)}()
            end
        end
    end
    IndicesInfo{N}(@nospecialize(t::Tuple)) where {N} = IndicesInfo{N}(typeof(t))
    function IndicesInfo(@nospecialize(T::Type{<:Tuple}))
        ndi = map_tuple_type(ndims_index, T)
        nds = map_tuple_type(ndims_shape, T)
        ndi_sum = cumsum(ndi)
        nds_sum = cumsum(nds)
        nf = nfields(ndi_sum)
        pdims = _accum_dims(ndi_sum, ndi)
        cdims = _accum_dims(nds_sum, nds)
        new{getfield(ndi_sum, nf), pdims, cdims, getfield(nds_sum, nf)}()
    end
    IndicesInfo(@nospecialize t::Tuple) = IndicesInfo(typeof(t))
    @inline function IndicesInfo(@nospecialize T::Type{<:SubArray})
        IndicesInfo{ndims(parent_type(T))}(fieldtype(T, :indices))
    end
    IndicesInfo(x::SubArray) = IndicesInfo{ndims(parent(x))}(typeof(x.indices))
end

@inline function _map_splats(nsplat::Int, splat_index::Int, dims::Tuple{Vararg{Union{Int,StaticInt}}})
    ntuple(length(dims)) do i
        i === splat_index ? (nsplat * getfield(dims, i)) : getfield(dims, i)
    end
end

@inline function _accum_dims(csdims::NTuple{N, Int}, nd::NTuple{N, Int}) where {N}
    ntuple(N) do i
        nd_i = getfield(nd, i)
        if nd_i === 0
            0
        elseif nd_i === 1
            getfield(csdims, i)
        else
            ntuple(Base.Fix1(+, getfield(csdims, i) - nd_i), nd_i)
        end
    end
end

function _lower_info(::IndicesInfo{Np, pdims, cdims, Nc}) where {Np, pdims, cdims, Nc}
    Np, pdims, cdims, Nc
end

ndims_index(@nospecialize(info::IndicesInfo)) = getfield(_lower_info(info), 1)
ndims_shape(@nospecialize(info::IndicesInfo)) = getfield(_lower_info(info), 4)

"""
    parentdims(::IndicesInfo) -> Tuple

Returns the parent dimension mapping from `IndicesInfo`.

See also: [`IndicesInfo`](@ref), [`childdims`](@ref)
"""
parentdims(@nospecialize info::IndicesInfo) = getfield(_lower_info(info), 2)

"""
    childdims(::IndicesInfo) -> Tuple

Returns the child dimension mapping from `IndicesInfo`.

See also: [`IndicesInfo`](@ref), [`parentdims`](@ref)
"""
childdims(@nospecialize info::IndicesInfo) = getfield(_lower_info(info), 3)

@generated function merge_tuple_type(::Type{X}, ::Type{Y}) where {X <: Tuple, Y <: Tuple}
    Tuple{X.parameters..., Y.parameters...}
end

Base.size(A::AbstractArray2) = map(Int, static_size(A))
Base.size(A::AbstractArray2, dim) = Int(static_size(A, dim))

function Base.axes(A::AbstractArray2)
    is_forwarding_wrapper(A) && return static_axes(parent(A))
    throw(ArgumentError("Subtypes of `AbstractArray2` must define an axes method"))
end
function Base.axes(A::AbstractArray2, dim::Union{Symbol, StaticSymbol})
    static_axes(A, to_dims(A, dim))
end

function Base.strides(A::AbstractArray2)
    defines_strides(A) && return map(Int, static_strides(A))
    throw(MethodError(Base.strides, (A,)))
end
Base.strides(A::AbstractArray2, dim) = Int(static_strides(A, dim))

function Base.IndexStyle(::Type{T}) where {T <: AbstractArray2}
    is_forwarding_wrapper(T) ? IndexStyle(parent_type(T)) : IndexCartesian()
end

function Base.length(A::AbstractArray2)
    len = known_length(A)
    if len === nothing
        return Int(prod(static_size(A)))
    else
        return Int(len)
    end
end

@propagate_inbounds Base.getindex(A::AbstractArray2, args...) = static_getindex(A, args...)
@propagate_inbounds Base.getindex(A::AbstractArray2; kwargs...) = static_getindex(A; kwargs...)

@propagate_inbounds function Base.setindex!(A::AbstractArray2, val, args...)
    return setindex!(A, val, args...)
end
@propagate_inbounds function Base.setindex!(A::AbstractArray2, val; kwargs...)
    return setindex!(A, val; kwargs...)
end

@inline static_first(x) = Static.maybe_static(known_first, first, x)
@inline static_last(x) = Static.maybe_static(known_last, last, x)
@inline static_step(x) = Static.maybe_static(known_step, step, x)

@inline function _to_cartesian(a, i::IntType)
    @inbounds(CartesianIndices(ntuple(dim -> indices(a, dim), Val(ndims(a))))[i])
end
@inline function _to_linear(a, i::Tuple{IntType, Vararg{IntType}})
    _strides2int(offsets(a), size_to_strides(static_size(a), static(1)), i) + static(1)
end

"""
    has_parent(::Type{T}) -> StaticBool

Returns `static(true)` if `parent_type(T)` a type unique to `T`.
"""
has_parent(x) = has_parent(typeof(x))
has_parent(::Type{T}) where {T} = _has_parent(parent_type(T), T)
_has_parent(::Type{T}, ::Type{T}) where {T} = False()
_has_parent(::Type{T1}, ::Type{T2}) where {T1, T2} = True()

"""
    is_lazy_conjugate(::AbstractArray) -> Bool

Determine if a given array will lazyily take complex conjugates, such as with `Adjoint`. This will work with
nested wrappers, so long as there is no type in the chain of wrappers such that `parent_type(T) == T`

Examples

    julia> a = transpose([1 + im, 1-im]')
    2×1 transpose(adjoint(::Vector{Complex{Int64}})) with eltype Complex{Int64}:
     1 - 1im
     1 + 1im

    julia> is_lazy_conjugate(a)
    True()

    julia> b = a'
    1×2 adjoint(transpose(adjoint(::Vector{Complex{Int64}}))) with eltype Complex{Int64}:
     1+1im  1-1im

    julia> is_lazy_conjugate(b)
    False()

"""
is_lazy_conjugate(::T) where {T <: AbstractArray} = _is_lazy_conjugate(T, False())
is_lazy_conjugate(::AbstractArray{<:Real}) = False()

function _is_lazy_conjugate(::Type{T}, isconj) where {T <: AbstractArray}
    Tp = parent_type(T)
    if T !== Tp
        _is_lazy_conjugate(Tp, isconj)
    else
        isconj
    end
end

function _is_lazy_conjugate(::Type{T}, isconj) where {T <: Adjoint}
    Tp = parent_type(T)
    if T !== Tp
        _is_lazy_conjugate(Tp, !isconj)
    else
        !isconj
    end
end

"""
    insert(collection, index, item)

Returns a new instance of `collection` with `item` inserted into at the given `index`.
"""
Base.@propagate_inbounds function insert(collection, index, item)
    @boundscheck checkbounds(collection, index)
    ret = similar(collection, static_length(collection) + 1)
    @inbounds for i in firstindex(ret):(index - 1)
        ret[i] = collection[i]
    end
    @inbounds ret[index] = item
    @inbounds for i in (index + 1):lastindex(ret)
        ret[i] = collection[i - 1]
    end
    return ret
end

function insert(x::Tuple{Vararg{Any, N}}, index, item) where {N}
    @boundscheck if !checkindex(Bool, StaticInt{1}():StaticInt{N}(), index)
        throw(BoundsError(x, index))
    end
    return unsafe_insert(x, Int(index), item)
end

@inline function unsafe_insert(x::Tuple, i::Int, item)
    if i === 1
        return (item, x...)
    else
        return (first(x), unsafe_insert(tail(x), i - 1, item)...)
    end
end

"""
    deleteat(collection, index)

Returns a new instance of `collection` with the item at the given `index` removed.
"""
Base.@propagate_inbounds function deleteat(collection::AbstractVector, index)
    @boundscheck if !checkindex(Bool, eachindex(collection), index)
        throw(BoundsError(collection, index))
    end
    return unsafe_deleteat(collection, index)
end
Base.@propagate_inbounds function deleteat(collection::Tuple{Vararg{Any, N}},
                                           index) where {N}
    @boundscheck if !checkindex(Bool, StaticInt{1}():StaticInt{N}(), index)
        throw(BoundsError(collection, index))
    end
    return unsafe_deleteat(collection, index)
end

function unsafe_deleteat(src::AbstractVector, index)
    dst = similar(src, static_length(src) - 1)
    @inbounds for i in indices(dst)
        if i < index
            dst[i] = src[i]
        else
            dst[i] = src[i + 1]
        end
    end
    return dst
end

@inline function unsafe_deleteat(src::AbstractVector, inds::AbstractVector)
    dst = similar(src, static_length(src) - static_length(inds))
    dst_index = firstindex(dst)
    @inbounds for src_index in indices(src)
        if !in(src_index, inds)
            dst[dst_index] = src[src_index]
            dst_index += one(dst_index)
        end
    end
    return dst
end

@inline function unsafe_deleteat(src::Tuple, inds::AbstractVector)
    dst = Vector{eltype(src)}(undef, static_length(src) - static_length(inds))
    dst_index = firstindex(dst)
    @inbounds for src_index in static(1):static_length(src)
        if !in(src_index, inds)
            dst[dst_index] = src[src_index]
            dst_index += one(dst_index)
        end
    end
    return Tuple(dst)
end

@inline unsafe_deleteat(x::Tuple{T}, i) where {T} = ()
@inline unsafe_deleteat(x::Tuple{T1, T2}, i) where {T1, T2} = isone(i) ? (x[2],) : (x[1],)
@inline function unsafe_deleteat(x::Tuple, i)
    if i === one(i)
        return tail(x)
    elseif i == static_length(x)
        return Base.front(x)
    else
        return (first(x), unsafe_deleteat(tail(x), i - one(i))...)
    end
end

include("array_index.jl")
include("ranges.jl")
include("axes.jl")
include("size.jl")
include("dimensions.jl")
include("indexing.jl")
include("stridelayout.jl")
include("broadcast.jl")

## Precompilation

@setup_workload begin
    # Putting some things in `setup` can reduce the size of the
    # precompile file and potentially make loading faster.
    arrays = [rand(4), Base.oneto(5)]
    @compile_workload begin for x in arrays
        known_first(x)
        known_step(x)
        known_last(x)
    end end
end

end
