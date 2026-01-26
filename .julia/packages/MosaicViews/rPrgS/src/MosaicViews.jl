module MosaicViews

using PaddedViews
using OffsetArrays
using MappedArrays: of_eltype
using StackViews

export

    MosaicView,
    mosaicview,
    mosaic

"""
    MosaicView(A::AbstractArray)

Create a two dimensional "view" of the three or four dimensional
array `A`. The resulting `MosaicView` will display the data in
`A` such that it emulates using `vcat` for all elements in the
third dimension of `A`, and `hcat` for all elements in the fourth
dimension of `A`.

For example, if `size(A)` is `(2,3,4)`, then the resulting
`MosaicView` will have the size `(2*4,3)` which is `(8,3)`.
Alternatively, if `size(A)` is `(2,3,4,5)`, then the resulting
size will be `(2*4,3*5)` which is `(8,15)`.

Another way to think about this is that `MosaicView` creates a
mosaic of all the individual matrices enumerated in the third
(and optionally fourth) dimension of the given 3D or 4D array
`A`. This can be especially useful for creating a single
composite image from a set of equally sized images.

```@jldoctest
julia> using MosaicViews

julia> A = [(k+1)*l-1 for i in 1:2, j in 1:3, k in 1:2, l in 1:2]
2×3×2×2 Array{Int64,4}:
[:, :, 1, 1] =
 1  1  1
 1  1  1

[:, :, 2, 1] =
 2  2  2
 2  2  2

[:, :, 1, 2] =
 3  3  3
 3  3  3

[:, :, 2, 2] =
 5  5  5
 5  5  5

julia> MosaicView(A)
4×6 MosaicViews.MosaicView{Int64,4,Array{Int64,4}}:
 1  1  1  3  3  3
 1  1  1  3  3  3
 2  2  2  5  5  5
 2  2  2  5  5  5
```
"""
struct MosaicView{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,2}
    parent::A
    dims::Tuple{Int,Int}
    pdims::NTuple{N,Int}

    function MosaicView{T,N}(A::AbstractArray{T,N}, dims) where {T,N}
        1 <= N <= 4 || throw(ArgumentError("The given array must have dimensionality N=3 or N=4"))
        # unless we store the axes in the struct, we can't support offset indices
        # N=2 is a special case that we can provide a specialization on `axes`
        # but for consistency with cases of other dimensions, disable it as well
        require_one_based_indexing(A)
        new{T,N,typeof(A)}(A, dims, size(A))
    end
end

function MosaicView(A::AbstractArray{T,N}) where {T,N}
    # vectors are lifted to 2D matrix, 3d/4d arrays are reshaped to 2D matrix
    dims = (size(A,1) * size(A,3), size(A,2) * size(A,4))
    MosaicView{T,N}(A, dims)
end

Base.parent(mva::MosaicView) = mva.parent
Base.size(mva::MosaicView) = mva.dims

# fallback for 1d/2d case
@inline Base.getindex(mva::MosaicView, ind::Int...) = mva.parent[ind...]

@inline function Base.getindex(mva::MosaicView{T,3,A}, i::Int, j::Int) where {T,A}
    @boundscheck checkbounds(mva, i, j)
    pdims = mva.pdims
    parent = mva.parent
    idx1 = (i-1) % pdims[1] + 1
    idx2 = (j-1) % pdims[2] + 1
    idx3 = (i-1) ÷ pdims[1] + 1
    @inbounds res = parent[idx1, idx2, idx3]
    res
end

# FIXME: we need the annotation T because mosaicview + StackView is currently not type stable
@inline function Base.getindex(mva::MosaicView{T,4,A}, i::Int, j::Int)::T where {T,A}
    @boundscheck checkbounds(mva, i, j)
    pdims = mva.pdims
    parent = mva.parent
    idx1 = (i-1) % pdims[1] + 1
    idx2 = (j-1) % pdims[2] + 1
    idx3 = (i-1) ÷ pdims[1] + 1
    idx4 = (j-1) ÷ pdims[2] + 1
    @inbounds res = parent[idx1, idx2, idx3, idx4]
    res
end

"""
    mosaicview(A::AbstractArray;
               [fillvalue=<zero unit>], [npad=0],
               [nrow], [ncol], [rowmajor=false]) -> MosaicView
    mosaicview(As::AbstractArray...; kwargs...)
    mosaicview(As::Union{Tuple, AbstractVector}; kwargs...)

Create a two dimensional "view" from array `A`.

The resulting [`MosaicView`](@ref) will display all the matrix
slices of the first two dimensions of `A` arranged as a single
large mosaic (in the form of a matrix).

# Arguments

In contrast to using the constructor of [`MosaicView`](@ref)
directly, the function `mosaicview` also allows for a couple of
convenience keywords. A typical use case would be to create an
image mosaic from a set of input images.

- The parameter `fillvalue` defines the value that
  that should be used for empty space. This can be padding caused
  by `npad`, or empty mosaic tiles in case the number of matrix
  slices in `A` is smaller than `nrow*ncol`.

- The parameter `npad` defines the empty padding space between
  adjacent mosaic tiles. This can be especially useful if the
  individual tiles (i.e. matrix slices in `A`) are images that
  should be visually separated by some grid lines.

- The parameters `nrow` and `ncol` can be used to choose the
  number of tiles in row and/or column direction the mosaic should
  be arranged in. Note that it suffices to specify one of the
  two parameters, as the other one can be inferred accordingly.
  The default in case none of the two are specified is `nrow = size(A,3)`.

- If `rowmajor` is set to `true`, then the slices will be
  arranged left-to-right-top-to-bottom, instead of
  top-to-bottom-left-to-right (default). The layout only differs
  in non-trivial cases, i.e., when `nrow != 1` and `ncol != 1`.

!!! tip
    This function is not type stable and should only be used if
    performance is not a priority. To achieve optimized performance,
    you need to manually construct a [`MosaicView`](@ref).

# Examples

The simplest usage is to `cat` two arrays of the same dimension.

```julia-repl
julia> A1 = fill(1, 3, 1)
3×1 Array{Int64,2}:
 1
 1
 1

julia> A2 = fill(2, 1, 3)
1×3 Array{Int64,2}:
 2  2  2

julia> mosaicview(A1, A2)
6×3 MosaicView{Int64,4, ...}:
 0  1  0
 0  1  0
 0  1  0
 0  0  0
 2  2  2
 0  0  0

julia> mosaicview(A1, A2; center=false)
 6×3 MosaicView{Int64,4, ...}:
  1  0  0
  1  0  0
  1  0  0
  2  2  2
  0  0  0
  0  0  0
```

Other keyword arguments can be useful to get a nice looking results.

```julia-repl
julia> using MosaicViews

julia> A = [k for i in 1:2, j in 1:3, k in 1:5]
2×3×5 Array{Int64,3}:
[:, :, 1] =
 1  1  1
 1  1  1

[:, :, 2] =
 2  2  2
 2  2  2

[:, :, 3] =
 3  3  3
 3  3  3

[:, :, 4] =
 4  4  4
 4  4  4

[:, :, 5] =
 5  5  5
 5  5  5

julia> mosaicview(A, ncol=2)
6×6 MosaicViews.MosaicView{Int64,4,...}:
 1  1  1  4  4  4
 1  1  1  4  4  4
 2  2  2  5  5  5
 2  2  2  5  5  5
 3  3  3  0  0  0
 3  3  3  0  0  0

julia> mosaicview(A, nrow=2)
4×9 MosaicViews.MosaicView{Int64,4,...}:
 1  1  1  3  3  3  5  5  5
 1  1  1  3  3  3  5  5  5
 2  2  2  4  4  4  0  0  0
 2  2  2  4  4  4  0  0  0

julia> mosaicview(A, nrow=2, rowmajor=true)
4×9 MosaicViews.MosaicView{Int64,4,...}:
 1  1  1  2  2  2  3  3  3
 1  1  1  2  2  2  3  3  3
 4  4  4  5  5  5  0  0  0
 4  4  4  5  5  5  0  0  0

julia> mosaicview(A, nrow=2, npad=1, rowmajor=true)
5×11 MosaicViews.MosaicView{Int64,4,...}:
 1  1  1  0  2  2  2  0  3  3  3
 1  1  1  0  2  2  2  0  3  3  3
 0  0  0  0  0  0  0  0  0  0  0
 4  4  4  0  5  5  5  0  0  0  0
 4  4  4  0  5  5  5  0  0  0  0

julia> mosaicview(A, fillvalue=-1, nrow=2, npad=1, rowmajor=true)
5×11 MosaicViews.MosaicView{Int64,4,...}:
  1   1   1  -1   2   2   2  -1   3   3   3
  1   1   1  -1   2   2   2  -1   3   3   3
 -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1
  4   4   4  -1   5   5   5  -1  -1  -1  -1
  4   4   4  -1   5   5   5  -1  -1  -1  -1
```
"""
function mosaicview(A::AbstractArray{T,3};
                    fillvalue = zero(T),
                    npad = 0,
                    nrow = -1,
                    ncol = -1,
                    rowmajor = false,
                    kwargs...) where T   # delete `kwargs...` when we delete the `center` depwarn
    nrow == -1 || nrow > 0 || throw(ArgumentError("The parameter \"nrow\" must be greater than 0"))
    ncol == -1 || ncol > 0 || throw(ArgumentError("The parameter \"ncol\" must be greater than 0"))
    npad >= 0 || throw(ArgumentError("The parameter \"npad\" must be greater than or equal to 0"))
    if !isempty(kwargs)
        if haskey(kwargs, :center)
            Base.depwarn("use of `center` in `mosaicview` is deprecated; see `mosaic`", :mosaicview)
        else
            Base.depwarn("passing extraneous keyword arguments $(values(kwargs)) to `mosaicview` is deprecated", :mosaicview)
        end
    end
    ntile = size(A,3)
    ntile_ceil = ntile # ntile need not be integer divideable
    if nrow == -1 && ncol == -1
        # automatically choose nrow to reflect what MosaicView does
        nrow = ntile
        ncol = 1
    elseif nrow == -1
        # compute nrow based on ncol
        nrow = ceil(Int, ntile / ncol)
        ntile_ceil = nrow * ncol
    elseif ncol == -1
        # compute ncol based on nrow
        ncol = ceil(Int, ntile / nrow)
        ntile_ceil = nrow * ncol
    else
        # accept nrow and ncol as is if it covers at least all
        # existing tiles
        ntile_ceil = nrow * ncol
        ntile_ceil < ntile && throw(ArgumentError("The product of the parameters \"ncol\" (value: $ncol) and \"nrow\" (value: $nrow) must be greater than or equal to $ntile"))
    end
    # we pad size(A,3) to nrow*ncol. we also pad the first two
    # dimensions according to npad. think of this as border
    # between tiles (useful for images)
    pad_dims = (size(A,1) + npad, size(A,2) + npad, ntile_ceil)
    A_pad = PaddedView(convert(T, fillvalue), A, pad_dims)
    # next we reshape the image such that it reflects the
    # specified nrow and ncol
    A_new = if !rowmajor
        res_dims = (size(A_pad,1), size(A_pad,2), nrow, ncol)
        reshape(A_pad, res_dims)
    else
        # same as above but we additionally permute dimensions
        # to mimic row first layout for the tiles. this is useful
        # for images since user often reads these from left-to-right
        # before top-to-bottom. (note the swap of "ncol" and "nrow")
        res_dims = (size(A_pad,1), size(A_pad,2), ncol, nrow)
        A_tp = reshape(A_pad, res_dims)
        PermutedDimsArray{eltype(A_tp), 4, (1, 2, 4, 3), (1, 2, 4, 3), typeof(A_tp)}(A_tp)
    end
    # decrease size of the resulting MosaicView by npad to not have
    # a border on the right side and bottom side of the final mosaic.
    dims = (size(A_new,1) * size(A_new,3) - npad, size(A_new,2) * size(A_new,4) - npad)
    MosaicView{T,4}(A_new, dims)
end

function mosaicview(A::AbstractArray{T,N};
                    nrow = -1,
                    ncol = -1,
                    kwargs...) where {T,N}
    # if neither nrow nor ncol is provided then automatically choose
    # nrow and ncol to reflect what MosaicView does (i.e. use size)
    if nrow == -1 && ncol == -1
        nrow = size(A, 3)
        # ncol = size(A, 4)
    end
    mosaicview(reshape(A, (size(A,1), size(A,2), :));
               nrow=nrow, ncol=ncol, kwargs...)
end

"""
    mosaic(A1, A2...; center=true, kwargs...)
    mosaic([A1, A2, ...]; center=true, kwargs...)

Create a mosaic out of input arrays `A1`, `A2`, .... `mosaic` is essentially
a more flexible version of `cat` or `hvcat`; like them it makes a copy of
the inputs rather than returning a "view."

If `center` is set to `true`, then the padded arrays will be shifted
to the center; if set to false, they shift to the top-left corner. This
parameter is only useful when arrays are of different sizes.

All the keyword arguments of [`mosaicview`](@ref) are also supported.
"""
@inline mosaic(As::AbstractArray...; kwargs...) = mosaic(As; kwargs...)

function mosaic(As::AbstractVector{<:AbstractArray};
                fillvalue=zero(_filltype(As)),
                center::Bool=true,
                kwargs...)
    length(As) == 0 && throw(ArgumentError("The given vector should not be empty"))
    nd = ndims(As[1])
    all(A->ndims(A)==nd, As) || throw(ArgumentError("All arrays should have the same dimension"))
    T = _filltype(As)
    fillvalue = convert(T, fillvalue)
    mosaicview(_padded_cat(As; center=center, fillvalue=fillvalue, dims=valdim(first(As)));
               fillvalue=fillvalue, kwargs...)
end

function mosaic(As::Tuple;
                fillvalue=zero(_filltype(As)),
                center::Bool=true,
                kwargs...)
    length(As) == 0 && throw(ArgumentError("The given tuple should not be empty"))
    nd = ndims(As[1])
    all(A->ndims(A)==nd, As) || throw(ArgumentError("All arrays should have the same dimension"))
    T = _filltype(As)
    fillvalue = convert(T, fillvalue)
    vd = valdim(first(As))
    if isconcretetype(eltype(As)) || VERSION < v"1.2.0"
        # Base.inferencebarrier requires Julia at least v1.2.0
        mosaicview(_padded_cat(As; center=center, fillvalue=fillvalue, dims=vd);
                fillvalue=fillvalue, kwargs...)
    else
        # Reduce latency by despecializing calls with heterogeneous array types
        mosaicview(_padded_cat(Base.inferencebarrier(As); center=center, fillvalue=Base.inferencebarrier(fillvalue), dims=Base.inferencebarrier(vd));
                fillvalue=fillvalue, kwargs...)
    end
end

valdim(A::AbstractArray{T,0}) where T     = Val(3)
valdim(A::AbstractVector)                 = Val(3)
valdim(A::AbstractMatrix)                 = Val(3)
valdim(A::AbstractArray{T,N}) where {T,N} = Val(N+1)

function _padded_cat(imgs, center::Bool, fillvalue, dims)
    @nospecialize # because this is frequently called with heterogeneous inputs, we @nospecialize it
    pv(@nospecialize(imgs::AbstractVector{<:AbstractArray})) = PaddedViews.paddedviews_itr(fillvalue, imgs)
    pv(@nospecialize(imgs)) = paddedviews(fillvalue, imgs...)
    sym_pv(@nospecialize(imgs::AbstractVector{<:AbstractArray})) = PaddedViews.sym_paddedviews_itr(fillvalue, imgs)
    sym_pv(@nospecialize(imgs)) = sym_paddedviews(fillvalue, imgs...)

    pv_fn = center ? sym_pv : pv
    return StackView{_filltype(imgs)}(pv_fn(imgs), dims)
end
# compat: some packages uses this method
_padded_cat(imgs; center::Bool, fillvalue, dims) = _padded_cat(imgs, center, fillvalue, dims)

has_common_axes(@nospecialize(imgs)) = isempty(imgs) || all(isequal(axes(first(imgs))) ∘ axes, imgs)


# This uses Union{} as a sentinel eltype (all other types "beat" it),
# and Bool as a near-neutral fill type.
_filltype(As) = PaddedViews.filltype(Bool, _filltypeT(Union{}, As...))

@inline _filltypeT(::Type{T}, A, tail...) where T = _filltypeT(promote_wrapped_type(T, _gettype(A)), tail...)
_filltypeT(::Type{T}) where T = T

# When the inputs are homogenous we can circumvent varargs despecialization
# This also handles the case of empty `As` but concrete `T`.
function _filltype(As::AbstractVector{A}) where A<:AbstractArray{T} where T
    # (!@isdefined(T) || T === Any) && return invoke(_filltype, Tuple{Any}, As)
    T === Any && return invoke(_filltype, Tuple{Any}, As)
    return PaddedViews.filltype(Bool, T)
end

_gettype(A::AbstractArray{T}) where T = T === Any ? typeof(first(A)) : T

"""
    promote_wrapped_type(S, T)

Similar to `promote_type`, except designed to be extensible to cases where you want to promote should occur through a wrapper type.

`promote_wrapped_type` is used by `_filltype` to compute the common element type for handling heterogeneous types when building the mosaic.
It does not have the order-independence of `promote_type`, and you should extend it directly rather than via a `promote_rule`-like mechanism.

# Example

Suppose you have
```
struct MyWrapper{T}
    x::T
end
```
and you don't want to define `promote_type(MyWrapper{Int},Float32)` generally as anything other than `Any`,
but for the purpose of building mosaics a `MyWrapper{Float32}` would be a valid common type.
Then you could define

```
MosaicViews.promote_wrapped_type(::Type{MyWrapper{S}}, ::Type{MyWrapper{T}}) where {S,T} = MyWrapper{MosaicViews.promote_wrapped_type(S,T)}
MosaicViews.promote_wrapped_type(::Type{MyWrapper{S}}, ::Type{T}) where {S,T} = MyWrapper{MosaicViews.promote_wrapped_type(S,T)}
MosaicViews.promote_wrapped_type(::Type{S}, ::Type{MyWrapper{T}}) where {S,T} = MosaicViews.promote_wrapped_type(MyWrapper{T}, S)
```
"""
promote_wrapped_type(::Type{S}, ::Type{T}) where {S, T} = promote_type(S, T)

### compat
if VERSION < v"1.2"
    require_one_based_indexing(A...) = !Base.has_offset_axes(A...) || throw(ArgumentError("offset arrays are not supported but got an array with index other than 1"))
else
    const require_one_based_indexing = Base.require_one_based_indexing
end

### deprecations

@deprecate mosaicview(A1::AbstractArray, A2::AbstractArray; kwargs...) mosaic(A1, A2; kwargs...) # prevent A2 from being interpreted as fillvalue
@deprecate mosaicview(As::AbstractArray...; kwargs...) mosaic(As...; kwargs...)
@deprecate mosaicview(As::AbstractVector{<:AbstractArray};
                      fillvalue=zero(_filltype(As)),
                      center::Bool=true,
                      kwargs...) mosaic(As; fillvalue=fillvalue, center=center, kwargs...)
@deprecate mosaicview(As::Tuple;
                      fillvalue=zero(_filltype(As)),
                      center::Bool=true,
                      kwargs...) mosaic(As; fillvalue=fillvalue, center=center, kwargs...)

end # module
