"""
    restrict(img[, dims]) -> imgr

Reduce the size of `img` by approximately two-fold along the dimensions listed in
`dims`, or all spatial coordinates if `dims` is not specified.

# Output

The type of output array `imgr` depends on the input type:

- If `img` is an `OffsetArray`, then output array `imgr` will also be an `OffsetArray`.
- If `img` is not an `OffsetArray`, then output array `imgr` will be an `Array` type
  even if it has offset indices.

The size of `imgr` is approximately `1/2` of the original size. More specifically:

- if `Nₖ = size(img, k)` is odd, then `size(imgr, k) == (Nₖ+1) ÷ 2`.
- if `Nₖ = size(img, k)` is even, then `size(imgr, k) == (Nₖ÷2) + 1`.

# Examples

The optional argument `dims` can be a `Tuple` or `Integer`:

```julia
A = rand(5, 5) # size: (5, 5)

restrict(A) # size: (3, 3)

restrict(A, 1) # size: (3, 5)
restrict(A, 2) # size: (5, 3)

restrict(A, (1, )) # size: (3, 5)
restrict(A, (1, 2)) # size: (3, 3)
```

Unless the input array is 1-based, the origin will be halfed:

```julia
julia> using ImageBase, OffsetArrays

julia> Ao = OffsetArray(rand(5, 4), 5, 6);

julia> Ar = restrict(Ao);

julia> axes(Ao)
(OffsetArrays.IdOffsetRange(values=6:10, indices=6:10), OffsetArrays.IdOffsetRange(values=7:10, indices=7:10))

julia> axes(Ar)
(OffsetArrays.IdOffsetRange(values=3:5, indices=3:5), OffsetArrays.IdOffsetRange(values=4:6, indices=4:6))
```

# Extended help

The term `restrict` is taken from the coarsening operation of algebraic multigrid
methods; it is the adjoint of "prolongation" (which is essentially interpolation).
`restrict` anti-aliases the image as it goes, so is better than a naive summation
over 2x2 blocks.
The implementation of `restrict` has been tuned for performance, and should be
a fast method for constructing [pyramids](https://en.wikipedia.org/wiki/Pyramid_(image_processing)).

If `l` is the size of `img` along a particular dimension, `restrict` produces an
array of size `(l+1)÷2` for odd `l`,
and `l÷2 + 1` for even `l`. See the example below for an explanation.

See also `ImageTransformations.imresize`.

# Example

```julia
a_coarse = [0, 1, 0.3]
```
If we were to interpolate this at the halfway points, we'd get
```julia
a_fine = [0, 0.5, 1, 0.65, 0.3]
```
Note that `a_fine` is obtained from `a_coarse` via the *prolongation* operator `P` as
`P*a_coarse`, where
```julia
P = [1   0   0;      # this line "copies over" the first point
     0.5 0.5 0;      # this line takes the mean of the first and second point
     0   1   0;      # copy the second point
     0   0.5 0.5;    # take the mean of the second and third
     0   0   1]      # copy the third
```
`restrict` is the adjoint of prolongation. Consequently,
```julia
julia> restrict(a_fine)
3-element Array{Float64,1}:
 0.125
 0.7875
 0.3125

julia> (P'*a_fine)/2
3-element Array{Float64,1}:
 0.125
 0.7875
 0.3125
```
where the division by 2 approximately preserves the mean intensity of the input.

As we see here, for odd-length `a_fine`, restriction is the adjoint of interpolation at half-grid points.
When `length(a_fine)` is even, restriction is the adjoint of interpolation at 1/4 and 3/4-grid points.
This turns out to be the origin of the `l->l÷2 + 1` behavior.

One consequence of this definition is that the edges move towards zero:
```julia
julia> restrict(ones(11))
6-element Array{Float64,1}:
 0.75
 1.0
 1.0
 1.0
 1.0
 0.75
```
In some applications (e.g., image registration), you may find it useful to trim the edges.
"""
restrict(img::AbstractArray, ::Tuple{}) = img

restrict(A::AbstractArray) = restrict(A, coords_spatial(A))
function restrict(A::AbstractArray, dims::Dims)
    restrict(restrict(A, dims[1]), Base.tail(dims))
end

function restrict(A::AbstractArray{T,N}, dim::Integer) where {T,N}
    if Base.has_offset_axes(A)
        # For type stability, we cannot return `OffsetArray` in this method
        A = OffsetArrays.no_offset_view(A)
    end

    indsA = axes(A)
    newinds = ntuple(i->i==dim ? restrict_indices(indsA[dim]) : indsA[i], Val(N))
    out = Array{restrict_eltype(first(A)), N}(undef, last.(newinds))
    restrict!(out, A, dim)
    out
end
function restrict(A::OffsetArray{T,N}, dim::Integer) where {T,N}
    indsA = axes(A)
    newinds = map(UnitRange, ntuple(i->i==dim ? restrict_indices(indsA[dim]) : indsA[i], Val(N)))
    # By shifting it back to normal array, the internal for loop becomes faster because
    # it requires less indexing computation
    OffsetArray(restrict(A.parent, dim), newinds)
end

function restrict_eltype(A::AbstractArray)
    # infer the restrict_eltype on `eltype(A)` while preserving the container type
    # TODO: maybe there's more efficient way than the `similar` here..
    typeof(similar(A, _restrict_eltype(eltype(A)), ntuple(_->1, ndims(A))))
end
restrict_eltype(x) = _restrict_eltype(typeof(x))

for CT in (:AbstractGray, :AbstractRGB, :TransparentGray, :TransparentRGB)
    @eval _restrict_eltype(::Type{C}) where C<:$CT = __restrict_eltype(C)
end
_restrict_eltype(::Type{T}) where T = typeof(one(T)/4 + one(T)/2)
_restrict_eltype(::Type{C}) where C<:Color = __restrict_eltype(RGB{eltype(C)})
_restrict_eltype(::Type{C}) where C<:Colorant = __restrict_eltype(ARGB{eltype(C)})
function __restrict_eltype(::Type{C}) where C
    BT = base_colorant_type(C)
    isconcretetype(BT) && return floattype(BT)
    BT{_restrict_eltype(eltype(C))}
end

function restrict!(out::AbstractArray{T,N}, A::AbstractArray, dim) where {T,N}
    # a no-op for singleton dimension
    if size(A, dim) == 1 # this includes `dim > N` cases
        return copyto!(out, A)
    end
    indsout, indsA = axes(out), axes(A)
    ndims(out) == ndims(A) || throw(DimensionMismatch("input and output must have the same number of dimensions"))
    for d = 1:length(indsA)
        target = d==dim ? restrict_indices(indsA[d]) : indsA[d]
        indsout[d] == target || error("input and output must have corresponding indices; to be consistent with the input indices,\ndimension $d should be $target, got $(indsout[d])")
    end
    indspre, indspost = indsA[1:dim-1], indsA[dim+1:end]
    _restrict!(out, indsout[dim], A, indspre, indsA[dim], indspost)
end

@generated function _restrict!(out, indout, A,
                               indspre::NTuple{Npre,AbstractUnitRange},
                               indA,
                               indspost::NTuple{Npost,AbstractUnitRange}) where {Npre,Npost}
    Ipre = [Symbol(:ipre_, d) for d = 1:Npre]
    Ipost = [Symbol(:ipost_, d) for d = 1:Npost]
    quote
        $(Expr(:meta, :noinline))
        T = eltype(out)
        if isodd(length(indA))
            half = convert(eltype(T), 0.5)
            quarter = convert(eltype(T), 0.25)
            @nloops $Npost ipost d->indspost[d] begin
                iout = first(indout)
                @nloops $Npre ipre d->indspre[d] begin
                    out[$(Ipre...), iout, $(Ipost...)] = zero(T)
                end
                ispeak = true
                for iA in indA
                    if ispeak
                        @inbounds @nloops $Npre ipre d->indspre[d] begin
                            out[$(Ipre...), iout, $(Ipost...)] +=
                                half*convert(T, A[$(Ipre...), iA, $(Ipost...)])
                        end
                    else
                        @inbounds @nloops $Npre ipre d->indspre[d] begin
                            tmp = quarter*convert(T, A[$(Ipre...), iA, $(Ipost...)])
                            out[$(Ipre...), iout, $(Ipost...)]   += tmp
                            out[$(Ipre...), iout+1, $(Ipost...)] = tmp
                        end
                    end
                    ispeak = !ispeak
                    iout += ispeak
                end
            end
        else
            threeeighths = convert(eltype(T), 0.375)
            oneeighth = convert(eltype(T), 0.125)
            z = zero(T)
            fill!(out, zero(T))
            @nloops $Npost ipost d->indspost[d] begin
                peakfirst = true
                iout = first(indout)
                for iA in indA
                    if peakfirst
                        @inbounds @nloops $Npre ipre d->indspre[d] begin
                            tmp = convert(T, A[$(Ipre...), iA, $(Ipost...)])
                            out[$(Ipre...), iout, $(Ipost...)] += threeeighths*tmp
                            out[$(Ipre...), iout+1, $(Ipost...)] += oneeighth*tmp
                        end
                    else
                        @inbounds @nloops $Npre ipre d->indspre[d] begin
                            tmp = convert(T, A[$(Ipre...), iA, $(Ipost...)])
                            out[$(Ipre...), iout, $(Ipost...)]   += oneeighth*tmp
                            out[$(Ipre...), iout+1, $(Ipost...)] += threeeighths*tmp
                        end
                    end
                    peakfirst = !peakfirst
                    iout += peakfirst
                end
            end
        end
        out
    end
end

# If we're restricting along dimension 1, there are some additional efficiencies possible
@generated function _restrict!(out, indout, A, ::NTuple{0,AbstractUnitRange},
                               indA, indspost::NTuple{Npost,AbstractUnitRange}) where Npost
    Ipost = [Symbol(:ipost_, d) for d = 1:Npost]
    quote
        $(Expr(:meta, :noinline))
        T = eltype(out)
        if isodd(length(indA))
            half = convert(eltype(T), 0.5)
            quarter = convert(eltype(T), 0.25)
            @inbounds @nloops $Npost ipost d->indspost[d] begin
                iout, iA = first(indout), first(indA)
                nxt = convert(T, A[iA+1, $(Ipost...)])
                out[iout, $(Ipost...)] = half*convert(T, A[iA, $(Ipost...)]) + quarter*nxt
                for iA in first(indA)+2:2:last(indA)-2
                    prv = nxt
                    nxt = convert(T, A[iA+1, $(Ipost...)])
                    out[iout+=1, $(Ipost...)] = quarter*(prv+nxt) + half*convert(T, A[iA, $(Ipost...)])
                end
                out[iout+1, $(Ipost...)] = quarter*nxt + half*convert(T, A[last(indA), $(Ipost...)])
            end
        else
            threeeighths = convert(eltype(T), 0.375)
            oneeighth = convert(eltype(T), 0.125)
            z = zero(T)
            @inbounds @nloops $Npost ipost d->indspost[d] begin
                c = d = z
                iA = first(indA)
                for iout = first(indout):last(indout)-1
                    a, b = c, d
                    c, d = convert(T, A[iA, $(Ipost...)]), convert(T, A[iA+1, $(Ipost...)])
                    iA += 2
                    out[iout, $(Ipost...)] = oneeighth*(a+d) + threeeighths*(b+c)
                end
                out[last(indout), $(Ipost...)] = oneeighth*c + threeeighths*d
            end
        end
        out
    end
end

restrict_size(len::Integer) = isodd(len) ? (len+1)>>1 : (len>>1)+1
function restrict_indices(r::AbstractUnitRange)
    f, l = first(r), last(r)
    isodd(f) && return (f+1)>>1:restrict_size(l)
    f>>1 : (isodd(l) ? (l+1)>>1 : l>>1)
end
restrict_indices(r::Base.OneTo) = Base.OneTo(restrict_size(length(r)))
function restrict_indices(r::UnitRange)
    f, l = first(r), last(r)
    isodd(f) && return (f+1)>>1:restrict_size(l)
    f>>1 : (isodd(l) ? (l+1)>>1 : l>>1)
end
