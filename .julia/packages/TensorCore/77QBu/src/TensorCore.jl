module TensorCore

using LinearAlgebra

export ⊙, hadamard, hadamard!
export ⊗, tensor, tensor!
export ⊡, boxdot, boxdot!

"""
    hadamard(a, b)
    a ⊙ b

For arrays `a` and `b`, perform elementwise multiplication.
`a` and `b` must have identical `axes`.

`⊙` can be passed as an operator to higher-order functions.

# Examples
```jldoctest; setup=:(using TensorCore)
julia> a = [2, 3]; b = [5, 7];

julia> a ⊙ b
2-element Array{$Int,1}:
 10
 21

julia> a ⊙ [5]
ERROR: DimensionMismatch("Axes of `A` and `B` must match, got (Base.OneTo(2),) and (Base.OneTo(1),)")
[...]
```

See also `hadamard!(y, a, b)`.
"""
function hadamard(A::AbstractArray, B::AbstractArray)
    @noinline throw_dmm(axA, axB) = throw(DimensionMismatch("Axes of `A` and `B` must match, got $axA and $axB"))

    axA, axB = axes(A), axes(B)
    axA == axB || throw_dmm(axA, axB)
    return map(*, A, B)
end
const ⊙ = hadamard

"""
    hadamard!(dest, A, B)

Similar to `hadamard(A, B)` (which can also be written `A ⊙ B`), but stores its results in
the pre-allocated array `dest`.
"""
function hadamard!(dest::AbstractArray, A::AbstractArray, B::AbstractArray)
    @noinline function throw_dmm(axA, axB, axdest)
        throw(DimensionMismatch("`axes(dest) = $axdest` must be equal to `axes(A) = $axA` and `axes(B) = $axB`"))
    end

    axA, axB, axdest = axes(A), axes(B), axes(dest)
    ((axdest == axA) & (axdest == axB)) || throw_dmm(axA, axB, axdest)
    @simd for I in eachindex(dest, A, B)
        @inbounds dest[I] = A[I] * B[I]
    end
    return dest
end

"""
    tensor(A, B)
    A ⊗ B

Compute the tensor product of `A` and `B`.
If `C = A ⊗ B`, then `C[i1, ..., im, j1, ..., jn] = A[i1, ... im] * B[j1, ..., jn]`.

For vectors `v` and `w`, the Kronecker product is related to the tensor product by
`kron(v,w) == vec(w ⊗ v)` or `w ⊗ v == reshape(kron(v,w), (length(w), length(v)))`.

# Examples
```jldoctest; setup=:(using TensorCore)
julia> a = [2, 3]; b = [5, 7, 11];

julia> a ⊗ b
2×3 Array{$Int,2}:
 10  14  22
 15  21  33
```
See also `tensor!(Y,A,B)`.
"""
tensor(A::AbstractArray, B::AbstractArray) = [a*b for a in A, b in B]
const ⊗ = tensor

const CovectorLike{T} = Union{Adjoint{T,<:AbstractVector},Transpose{T,<:AbstractVector}}
function tensor(u::AbstractArray, v::CovectorLike)
    # If `v` is thought of as a covector, you might want this to be two-dimensional,
    # but thought of as a matrix it should be three-dimensional.
    # The safest is to avoid supporting it at all. See discussion in #35150.
    error("`tensor` is not defined for co-vectors, perhaps you meant `*`?")
end
function tensor(u::CovectorLike, v::AbstractArray)
    error("`tensor` is not defined for co-vectors, perhaps you meant `*`?")
end
function tensor(u::CovectorLike, v::CovectorLike)
    error("`tensor` is not defined for co-vectors, perhaps you meant `*`?")
end

"""
    tensor!(dest, A, B)

Similar to `tensor(A, B)` (which can also be written `A ⊗ B`), but stores its results in
the pre-allocated array `dest`.
"""
function tensor!(dest::AbstractArray, A::AbstractArray, B::AbstractArray)
    @noinline function throw_dmm(axA, axB, axdest)
        throw(DimensionMismatch("`axes(dest) = $axdest` must concatenate `axes(A) = $axA` and `axes(B) = $axB`"))
    end

    axA, axB, axdest = axes(A), axes(B), axes(dest)
    axes(dest) == (axA..., axB...) || throw_dmm(axA, axB, axdest)
    if IndexStyle(dest) === IndexCartesian()
        for IB in CartesianIndices(axB)
            @inbounds b = B[IB]
            @simd for IA in CartesianIndices(axA)
                @inbounds dest[IA,IB] = A[IA]*b
            end
        end
    else
        i = firstindex(dest)
        @inbounds for b in B
            @simd for a in A
                dest[i] = a*b
                i += 1
            end
        end
    end
    return dest
end

export boxdot, ⊡, boxdot!

"""
    boxdot(A,B) = A ⊡ B    # \\boxdot

Generalised matrix multiplication: Contracts the last dimension of `A` with
the first dimension of `B`, for any `ndims(A)` & `ndims(B)`.
If both are vectors, then it returns a scalar `== sum(A .* B)`.

# Examples
```jldoctest; setup=:(using TensorCore)
julia> A = rand(3,4,5); B = rand(5,6,7);

julia> size(A ⊡ B)
(3, 4, 6, 7)

julia> typeof(rand(5) ⊡ rand(5))
Float64

julia> try B ⊡ A catch err println(err) end
DimensionMismatch("neighbouring axes of `A` and `B` must match, got Base.OneTo(7) and Base.OneTo(3)")
```
This is the same behaviour as Mathematica's function `Dot[A, B]`.
It is not identicaly to Python's `numpy.dot(A, B)`, which contracts with the second-last
dimension of `B` instead of the first, but both keep all the other dimensions.
Unlike Julia's `LinearAlgebra.dot`, it does not conjugate `A`, so these two agree only
for real-valued vectors.

When interacting with `Adjoint` vectors, this always obeys `(x ⊡ y)' == y' ⊡ x'`,
and hence may sometimes return another `Adjoint` vector. (And similarly for `Transpose`.)

```jldoctest; setup=:(using TensorCore)
julia> M = rand(5,5); v = rand(5);

julia> typeof(v ⊡ M')
Array{Float64,1}

julia> typeof(M ⊡ v')  # adjoint of the previous line
Adjoint{Float64,Array{Float64,1}}

julia> typeof(v' ⊡ M')  # same as *, and equal to adjoint(M ⊡ v)
Adjoint{Float64,Array{Float64,1}}

julia> typeof(v' ⊡ v)
Float64
```
See also `boxdot!(Y,A,B)`, which is to `⊡` as `mul!` is to `*`.
"""
function boxdot(A::AbstractArray, B::AbstractArray)
    Amat = _squash_left(A)
    Bmat = _squash_right(B)

    axA, axB = axes(Amat,2), axes(Bmat,1)
    axA == axB || _throw_dmm(axA, axB)

    return _boxdot_reshape(Amat * Bmat, A, B)
end

const ⊡ = boxdot

@noinline _throw_dmm(axA, axB) = throw(DimensionMismatch("neighbouring axes of `A` and `B` must match, got $axA and $axB"))

_squash_left(A::AbstractArray) = reshape(A, :,size(A,ndims(A)))
_squash_left(A::AbstractMatrix) = A

_squash_right(B::AbstractArray) = reshape(B, size(B,1),:)
_squash_right(B::AbstractVecOrMat) = B

function _boxdot_reshape(AB::AbstractArray, A::AbstractArray{T,N}, B::AbstractArray{S,M}) where {T,N,S,M}
    ax = ntuple(i -> i<N ? axes(A, i) : axes(B, i-N+2), Val(N+M-2))
    reshape(AB, ax) # some cases don't come here, so this doesn't really support OffsetArrays
end

# These can skip final reshape:
_boxdot_reshape(AB::AbstractVecOrMat, A::AbstractMatrix, B::AbstractVecOrMat) = AB

# These produce scalar output:
function boxdot(A::AbstractVector, B::AbstractVector)
    axA, axB = axes(A,1), axes(B,1)
    axA == axB || _throw_dmm(axA, axB)
    if eltype(A) <: Number
        return transpose(A)*B
    else
        return sum(a*b for (a,b) in zip(A,B))
    end
end

# Multiplication by a scalar:
boxdot(A::AbstractArray, b::Number) = A*b
boxdot(a::Number, B::AbstractArray) = a*B
boxdot(a::Number, b::Number) = a*b

using LinearAlgebra: AdjointAbsVec, TransposeAbsVec, AdjOrTransAbsVec

# Adjont and Transpose, vectors or almost (returning a scalar)
boxdot(A::AdjointAbsVec, B::AbstractVector) = A * B
boxdot(A::TransposeAbsVec, B::AbstractVector) = A * B

boxdot(A::AbstractVector, B::AdjointAbsVec) = A ⊡ vec(B)
boxdot(A::AbstractVector, B::TransposeAbsVec) = A ⊡ vec(B)

boxdot(A::AdjointAbsVec, B::AdjointAbsVec) = adjoint(adjoint(B) ⊡ adjoint(A))
boxdot(A::AdjointAbsVec, B::TransposeAbsVec) = vec(A) ⊡ vec(B)
boxdot(A::TransposeAbsVec, B::AdjointAbsVec) = vec(A) ⊡ vec(B)
boxdot(A::TransposeAbsVec, B::TransposeAbsVec) = transpose(transpose(B) ⊡ transpose(A))

# ... with a matrix (returning another such)
boxdot(A::AdjointAbsVec, B::AbstractMatrix) = A * B
boxdot(A::TransposeAbsVec, B::AbstractMatrix) = A * B

boxdot(A::AbstractMatrix, B::AdjointAbsVec) = (B' ⊡ A')'
boxdot(A::AbstractMatrix, B::TransposeAbsVec) = transpose(transpose(B) ⊡ transpose(A))

# ... and with higher-dim (returning a plain array)
boxdot(A::AdjointAbsVec, B::AbstractArray) = vec(A) ⊡ B
boxdot(A::TransposeAbsVec, B::AbstractArray) = vec(A) ⊡ B

boxdot(A::AbstractArray, B::AdjointAbsVec) = A ⊡ vec(B)
boxdot(A::AbstractArray, B::TransposeAbsVec) = A ⊡ vec(B)


"""
    boxdot!(Y, A, B, α=1, β=0)

In-place version of `boxdot`, i.e. `Y .= (A ⊡ B) .* β .+ Y .* α`.
Like 5-argument `mul!`, the use of `α, β` here requires Julia 1.3 or later.
"""
function boxdot! end

if VERSION < v"1.3" # Then 5-arg mul! isn't defined

    function boxdot!(Y::AbstractArray, A::AbstractArray, B::AbstractArray)
        szY = prod(size(A)[1:end-1]), prod(size(B)[2:end])
        mul!(reshape(Y, szY), _squash_left(A), _squash_right(B))
        Y
    end

    boxdot!(Y::AbstractArray, A::AbstractArray, B::AdjOrTransAbsVec) = boxdot!(Y, A, vec(B))

else

    function boxdot!(Y::AbstractArray, A::AbstractArray, B::AbstractArray, α::Number=true, β::Number=false)
        szY = prod(size(A)[1:end-1]), prod(size(B)[2:end])
        mul!(reshape(Y, szY), _squash_left(A), _squash_right(B), α, β)
        Y
    end

    # For boxdot!, only where mul! behaves differently:
    boxdot!(Y::AbstractArray, A::AbstractArray, B::AdjOrTransAbsVec,
        α::Number=true, β::Number=false) = boxdot!(Y, A, vec(B))

end

"""
    TensorCore._adjoint(A)

This extends `adjoint` to understand higher-dimensional arrays, always reversing the
order of dimensions. On Julia 1.5 and later, the symbol `'` can be overloaded locally
as `var"'"`, as shown below.

Then `(x ⊡ y)' == y' ⊡ x'` holds for `x` and `y` arrays of any dimension.

# Examples
```jldoctest; setup=:(using TensorCore)
julia> T3 = rand(3,4,5); v = rand(5);

julia> size(T3 ⊡ v')
(3, 4)

julia> let var"'" = TensorCore._adjoint
         v ⊡ T3' ≈ (T3 ⊡ v')'
       end
true
```
"""
_adjoint(x) = adjoint(x)
_adjoint(x::AbstractVecOrMat) = adjoint(x)
_adjoint(x::AbstractArray{T,N}) where {T<:Number,N} = conj(PermutedDimsArray(x, ntuple(i -> N-i+1, N)))
_adjoint(x::AbstractArray{T,N}) where {T,N} = adjoint.(PermutedDimsArray(x, ntuple(i -> N-i+1, N)))

_transpose(x) = transpose(x)
_transpose(x::AbstractVecOrMat) = transpose(x)
_transpose(x::AbstractArray{T,N}) where {T<:Number,N} = PermutedDimsArray(x, ntuple(i -> N-i+1, N))
_transpose(x::AbstractArray{T,N}) where {T,N} = transpose.(PermutedDimsArray(x, ntuple(i -> N-i+1, N)))

end
