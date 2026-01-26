struct Mul{StyleA, StyleB, AA, BB}
    A::AA
    B::BB
end

@inline Mul{StyleA,StyleB}(A::AA, B::BB) where {StyleA,StyleB,AA,BB} =
    Mul{StyleA,StyleB,AA,BB}(A,B)

@inline Mul(A, B) = Mul{typeof(MemoryLayout(A)), typeof(MemoryLayout(B))}(A, B)

@inline _mul_eltype(A) = A
@inline _mul_eltype(A, B) = Base.promote_op(*, A, B)
@inline _mul_eltype(A, B, C, D...) = _mul_eltype(Base.promote_op(*, A, B), C, D...)

@inline eltype(M::Mul) = _mul_eltype(eltype(M.A), eltype(M.B))

size(M::Mul, p::Int) = size(M)[p]
axes(M::Mul, p::Int) = axes(M)[p]
length(M::Mul) = prod(size(M))
size(M::Mul) = map(length,axes(M))

_mul_axes(A::Tuple{<:Any,<:Any}, ::Tuple{<:Any}) = (A[1],) # matrix * vector
_mul_axes(A::Tuple{<:Any}, B::Tuple{<:Any,<:Any}) = (A[1],B[2]) # vector * matrix
_mul_axes(A::Tuple{<:Any,<:Any}, B::Tuple{<:Any,<:Any}) = (A[1],B[2]) # matrix * matrix
_mul_axes(::Tuple{}, ::Tuple{}) = () # scalar * scalar
_mul_axes(::Tuple{}, B::Tuple) = B # scalar * B
_mul_axes(A::Tuple, B::Tuple{}) = A # A * scalar

axes(M::Mul) = _mul_axes(axes(M.A), axes(M.B))

# The following design is to support QuasiArrays.jl where indices
# may not be `Int`

zeroeltype(M) = zero(eltype(M)) # allow special casing where we know more about zero
zeroeltype(M::Mul{<:Any,<:Any,<:SubArray}) = zeroeltype(Mul(parent(M.A), M.B))

function _getindex(::Type{Tuple{AA}}, M::Mul, (k,)::Tuple{AA}) where AA
    A,B = M.A, M.B
    ret = zeroeltype(M)
    for j = rowsupport(A, k) ∩ colsupport(B,1)
        ret += A[k,j] * B[j]
    end
    ret
end

function _getindex(::Type{Tuple{AA,BB}}, M::Mul, (k, j)::Tuple{AA,BB}) where {AA,BB}
    A,B = M.A,M.B
    ret = zeroeltype(M)
    @inbounds for ℓ in (rowsupport(A,k) ∩ colsupport(B,j))
        ret += A[k,ℓ] * B[ℓ,j]
    end
    ret
end

# linear indexing
function _getindex(::Type{NTuple{2,Int}}, M, k::Tuple{Int})
    # convert from linear to CartesianIndex
    CartInd = CartesianIndices(axes(M))[k...]
    M[Tuple(CartInd)...]
end

_getindex(::Type{Tuple{Int}}, M, (k,)::Tuple{CartesianIndex{1}}) = M[convert(Int, k)]
_getindex(::Type{NTuple{2,Int}}, M, (kj,)::Tuple{CartesianIndex{2}}) = M[kj[1], kj[2]]

"""
   indextype(A)

gives the expected index type for an array, or array-like object. For example,
if it is vector-like it will return `Tuple{Int}`, or if it is matrix-like it will
return `Tuple{Int,Int}`. Other types may have non-integer based indexing.
"""
indextype(A) = indextype(axes(A))

indextype(axes::Tuple) = Tuple{map(eltype, axes)...}

getindex(M::Mul, k...) = _getindex(indextype(M), M, k)


"""
   mulreduce(M::Mul)

returns a lower level lazy multiplication object such as `MulAdd`, `Lmul` or `Rmul`.
The Default is `MulAdd`. Note that the lowered type must overload `copyto!` and `copy`.
"""
mulreduce(M::Mul) = MulAdd(M)

similar(M::Mul, ::Type{T}, axes) where {T} = similar(mulreduce(M), T, axes)
similar(M::Mul, ::Type{T}) where T = similar(mulreduce(M), T)
similar(M::Mul) = similar(mulreduce(M))

check_mul_axes(A) = nothing
_check_mul_axes(::Number, ::Number) = nothing
_check_mul_axes(::Number, _) = nothing
_check_mul_axes(_, ::Number) = nothing
_check_mul_axes(A, B) = axes(A, 2) == axes(B, 1) || throw_mul_axes_err(axes(A,2), axes(B,1))
@noinline function throw_mul_axes_err(axA2, axB1)
    throw(
        DimensionMismatch(
            LazyString("second axis of A, ", axA2, ", and first axis of B, ", axB1, ", must match")))
end
@noinline function throw_mul_axes_err(axA2::Base.OneTo, axB1::Base.OneTo)
    throw(
        DimensionMismatch(
            LazyString("second dimension of A, ", length(axA2), ", does not match length of x, ", length(axB1))))
end
# we need to special case AbstractQ as it allows non-compatiple multiplication
const FlexibleLeftQs = Union{QRCompactWYQ,QRPackedQ,HessenbergQ}
_check_mul_axes(::FlexibleLeftQs, ::Number) = nothing
_check_mul_axes(Q::FlexibleLeftQs, B) =
    axes(Q.factors, 1) == axes(B, 1) || axes(Q.factors, 2) == axes(B, 1) ||
        throw(DimensionMismatch(LazyString("First axis of B, ", axes(B,1), " must match either axes of A, ", axes(Q.factors))))
_check_mul_axes(::Number, ::AdjointQtype{<:Any,<:FlexibleLeftQs}) = nothing
function _check_mul_axes(A, adjQ::AdjointQtype{<:Any,<:FlexibleLeftQs})
    Q = parent(adjQ)
    axes(A, 2) == axes(Q.factors, 1) || axes(A, 2) == axes(Q.factors, 2) ||
        throw(DimensionMismatch(LazyString("Second axis of A, ", axes(A,2), " must match either axes of B, ", axes(Q.factors))))
end
_check_mul_axes(Q::FlexibleLeftQs, adjQ::AdjointQtype{<:Any,<:FlexibleLeftQs}) =
    invoke(_check_mul_axes, Tuple{Any,Any}, Q, adjQ)
function check_mul_axes(A, B, C...)
    _check_mul_axes(A, B)
    check_mul_axes(B, C...)
end

# we need to special case AbstractQ as it allows non-compatiple multiplication
function check_mul_axes(A::Union{QRCompactWYQ,QRPackedQ}, B, C...)
    axes(A.factors, 1) == axes(B, 1) || axes(A.factors, 2) == axes(B, 1) ||
        throw(DimensionMismatch(LazyString("First axis of B, ", axes(B,1), " must match either axes of A, ", axes(A))))
    check_mul_axes(B, C...)
end

@inline function instantiate(M::Mul)
    @boundscheck check_mul_axes(M.A, M.B)
    M
end

@inline materialize(M::Mul) = copy(instantiate(M))
@inline mul(A, B) = materialize(Mul(A,B))

@inline copy(M::Mul) = copy(mulreduce(M))
@inline copyto!(dest, M::Mul) = copyto!(dest, mulreduce(M))
@inline copyto!(dest::AbstractArray, M::Mul) = copyto!(dest, mulreduce(M))
mul!(dest::AbstractArray, A::AbstractArray, B::AbstractArray) = copyto!(dest, Mul(A,B))
mul!(dest::AbstractArray, A::AbstractArray, B::AbstractArray, α::Number, β::Number) = muladd!(α, A, B, β, dest)


broadcastable(M::Mul) = M

macro veclayoutmul(Typ)
    ret = quote
        (*)(A::AbstractMatrix, B::$Typ) = ArrayLayouts.mul(A,B)
        (*)(A::Adjoint{<:Any,<:AbstractMatrix{T}}, B::$Typ{S}) where {T,S} = ArrayLayouts.mul(A,B)
        (*)(A::Transpose{<:Any,<:AbstractMatrix{T}}, B::$Typ{S}) where {T,S} = ArrayLayouts.mul(A,B)
        (*)(A::LinearAlgebra.AdjointAbsVec, B::$Typ) = ArrayLayouts.mul(A,B)
        (*)(A::LinearAlgebra.TransposeAbsVec, B::$Typ) = ArrayLayouts.mul(A,B)
        (*)(A::LinearAlgebra.AdjointAbsVec{<:Number}, B::$Typ{<:Number}) = ArrayLayouts.mul(A,B)
        (*)(A::LinearAlgebra.TransposeAbsVec{T}, B::$Typ{T}) where T<:Real = ArrayLayouts.mul(A,B)
        (*)(A::LinearAlgebra.TransposeAbsVec{T,<:$Typ{T}}, B::$Typ{T}) where T<:Real = ArrayLayouts.mul(A,B)
        (*)(A::LinearAlgebra.TransposeAbsVec{T,<:$Typ{T}}, B::AbstractVector{T}) where T<:Real = ArrayLayouts.mul(A,B)

        (*)(A::LinearAlgebra.AbstractQ, B::$Typ) = ArrayLayouts.mul(A,B)
        (*)(A::$Typ, B::LinearAlgebra.LQPackedQ) = ArrayLayouts.mul(A,B)
    end
    if isdefined(LinearAlgebra, :AdjointQ)
        ret = quote
            $ret

            const FlexibleLeftQs = Union{LinearAlgebra.HessenbergQ, LinearAlgebra.QRCompactWYQ, LinearAlgebra.QRPackedQ}
            # disambiguation for flexible left-mul Qs
            (*)(A::FlexibleLeftQs, B::$Typ) = ArrayLayouts.mul(A,B)
            # flexible right-mul/adjoint left-mul Qs
            (*)(A::LinearAlgebra.AdjointQ{<:Any,<:LinearAlgebra.LQPackedQ}, B::$Typ) = ArrayLayouts.mul(A,B)
        end
    end
    for Struc in (:AbstractTriangular, :Diagonal)
        ret = quote
            $ret

            (*)(A::LinearAlgebra.$Struc, B::$Typ) = ArrayLayouts.mul(A,B)
        end
    end
    for Mod in (:AdjointAbsVec, :TransposeAbsVec)
        ret = quote
            $ret

            LinearAlgebra.mul!(dest::AbstractVector, A::$Mod{<:Any,<:$Typ}, b::AbstractVector) =
                ArrayLayouts.mul!(dest,A,b)

            (*)(A::$Mod{<:Any,<:$Typ}, B::ArrayLayouts.LayoutMatrix) = ArrayLayouts.mul(A,B)
            (*)(A::$Mod{<:Any,<:$Typ}, B::AbstractMatrix) = ArrayLayouts.mul(A,B)
            (*)(A::$Mod{<:Any,<:$Typ}, B::AbstractVector) = ArrayLayouts.mul(A,B)
            (*)(A::$Mod{<:Number,<:$Typ}, B::AbstractVector{<:Number}) = ArrayLayouts.mul(A,B)
            (*)(A::$Mod{<:Any,<:$Typ}, B::$Typ) = ArrayLayouts.mul(A,B)
            (*)(A::$Mod{<:Number,<:$Typ}, B::$Typ{<:Number}) = ArrayLayouts.mul(A,B)
            (*)(A::$Mod{<:Any,<:$Typ}, B::Diagonal) = ArrayLayouts.mul(A,B)
            (*)(A::$Mod{<:Any,<:$Typ}, B::LinearAlgebra.AbstractTriangular) = ArrayLayouts.mul(A,B)
        end
    end

    esc(ret)
end

macro layoutmul(Typ)
    ret = quote
        LinearAlgebra.mul!(dest::AbstractVector, A::$Typ, b::AbstractVector) =
            ArrayLayouts.mul!(dest,A,b)
        LinearAlgebra.mul!(dest::AbstractVector, A::$Typ, b::AbstractVector, α::Number, β::Number) =
            ArrayLayouts.mul!(dest,A,b,α,β)

        LinearAlgebra.mul!(dest::AbstractMatrix, A::$Typ, B::AbstractMatrix) =
            ArrayLayouts.mul!(dest,A,B)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::AbstractMatrix, B::$Typ) =
            ArrayLayouts.mul!(dest,A,B)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::$Typ, B::$Typ) =
            ArrayLayouts.mul!(dest,A,B)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::$Typ, B::AbstractMatrix, α::Number, β::Number) =
            ArrayLayouts.mul!(dest,A,B,α,β)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::AbstractMatrix, B::$Typ, α::Number, β::Number) =
            ArrayLayouts.mul!(dest,A,B,α,β)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::$Typ, B::$Typ, α::Number, β::Number) =
            ArrayLayouts.mul!(dest,A,B,α,β)

        (*)(A::$Typ, B::$Typ) = ArrayLayouts.mul(A,B)
        (*)(A::$Typ, B::AbstractMatrix) = ArrayLayouts.mul(A,B)
        (*)(A::$Typ, B::AbstractVector) = ArrayLayouts.mul(A,B)
        (*)(A::$Typ, B::ArrayLayouts.LayoutVector) = ArrayLayouts.mul(A,B)
        (*)(A::AbstractMatrix, B::$Typ) = ArrayLayouts.mul(A,B)
        (*)(A::LinearAlgebra.AdjointAbsVec, B::$Typ) = ArrayLayouts.mul(A,B)
        (*)(A::LinearAlgebra.TransposeAbsVec, B::$Typ) = ArrayLayouts.mul(A,B)
        (*)(A::LinearAlgebra.AdjointAbsVec{<:Any,<:Zeros{<:Any,1}}, B::$Typ) = ArrayLayouts.mul(A,B)
        (*)(A::LinearAlgebra.TransposeAbsVec{<:Any,<:Zeros{<:Any,1}}, B::$Typ) = ArrayLayouts.mul(A,B)
        (*)(A::LinearAlgebra.Transpose{T,<:$Typ}, B::Zeros{T,1}) where T<:Real = ArrayLayouts.mul(A,B)

        (*)(A::LinearAlgebra.AbstractQ, B::$Typ) = ArrayLayouts.mul(A,B)
        (*)(A::$Typ, B::LinearAlgebra.AbstractQ) = ArrayLayouts.mul(A,B)
    end
    if isdefined(LinearAlgebra, :AdjointQ)
        ret = quote
            $ret

            const FlexibleLeftQs = Union{LinearAlgebra.HessenbergQ, LinearAlgebra.QRCompactWYQ, LinearAlgebra.QRPackedQ}
            # disambiguation for flexible left-mul/adjoint right-mul Qs
            (*)(A::FlexibleLeftQs, B::$Typ) = ArrayLayouts.mul(A,B)
            (*)(A::$Typ, B::LinearAlgebra.AdjointQ{<:Any,<:FlexibleLeftQs}) = ArrayLayouts.mul(A,B)
            # disambiguation for flexible right-mul/adjoint left-mul Qs
            (*)(A::$Typ, B::LinearAlgebra.LQPackedQ) = ArrayLayouts.mul(A,B)
            (*)(A::LinearAlgebra.AdjointQ{<:Any,<:LinearAlgebra.LQPackedQ}, B::$Typ) = ArrayLayouts.mul(A,B)
        end
    end
    for Struc in (:AbstractTriangular, :Diagonal, :Bidiagonal, :SymTridiagonal, :Tridiagonal)
        # starting from Julia v1.10, the last four could be put into a single Union to
        # reduce the number of mul! methods; or perhaps addressed as some common supertype
        ret = quote
            $ret

            LinearAlgebra.mul!(dest::AbstractMatrix, A::$Typ, B::LinearAlgebra.$Struc, α::Number, β::Number) =
                ArrayLayouts.mul!(dest,A,B,α,β)
            LinearAlgebra.mul!(dest::AbstractMatrix, A::LinearAlgebra.$Struc, B::$Typ, α::Number, β::Number) =
                ArrayLayouts.mul!(dest,A,B,α,β)

            (*)(A::LinearAlgebra.$Struc, B::$Typ) = ArrayLayouts.mul(A,B)
            (*)(A::$Typ, B::LinearAlgebra.$Struc) = ArrayLayouts.mul(A,B)
        end
    end
    for Mod in (:Adjoint, :Transpose, :Symmetric, :Hermitian)
        ret = quote
            $ret

            LinearAlgebra.mul!(dest::AbstractMatrix, A::$Typ, b::$Mod{<:Any,<:AbstractMatrix}) =
                ArrayLayouts.mul!(dest,A,b)
            LinearAlgebra.mul!(dest::AbstractVector, A::$Mod{<:Any,<:$Typ}, b::AbstractVector) =
                ArrayLayouts.mul!(dest,A,b)
            LinearAlgebra.mul!(dest::AbstractVector, A::$Mod{<:Any,<:$Typ}, b::AbstractVector, α::Number, β::Number) =
                ArrayLayouts.mul!(dest,A,b,α,β)
            LinearAlgebra.mul!(dest::AbstractMatrix, A::$Mod{<:Any,<:$Typ}, B::AbstractMatrix, α::Number, β::Number) =
                ArrayLayouts.mul!(dest,A,B,α,β)
            LinearAlgebra.mul!(dest::AbstractMatrix, A::AbstractMatrix, B::$Mod{<:Any,<:$Typ}, α::Number, β::Number) =
                ArrayLayouts.mul!(dest,A,B,α,β)
            LinearAlgebra.mul!(dest::AbstractMatrix, A::$Mod{<:Any,<:$Typ}, B::$Typ, α::Number, β::Number) =
                ArrayLayouts.mul!(dest,A,B,α,β)
            LinearAlgebra.mul!(dest::AbstractMatrix, A::$Typ, B::$Mod{<:Any,<:$Typ}, α::Number, β::Number) =
                ArrayLayouts.mul!(dest,A,B,α,β)
            LinearAlgebra.mul!(dest::AbstractMatrix, A::$Mod{<:Any,<:AbstractVecOrMat}, B::$Typ, α::Number, β::Number) =
                ArrayLayouts.mul!(dest,A,B,α,β)
            LinearAlgebra.mul!(dest::AbstractMatrix, A::$Typ, B::$Mod{<:Any,<:AbstractVecOrMat}, α::Number, β::Number) =
                ArrayLayouts.mul!(dest,A,B,α,β)
            LinearAlgebra.mul!(dest::AbstractMatrix, A::$Mod{<:Any,<:$Typ}, B::$Mod{<:Any,<:$Typ}, α::Number, β::Number) =
                ArrayLayouts.mul!(dest,A,B,α,β)

            (*)(A::$Mod{<:Any,<:$Typ}, B::$Mod{<:Any,<:$Typ}) = ArrayLayouts.mul(A,B)
            (*)(A::$Mod{<:Any,<:$Typ}, B::AbstractMatrix) = ArrayLayouts.mul(A,B)
            (*)(A::AbstractMatrix, B::$Mod{<:Any,<:$Typ}) = ArrayLayouts.mul(A,B)
            (*)(A::LinearAlgebra.AdjointAbsVec, B::$Mod{<:Any,<:$Typ}) = ArrayLayouts.mul(A,B)
            (*)(A::LinearAlgebra.TransposeAbsVec, B::$Mod{<:Any,<:$Typ}) = ArrayLayouts.mul(A,B)
            (*)(A::$Mod{<:Any,<:$Typ}, B::AbstractVector) = ArrayLayouts.mul(A,B)
            (*)(A::$Mod{<:Any,<:$Typ}, B::ArrayLayouts.LayoutVector) = ArrayLayouts.mul(A,B)
            (*)(A::$Mod{<:Any,<:$Typ}, B::Zeros{<:Any,1}) = ArrayLayouts.mul(A,B)

            (*)(A::$Mod{<:Any,<:$Typ}, B::$Typ) = ArrayLayouts.mul(A,B)
            (*)(A::$Typ, B::$Mod{<:Any,<:$Typ}) = ArrayLayouts.mul(A,B)

            (*)(A::$Mod{<:Any,<:$Typ}, B::Diagonal) = ArrayLayouts.mul(A,B)
            (*)(A::Diagonal, B::$Mod{<:Any,<:$Typ}) = ArrayLayouts.mul(A,B)

            (*)(A::LinearAlgebra.AbstractTriangular, B::$Mod{<:Any,<:$Typ}) = ArrayLayouts.mul(A,B)
            (*)(A::$Mod{<:Any,<:$Typ}, B::LinearAlgebra.AbstractTriangular) = ArrayLayouts.mul(A,B)
        end
    end

    ret = quote
        $ret

        LinearAlgebra.mul!(dest::AbstractMatrix, A::$Typ, b::UpperOrLowerTriangular) =
            ArrayLayouts.mul!(dest,A,b)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::UpperOrLowerTriangular, b::$Typ) =
            ArrayLayouts.mul!(dest,A,b)
        LinearAlgebra.mul!(dest::AbstractVector, A::UpperOrLowerTriangular{<:Any,<:$Typ}, b::AbstractVector) =
            ArrayLayouts.mul!(dest,A,b)
        LinearAlgebra.mul!(dest::AbstractVector, A::UpperOrLowerTriangular{<:Any,<:$Typ}, b::AbstractVector, α::Number, β::Number) =
            ArrayLayouts.mul!(dest,A,b,α,β)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::UpperOrLowerTriangular{<:Any,<:$Typ}, B::AbstractMatrix, α::Number, β::Number) =
            ArrayLayouts.mul!(dest,A,B,α,β)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::UpperOrLowerTriangular{<:Any,<:$Typ}, B::UpperOrLowerTriangular, α::Number, β::Number) =
            ArrayLayouts.mul!(dest,A,B,α,β)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::AbstractMatrix, B::UpperOrLowerTriangular{<:Any,<:$Typ}, α::Number, β::Number) =
            ArrayLayouts.mul!(dest,A,B,α,β)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::UpperOrLowerTriangular, B::UpperOrLowerTriangular{<:Any,<:$Typ}, α::Number, β::Number) =
            ArrayLayouts.mul!(dest,A,B,α,β)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::UpperOrLowerTriangular{<:Any,<:$Typ}, B::$Typ, α::Number, β::Number) =
            ArrayLayouts.mul!(dest,A,B,α,β)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::$Typ, B::UpperOrLowerTriangular{<:Any,<:$Typ}, α::Number, β::Number) =
            ArrayLayouts.mul!(dest,A,B,α,β)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::UpperOrLowerTriangular, B::$Typ, α::Number, β::Number) =
            ArrayLayouts.mul!(dest,A,B,α,β)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::$Typ, B::UpperOrLowerTriangular, α::Number, β::Number) =
            ArrayLayouts.mul!(dest,A,B,α,β)
        LinearAlgebra.mul!(dest::AbstractMatrix, A::UpperOrLowerTriangular{<:Any,<:$Typ}, B::UpperOrLowerTriangular{<:Any,<:$Typ}, α::Number, β::Number) =
            ArrayLayouts.mul!(dest,A,B,α,β)

    end
    esc(ret)
end

@veclayoutmul LayoutVector
*(A::Adjoint{<:Any,<:LayoutVector}, B::Adjoint{<:Any,<:LayoutMatrix}) = mul(A,B)
*(A::Adjoint{<:Any,<:LayoutVector}, B::Transpose{<:Any,<:LayoutMatrix}) = mul(A,B)
*(A::Transpose{<:Any,<:LayoutVector}, B::Adjoint{<:Any,<:LayoutMatrix}) = mul(A,B)
*(A::Transpose{<:Any,<:LayoutVector}, B::Transpose{<:Any,<:LayoutMatrix}) = mul(A,B)

# Disambiguation with FillArrays
*(A::AbstractFill{<:Any,2}, B::LayoutVector) = invoke(*, Tuple{AbstractFill{<:Any,2}, AbstractVector}, A, B)
*(A::Adjoint{<:Any, <:LayoutVector}, B::AbstractFill{<:Any,2}) = invoke(*, Tuple{Adjoint{<:Any, <:AbstractVector}, AbstractFill{<:Any,2}}, A, B)
*(A::Transpose{<:Any, <:LayoutVector}, B::AbstractFill{<:Any,2}) = invoke(*, Tuple{Transpose{<:Any, <:AbstractVector}, AbstractFill{<:Any,2}}, A, B)

## special routines introduced in v0.9. We need to avoid these to support ∞-arrays

*(x::Adjoint{<:Any,<:LayoutVector},   D::Diagonal{<:Any,<:LayoutVector}) = mul(x, D)
*(x::Transpose{<:Any,<:LayoutVector},   D::Diagonal{<:Any,<:LayoutVector}) = mul(x, D)
*(x::AdjointAbsVec,   D::Diagonal, y::LayoutVector) = x * mul(D,y)
*(x::TransposeAbsVec, D::Diagonal, y::LayoutVector) = x * mul(D,y)
*(x::AdjointAbsVec{<:Any,<:Zeros{<:Any,1}},   D::Diagonal, y::LayoutVector) = FillArrays._triple_zeromul(x, D, y)
*(x::TransposeAbsVec{<:Any,<:Zeros{<:Any,1}}, D::Diagonal, y::LayoutVector) = FillArrays._triple_zeromul(x, D, y)


*(A::UpperOrLowerTriangular{<:Any,<:LayoutMatrix}, B::UpperOrLowerTriangular{<:Any,<:LayoutMatrix}) = mul(A, B)
*(A::UpperOrLowerTriangular{<:Any,<:AdjOrTrans{<:Any,<:LayoutMatrix}}, B::UpperOrLowerTriangular{<:Any,<:LayoutMatrix}) = mul(A, B)
*(A::UpperOrLowerTriangular{<:Any,<:LayoutMatrix}, B::UpperOrLowerTriangular{<:Any,<:AdjOrTrans{<:Any,<:LayoutMatrix}}) = mul(A, B)
*(A::UpperOrLowerTriangular{<:Any,<:AdjOrTrans{<:Any,<:LayoutMatrix}}, B::UpperOrLowerTriangular{<:Any,<:AdjOrTrans{<:Any,<:LayoutMatrix}}) = mul(A, B)

*(A::SymTridiagonal{<:Any,<:LayoutVector}, B::SymTridiagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::SymTridiagonal{<:Any,<:LayoutVector}, B::Tridiagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::SymTridiagonal{<:Any,<:LayoutVector}, B::Bidiagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::SymTridiagonal{<:Any,<:LayoutVector}, B::UpperOrLowerTriangular{<:Any,<:LayoutMatrix}) = mul(A, B)
*(A::SymTridiagonal, B::Diagonal{<:Any,<:LayoutVector}) = mul(A, B) # ambiguity
*(A::Tridiagonal{<:Any,<:LayoutVector}, B::Tridiagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::Tridiagonal{<:Any,<:LayoutVector}, B::SymTridiagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::Tridiagonal{<:Any,<:LayoutVector}, B::Bidiagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::Tridiagonal{<:Any,<:LayoutVector}, B::UpperOrLowerTriangular{<:Any,<:LayoutMatrix}) = mul(A, B)
*(A::Bidiagonal{<:Any,<:LayoutVector}, B::Bidiagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::Bidiagonal{<:Any,<:LayoutVector}, B::SymTridiagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::Bidiagonal{<:Any,<:LayoutVector}, B::Tridiagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::Bidiagonal{<:Any,<:LayoutVector}, B::UpperOrLowerTriangular{<:Any,<:LayoutMatrix}) = mul(A, B)
*(A::UpperOrLowerTriangular{<:Any,<:LayoutMatrix}, B::SymTridiagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::UpperOrLowerTriangular{<:Any,<:LayoutMatrix}, B::Tridiagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::UpperOrLowerTriangular{<:Any,<:LayoutMatrix}, B::Bidiagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::UpperOrLowerTriangular{<:Any,<:LayoutMatrix}, B::Diagonal{<:Any,<:LayoutVector}) = mul(A, B) # ambiguity
*(A::Diagonal{<:Any,<:LayoutVector}, B::SymTridiagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::Diagonal{<:Any,<:LayoutVector}, B::UpperOrLowerTriangular{<:Any,<:LayoutMatrix}) = mul(A, B)

# mul! for subarray of layout matrix
const SubLayoutMatrix = Union{SubArray{<:Any,2,<:LayoutMatrix}, SubArray{<:Any,2,<:AdjOrTrans{<:Any,<:LayoutMatrix}}}

*(A::Diagonal, B::SubLayoutMatrix) = mul(A, B)
*(A::SubLayoutMatrix, B::Diagonal) = mul(A, B)

LinearAlgebra.mul!(C::AbstractMatrix, A::SubLayoutMatrix, B::AbstractMatrix, α::Number, β::Number) =
    ArrayLayouts.mul!(C, A, B, α, β)
LinearAlgebra.mul!(C::AbstractMatrix, A::AbstractMatrix, B::SubLayoutMatrix, α::Number, β::Number) =
    ArrayLayouts.mul!(C, A, B, α, β)
LinearAlgebra.mul!(C::AbstractMatrix, A::Diagonal, B::SubLayoutMatrix, α::Number, β::Number) =
    ArrayLayouts.mul!(C, A, B, α, β)
LinearAlgebra.mul!(C::AbstractMatrix, A::SubLayoutMatrix, B::Diagonal, α::Number, β::Number) =
    ArrayLayouts.mul!(C, A, B, α, β)
LinearAlgebra.mul!(C::AbstractMatrix, A::SubLayoutMatrix, B::LayoutMatrix, α::Number, β::Number) =
    ArrayLayouts.mul!(C, A, B, α, β)
LinearAlgebra.mul!(C::AbstractMatrix, A::LayoutMatrix, B::SubLayoutMatrix, α::Number, β::Number) =
    ArrayLayouts.mul!(C, A, B, α, β)    
LinearAlgebra.mul!(C::AbstractMatrix, A::SubLayoutMatrix, B::SubLayoutMatrix, α::Number, β::Number) =
    ArrayLayouts.mul!(C, A, B, α, β)
LinearAlgebra.mul!(C::AbstractVector, A::SubLayoutMatrix, B::AbstractVector, α::Number, β::Number) =
    ArrayLayouts.mul!(C, A, B, α, β)


###
# Dot
###
"""
    Dot(A, B)

is a lazy version of `dot(A, B)`, designed to support
materializing based on `MemoryLayout`.
"""
struct Dot{StyleA,StyleB,ATyp,BTyp}
    A::ATyp
    B::BTyp
end

"""
    Dotu(A, B)

is a lazy version of `BLAS.dotu(A, B)`, designed to support
materializing based on `MemoryLayout`.
"""
struct Dotu{StyleA,StyleB,ATyp,BTyp}
    A::ATyp
    B::BTyp
end



for Dt in (:Dot, :Dotu)
    @eval begin
        @inline $Dt(A::ATyp,B::BTyp) where {ATyp,BTyp} = $Dt{typeof(MemoryLayout(ATyp)), typeof(MemoryLayout(BTyp)), ATyp, BTyp}(A, B)
        @inline materialize(d::$Dt) = copy(instantiate(d))
        @inline eltype(D::$Dt) = promote_type(eltype(D.A), eltype(D.B))
    end
end
@inline copy(d::Dot) = invoke(LinearAlgebra.dot, Tuple{AbstractArray,AbstractArray}, d.A, d.B)
@inline copy(d::Dotu{<:AbstractStridedLayout,<:AbstractStridedLayout,<:AbstractVector{T},<:AbstractVector{T}}) where T <: BlasComplex = BLAS.dotu(d.A, d.B)
@inline copy(d::Dotu{<:AbstractStridedLayout,<:AbstractStridedLayout,<:AbstractVector{T},<:AbstractVector{T}}) where T <: BlasReal = BLAS.dot(d.A, d.B)
@inline copy(d::Dotu) = LinearAlgebra._dot_nonrecursive(d.A, d.B)

@inline Dot(M::Mul{<:DualLayout,<:Any,<:AbstractMatrix,<:AbstractVector}) = Dot(M.A', M.B)
@inline Dotu(M::Mul{<:DualLayout,<:Any,<:AbstractMatrix,<:AbstractVector}) = Dotu(transpose(M.A), M.B)

@inline _dot_or_dotu(::typeof(transpose), ::Type{<:Complex}, M) = Dotu(M)
@inline _dot_or_dotu(_, _, M) = Dot(M)
@inline mulreduce(M::Mul{<:DualLayout,<:Any,<:AbstractMatrix,<:AbstractVector}) = _dot_or_dotu(dualadjoint(M.A), eltype(M.A), M)



dot(a, b) = materialize(Dot(a, b))
dotu(a, b) = materialize(Dotu(a, b))

@inline LinearAlgebra.dot(a::LayoutArray, b::LayoutArray) = dot(a,b)
@inline LinearAlgebra.dot(a::LayoutArray, b::AbstractArray) = dot(a,b)
@inline LinearAlgebra.dot(a::AbstractArray, b::LayoutArray) = dot(a,b)
@inline LinearAlgebra.dot(a::LayoutVector, b::AbstractFill{<:Any,1}) = FillArrays._fill_dot_rev(a,b)
@inline LinearAlgebra.dot(a::AbstractFill{<:Any,1}, b::LayoutVector) = FillArrays._fill_dot(a,b)

@inline LinearAlgebra.dot(a::SubArray{<:Any,N,<:LayoutArray}, b::AbstractArray) where N = dot(a,b)
@inline LinearAlgebra.dot(a::SubArray{<:Any,N,<:LayoutArray}, b::LayoutArray) where N = dot(a,b)
@inline LinearAlgebra.dot(a::AbstractArray, b::SubArray{<:Any,N,<:LayoutArray}) where N = dot(a,b)
@inline LinearAlgebra.dot(a::LayoutArray, b::SubArray{<:Any,N,<:LayoutArray}) where N = dot(a,b)
@inline LinearAlgebra.dot(a::SubArray{<:Any,N,<:LayoutArray}, b::SubArray{<:Any,N,<:LayoutArray}) where N = dot(a,b)

# Temporary until layout 3-arg dot is added.
# We go to generic fallback as layout-arrays are structured
dot(x, A, y) = dot(x, mul(A, y))
LinearAlgebra.dot(x::AbstractVector, A::LayoutMatrix, y::AbstractVector) = dot(x, A, y)
LinearAlgebra.dot(x::AbstractVector, A::Symmetric{<:Real,<:LayoutMatrix}, y::AbstractVector) = dot(x, A, y)



# allow overloading for infinite or lazy case
@inline _power_by_squaring(_, _, A, p) = invoke(Base.power_by_squaring, Tuple{AbstractMatrix,Integer}, A, p)
# TODO: Remove unnecessary _apply
_apply(_, _, op, Λ::UniformScaling, A::AbstractMatrix) = op(Diagonal(Fill(Λ.λ,(axes(A,1),))), A)
_apply(_, _, op, A::AbstractMatrix, Λ::UniformScaling) = op(A, Diagonal(Fill(Λ.λ,(axes(A,1),))))

for Typ in (:LayoutMatrix, :(Symmetric{<:Any,<:LayoutMatrix}), :(Hermitian{<:Any,<:LayoutMatrix}),
            :(Adjoint{<:Any,<:LayoutMatrix}), :(Transpose{<:Any,<:LayoutMatrix}))
    @eval begin
        @inline Base.power_by_squaring(A::$Typ, p::Integer) = _power_by_squaring(MemoryLayout(A), size(A), A, p)
        @inline +(A::$Typ, Λ::UniformScaling) = _apply(MemoryLayout(A), size(A), +, A, Λ)
        @inline +(Λ::UniformScaling, A::$Typ) = _apply(MemoryLayout(A), size(A), +, Λ, A)
        @inline -(A::$Typ, Λ::UniformScaling) = _apply(MemoryLayout(A), size(A), -, A, Λ)
        @inline -(Λ::UniformScaling, A::$Typ) = _apply(MemoryLayout(A), size(A), -, Λ, A)
    end
end
