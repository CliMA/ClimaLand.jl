module WoodburyMatrices

using LinearAlgebra
import LinearAlgebra: det, logdet, logabsdet, ldiv!, mul!, adjoint, transpose, diag, issymmetric
import Base: +, *, \, inv, /, convert, copy, show, similar, axes, size
using SparseArrays

export AbstractWoodbury, Woodbury, SymWoodbury

abstract type AbstractWoodbury{T} <: Factorization{T} end

safeinv(A) = inv(A)
safeinv(A::SparseMatrixCSC) = safeinv(Matrix(A))

include("woodbury.jl")
include("symwoodbury.jl")
include("sparsefactors.jl")

# Traits and algorithms expressible in terms of AbstractWoodbury

size(W::AbstractWoodbury) = size(W.A)
size(W::AbstractWoodbury, d) = size(W.A, d)
axes(W::AbstractWoodbury) = axes(W.A)
axes(W::AbstractWoodbury, d) = _axes(W.A, d)
_axes(F::Factorization, dim::Int) = hasmethod(axes, Tuple{typeof(F), Int}) ? axes(F, dim) : Base.OneTo(size(F, dim))
_axes(A::AbstractArray, dim::Int) = axes(A, dim)

Base.Matrix(W::AbstractWoodbury{T}) where {T} = Matrix{T}(W)
Base.Matrix{T}(W::AbstractWoodbury) where {T} = convert(Matrix{T}, W)
Base.Array(W::AbstractWoodbury) = Matrix(W)
Base.Array{T}(W::AbstractWoodbury) where T = Matrix{T}(W)

convert(::Type{Matrix{T}}, W::AbstractWoodbury) where {T} = convert(Matrix{T}, Matrix(W.A) + W.U*W.C*W.V)

# This is a slow hack, but generally these matrices aren't sparse.
SparseArrays.sparse(W::AbstractWoodbury) = sparse(Matrix(W))

# Multiplication

*(W::AbstractWoodbury, B::AbstractMatrix)=W.A*B + W.U*(W.C*(W.V*B))

function *(W::AbstractWoodbury{T}, x::AbstractVector{S}) where {T,S}
    TS = Base.promote_op(LinearAlgebra.matprod, T, S)
    mul!(similar(x, TS, size(W,1)), W, x)
end

function mul!(dest, W::AbstractWoodbury, x::AbstractVector)
    # A reduced-allocation optimization (using temp storage for the multiplications)
    if W.tmpN1 !== nothing
        mul!(W.tmpN1, W.A, x)
        mul!(W.tmpk1, W.V, x)
        mul!(W.tmpk2, W.C, W.tmpk1)
        mul!(W.tmpN2, W.U, W.tmpk2)
        dest .= W.tmpN1 .+ W.tmpN2
    else
        dest .= W.A * x + W.U * (W.C * (W.V * x))
    end
    return dest
end

# Division

function _ldiv(W::AbstractWoodbury, R::AbstractVecOrMat)
    AinvR = W.A\R
    return AinvR - W.A\(W.U*(W.Cp*(W.V*AinvR)))
end

\(W::AbstractWoodbury, R::AbstractMatrix) = _ldiv(W, R)
\(W::AbstractWoodbury{T}, R::Matrix{Complex{T}}) where T<:Union{Float32, Float64} = _ldiv(W, R)  # ambiguity resolution
\(W::AbstractWoodbury, D::Diagonal) = _ldiv(W, D)

ldiv!(W::AbstractWoodbury, B::AbstractVector) = ldiv!(B, W, B)

function ldiv!(dest::AbstractVector, W::AbstractWoodbury, B::AbstractVector)
    @noinline throwdmm(W, B) = throw(DimensionMismatch("Vector length $(length(B)) must match matrix size $(size(W,1))"))

    length(B) == size(W, 1) || throwdmm(W, B)
    _ldiv!(dest, W, W.A, B)
    return dest
end

function _ldiv!(dest, W::AbstractWoodbury, A::Union{Factorization,Diagonal}, B)
    if W.tmpN1 !== nothing
        ldiv!(W.tmpN1, A, B)
        mul!(W.tmpk1, W.V, W.tmpN1)
        mul!(W.tmpk2, W.Cp, W.tmpk1)
        mul!(W.tmpN2, W.U, W.tmpk2)
        ldiv!(A, W.tmpN2)
        for i in eachindex(W.tmpN2)
            dest[i] = W.tmpN1[i] - W.tmpN2[i]
        end
    else
        copyto!(dest, _ldiv(W, B))
    end
    return dest
end
_ldiv!(dest, W, A, B) = _ldiv!(dest, W, defaultfactor(W, A), B)

defaultfactor(::AbstractWoodbury, A) = lu(A)

function det(W::AbstractWoodbury)
    ret = det(W.A)
    if !isempty(W.C)
        ret *= det(W.C)/det(W.Cp)
    end
    return ret
end
function logdet(W::AbstractWoodbury)
    ret = logdet(W.A)
    if !isempty(W.C)
        ret += logdet(W.C) - logdet(W.Cp)
    end
    return ret
end
function logabsdet(W::AbstractWoodbury)
    lad_A = logabsdet(W.A)
    isempty(W.C) && return lad_A
    lad_C = logabsdet(W.C)
    lad_Cp = logabsdet(W.Cp)
    lad = lad_A[1] + lad_C[1] - lad_Cp[1]
    s = lad_A[2] * lad_C[2] / lad_Cp[2]
    return lad, s
end

*(W::AbstractWoodbury, α::Real) = α*W
/(W::AbstractWoodbury, α::Real) = (1/α)*W
+(M::AbstractMatrix, W::AbstractWoodbury) = W + M

end
