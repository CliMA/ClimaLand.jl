struct SymWoodbury{T,AType,BType,DType,DpType} <: AbstractWoodbury{T}
    A::AType
    B::BType
    D::DType
    Dp::DpType
    tmpN1::Union{Vector{T}, Nothing}
    tmpN2::Union{Vector{T}, Nothing}
    tmpk1::Union{Vector{T}, Nothing}
    tmpk2::Union{Vector{T}, Nothing}

    SymWoodbury{T}(A, B, D, Dp, tmpN1, tmpN2, tmpk1, tmpk2) where {T} =
        new{T,typeof(A),typeof(B),typeof(D),typeof(Dp)}(A, B, D, Dp, tmpN1, tmpN2, tmpk1, tmpk2)
end

"""
    W = SymWoodbury(A, B, D; allocatetmp::Bool=false)

Represent a matrix of the form `W = A + BDBᵀ`, where `A` and `D` are symmetric.
Equations `Wx = b` will be solved using the
[Woodbury matrix identity](https://en.wikipedia.org/wiki/Woodbury_matrix_identity).

If your main goal is to solve equations, it can be advantageous to supply
`A` as a factorization (e.g., `SymWoodbury(cholesky(A), B, D)` if `A` is a positive
semidefinite matrix).

This checks that `A` is symmetric; you can elide the check by passing a `Symmetric` matrix
or factorization.

See also [Woodbury](@ref), where `allocatetmp` is explained.
"""
function SymWoodbury(A, B::AbstractVecOrMat, D; allocatetmp::Bool=false)
    @noinline throwdmm(B, D, A) = throw(DimensionMismatch("Sizes of B ($(size(B))) and/or D ($(size(D))) are inconsistent with A ($(size(A)))"))

    n = size(A, 1)
    k = size(B, 2)
    if size(A, 2) != n || size(B, 1) != n || size(D,1) != k || size(D,2) != k
        throwdmm(B, D, A)
    end
    if !isa(A, Cholesky) && !isa(A, BunchKaufman) && !isa(A, LDLt)
        issymmetric(A) || throw(ArgumentError("A must be symmetric"))
    end
    issymmetric(D) || throw(ArgumentError("D must be symmetric"))
    Dp = safeinv(safeinv(D) .+ B'*(A\B))
    # temporary space for allocation-free solver (vector RHS only)
    T = typeof(float(zero(eltype(A)) * zero(eltype(B)) * zero(eltype(D))))
    if allocatetmp
        tmpN1 = Vector{T}(undef, n)
        tmpN2 = Vector{T}(undef, n)
        tmpk1 = Vector{T}(undef, k)
        tmpk2 = Vector{T}(undef, k)
    else
        tmpN1 = tmpN2 = tmpk1 = tmpk2 = nothing
    end

    SymWoodbury{T}(A, B, D, Dp, tmpN1, tmpN2, tmpk1, tmpk2)
end

convert(::Type{W}, O::SymWoodbury) where {W<:Woodbury} = Woodbury(O.A, O.B, O.D, O.B')

convert(::Type{Matrix{T}}, W::SymWoodbury) where {T} = convert(Matrix{T}, symmetrize!(Matrix(W.A) + W.B*W.D*W.B'))

# This allows a common interface with Woodbury
function Base.getproperty(O::SymWoodbury, name::Symbol)
    # SymWoodbury names
    name === :A && return getfield(O, :A)
    name === :B && return getfield(O, :B)
    name === :D && return getfield(O, :D)
    name === :Dp && return getfield(O, :Dp)
    name === :tmpN1 && return getfield(O, :tmpN1)
    name === :tmpN2 && return getfield(O, :tmpN2)
    name === :tmpk1 && return getfield(O, :tmpk1)
    name === :tmpk2 && return getfield(O, :tmpk2)
    # Woodbury names
    name === :U && return getfield(O, :B)
    name === :C && return getfield(O, :D)
    name === :Cp && return getfield(O, :Dp)
    name === :V && return getfield(O, :B)'
    error("property ", name, " not defined")
end

inv_invD_BtX(invD, B, X) = safeinv(invD - B'*X)
inv_invD_BtX(invD, B::AbstractVector, X) = safeinv(invD - dot(B,X))

function calc_inv(A, B, D)
    W = safeinv(A)
    X = W*B
    Z = inv_invD_BtX(-safeinv(D), B, X)
    SymWoodbury(symmetrize!(W), X, symmetrize!(Z))
end

inv(O::SymWoodbury{T,AType,BType,DType}) where {T<:Any, AType<:Any, BType<:AbstractVector, DType<:Real} =
  calc_inv(O.A, O.B, O.D)

inv(O::SymWoodbury{T,AType,BType,DType}) where {T<:Any, AType<:Any, BType<:Any, DType<:AbstractMatrix} =
  calc_inv(O.A, O.B, O.D)

# D is typically small, so this is acceptable.
inv(O::SymWoodbury{T,AType,BType,DType}) where {T<:Any, AType<:Any, BType<:Any, DType<:SparseMatrixCSC} =
  calc_inv(O.A, O.B, Matrix(O.D))

defaultfactor(::SymWoodbury, A) = bunchkaufman(A, true)
defaultfactor(::SymWoodbury, A::SymTridiagonal) = ldlt(A)
defaultfactor(::SymWoodbury, A::SparseMatrixCSC) = lu(A)

"""
    partialInv(O)

Get the factors `(X, Z)` in `Y + XZXᵀ == inv(A + BDBᵀ)`
"""
function partialInv(O::SymWoodbury)
    X = (O.A)\O.B
    invD = -1*inv(O.D)
    Z = inv_invD_BtX(invD, O.B, X)
    return (X,Z)
end

+(O::SymWoodbury, M::SymWoodbury)    = SymWoodbury(O.A + M.A, [O.B M.B],
                                                   cat(O.D,M.D; dims=(1,2)) )
*(α::Real, O::SymWoodbury)           = SymWoodbury(α*O.A, O.B, α*O.D)
+(O::SymWoodbury, M::AbstractMatrix) = SymWoodbury(O.A + M, O.B, O.D)

Base.copy(O::SymWoodbury{T}) where {T} = SymWoodbury(copy(O.A), copy(O.B), copy(O.D))

function square(O::SymWoodbury)
    A  = O.A^2
    AB = O.A*O.B
    Z  = [(AB + O.B) (AB - O.B)]
    R  = symmetrize!(O.D*(O.B'*O.B)*O.D/4)
    D  = [ O.D/2 + R  -R
          -R          -O.D/2 + R ]
    SymWoodbury(A, Z, D)
end

"""
The product of two SymWoodbury matrices will generally be a Woodbury Matrix,
except when they are the same, i.e. the user writes A'A or A*A' or A*A.

Z(A + B*D*Bᵀ) = ZA + ZB*D*Bᵀ

This package will not support support left multiplication by a generic
matrix, to keep return types consistent.
"""
function *(O1::SymWoodbury, O2::SymWoodbury)
    if (O1 === O2)
        return square(O1)
    elseif O1.A == O2.A && O1.B == O2.B && O1.D == O2.D
        return square(O1)
    else
        throw(MethodError("ERROR: To multiply two non-identical SymWoodbury matrices, first convert to Woodbury."))
    end
end

conjm(O::SymWoodbury, M) = SymWoodbury(symmetrize!(M*O.A*M'), M*O.B, symmetrize!(O.D))

Base.getindex(O::SymWoodbury, I::UnitRange, I2::UnitRange) =
  SymWoodbury(O.A[I,I], O.B[I,:], O.D);

adjoint(O::SymWoodbury{T}) where T<:Real = O   # it's Hermitian
transpose(O::SymWoodbury) = O

issymmetric(::SymWoodbury) = true

function show(io::IO, mime::MIME"text/plain", W::SymWoodbury)
    println(io, "Symmetric Woodbury factorization:\nA:")
    show(io, mime, W.A)
    print(io, "\nB:\n")
    Base.print_matrix(IOContext(io,:compact=>true), W.B)
    print(io, "\nD: ", W.D)
end
function show(io::IO, W::SymWoodbury)
    print(io, "SymWoodbury{$(eltype(W))}(A=")
    show(io, W.A)
    print(io, ", B=")
    show(io, W.B)
    print(io, ", D=")
    show(io, W.D)
    print(io, ")")
end


function symmetrize!(A::AbstractMatrix)
    axs = axes(A)
    axs[1] == axs[2] || error("matrix must be square")
    for j in axs[2], i in j+1:last(axs[1])
        a = (A[i,j] + A[j,i])/2
        A[i,j] = A[j,i] = a
    end
    return A
end

symmetrize!(x::Real) = x
