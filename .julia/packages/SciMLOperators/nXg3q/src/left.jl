#
###
# left multiplication
###

for op in (:*, :\)
    @eval function Base.$op(u::AbstractVecOrMat, L::AbstractSciMLOperator)
        oper = u isa Transpose ? transpose : adjoint
        $op(oper(L), oper(u)) |> oper
    end
end

function LinearAlgebra.mul!(v::AbstractVecOrMat,
        u::AbstractVecOrMat,
        L::AbstractSciMLOperator)
    op = (u isa Transpose) | (v isa Transpose) ? transpose : adjoint
    mul!(op(v), op(L), op(u))
    v
end

function LinearAlgebra.mul!(v::AbstractVecOrMat,
        u::AbstractVecOrMat,
        L::AbstractSciMLOperator,
        α,
        β)
    op = (u isa Transpose) | (v isa Transpose) ? transpose : adjoint
    mul!(op(v), op(L), op(u), α, β)
    v
end

function LinearAlgebra.ldiv!(v::AbstractVecOrMat,
        u::AbstractVecOrMat,
        L::AbstractSciMLOperator)
    op = (u isa Transpose) | (v isa Transpose) ? transpose : adjoint
    ldiv!(op(v), op(L), op(u))
    v
end

function LinearAlgebra.ldiv!(u::AbstractVecOrMat, L::AbstractSciMLOperator)
    op = (u isa Transpose) ? transpose : adjoint
    ldiv!(op(L), op(u))
    u
end

###
# fallback wrappers
###

"""
$TYPEDEF
"""
struct AdjointOperator{T, LType} <: AbstractSciMLOperator{T}
    L::LType

    function AdjointOperator(L::AbstractSciMLOperator{T}) where {T}
        new{T, typeof(L)}(L)
    end
end

"""
$TYPEDEF
"""
struct TransposedOperator{T, LType} <: AbstractSciMLOperator{T}
    L::LType

    function TransposedOperator(L::AbstractSciMLOperator{T}) where {T}
        new{T, typeof(L)}(L)
    end
end

AbstractAdjointVecOrMat = Adjoint{T, <:AbstractVecOrMat} where {T}
AbstractTransposedVecOrMat = Transpose{T, <:AbstractVecOrMat} where {T}

has_adjoint(::AdjointOperator) = true
has_adjoint(L::TransposedOperator) = isreal(L) & has_adjoint(L.L)

islinear(L::AdjointOperator) = islinear(L.L)
islinear(L::TransposedOperator) = islinear(L.L)

Base.transpose(L::AdjointOperator) = conj(L.L)
Base.adjoint(L::TransposedOperator) = conj(L.L)

function Base.show(io::IO, L::AdjointOperator)
    show(io, L.L)
    print(io, "'")
end

function Base.show(io::IO, L::TransposedOperator)
    print(io, "transpose(")
    show(io, L.L)
    print(io, ")")
end

for (op, LType, VType) in ((:adjoint, :AdjointOperator, :AbstractAdjointVecOrMat),
    (:transpose, :TransposedOperator, :AbstractTransposedVecOrMat))
    # constructor
    @eval Base.$op(L::AbstractSciMLOperator) = $LType(L)

    @eval Base.convert(::Type{AbstractMatrix}, L::$LType) = $op(convert(AbstractMatrix,
        L.L))

    # traits
    @eval Base.size(L::$LType) = size(L.L) |> reverse
    @eval Base.resize!(L::$LType, n::Integer) = (resize!(L.L, n); L)
    @eval Base.$op(L::$LType) = L.L

    @eval getops(L::$LType) = (L.L,)

    # Copy method to avoid aliasing
    @eval Base.copy(L::$LType) = $LType(copy(L.L))

    @eval @forward $LType.L (
        # LinearAlgebra
        LinearAlgebra.issymmetric,
        LinearAlgebra.ishermitian,
        LinearAlgebra.isposdef,
        LinearAlgebra.opnorm,

        # SciML
        isconstant,
        has_mul,
        has_mul!,
        has_ldiv,
        has_ldiv!)

    @eval function cache_internals(L::$LType, u::AbstractVecOrMat)
        @reset L.L = cache_operator(L.L, reshape(u, size(L, 1)))
        L
    end

    # operator application
    @eval Base.:*(u::$VType, L::$LType) = $op(L.L * u.parent)
    @eval Base.:/(u::$VType, L::$LType) = $op(L.L \ u.parent)

    # v' ← u' * A'
    # v  ← A  * u
    @eval function LinearAlgebra.mul!(v::$VType, u::$VType, L::$LType)
        mul!(v.parent, L.L, u.parent)
        v
    end

    # v' ← α * (u' * A') + β * v'
    # v  ← α * (A  * u ) + β * v
    @eval function LinearAlgebra.mul!(v::$VType, u::$VType, L::$LType, α, β)
        mul!(v.parent, L.L, u.parent, α, β)
        v
    end

    # v' ← u' / A'
    # v  ← A  \ u
    @eval function LinearAlgebra.ldiv!(v::$VType, u::$VType, L::$LType)
        ldiv!(v.parent, L.L, u.parent)
        v
    end

    # u' ← u' / A'
    # u  ← A  \ u
    @eval function LinearAlgebra.ldiv!(u::$VType, L::$LType)
        ldiv!(L.L, u.parent)
        u
    end
end

# For AdjointOperator
# Out-of-place: v is action vector, u is update vector
function (L::AdjointOperator)(v::AbstractVecOrMat, u, p, t; kwargs...)
    # Adjoint operator applied to v means L.L' * v
    # For matrices: (A')v = (v'A)'
    # This means we need to compute L.L(v', u, p, t)' 
    # Update the operator first, then apply adjoint operator
    L_updated = update_coefficients(L.L, u, p, t; kwargs...)
    # (A')v = (v'A)' where v'A is computed by A'*v'
    return (L_updated' * v')'
end

# In-place: w is destination, v is action vector, u is update vector
function (L::AdjointOperator)(w::AbstractVecOrMat, v::AbstractVecOrMat, u, p, t; kwargs...)
    # Update the operator in-place
    update_coefficients!(L.L, u, p, t; kwargs...)
    # Use direct in-place multiplicatieon for adjoints
    mul!(w', v', L.L)
    return w
end

# In-place with scaling: w = α*(L*v) + β*w
function (L::AdjointOperator)(
        w::AbstractVecOrMat, v::AbstractVecOrMat, u, p, t, α, β; kwargs...)
    # Update the operator in-place
    update_coefficients!(L.L, u, p, t; kwargs...)
    mul!(w', v', L.L, α, β)
    return w
end

# For TransposedOperator
# Out-of-place
function (L::TransposedOperator)(v::AbstractVecOrMat, u, p, t; kwargs...)
    L_updated = update_coefficients(L.L, u, p, t; kwargs...)
    # (A^T)v = (v'A)' where v'A is computed by A'*v'
    return (L_updated' * v')'
end

# In-place
function (L::TransposedOperator)(
        w::AbstractVecOrMat, v::AbstractVecOrMat, u, p, t; kwargs...)
    update_coefficients!(L.L, u, p, t; kwargs...)
    mul!(w', v', L.L)
    return w
end

# In-place with scaling
function (L::TransposedOperator)(
        w::AbstractVecOrMat, v::AbstractVecOrMat, u, p, t, α, β; kwargs...)
    update_coefficients!(L.L, u, p, t; kwargs...)
    mul!(w', v', L.L, α, β)
    return w
end
#
