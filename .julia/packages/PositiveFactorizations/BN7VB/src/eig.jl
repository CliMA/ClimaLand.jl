using LinearAlgebra: Eigen
import LinearAlgebra: eigen

function eigen(::Type{Positive{T}}, A::AbstractMatrix{T}, args...; tol=default_tol(A)) where {T}
    F = eigen(A, args...)
    for i = 1:size(A,1)
        tmp = abs(F.values[i])
        if tmp < tol
            tmp = one(tmp)
        end
        F.values[i] = tmp
    end
    F
end
eigen(::Type{Positive}, A::AbstractMatrix, args...; tol=default_tol(A)) = eigen(Positive{floattype(eltype(A))}, A, args...; tol=tol)
eigen(::Type{Positive{T}}, A::AbstractMatrix, args...; tol=default_tol(A)) where {T} = eigen(Positive{T}, convert(Matrix{T}, A), args...; tol=tol)
