function hcat(x::Vector{T}) where T<:AbstractVecOrMat
    l    = length(x)
    if l == 0
        throw(ArgumentError("cannot flatten empty vector"))
    else
        x1   = x[1]
        m, n = size(x1, 1), size(x1, 2)
        B    = similar(x1, eltype(x1), (m, l*n))
        for i = 1:l
            B[:, (i - 1)*n + 1:i*n] = x[i]
        end
        return B
    end
end

if VERSION < v"1.3.0-alpha.115"
mul!(y::StridedVector{T},
     A::StridedMatrix{T},
     x::StridedVector{T},
     α::Number,
     β::Number) where {T<:BlasFloat} =
    BLAS.gemv!('N', convert(T, α), A, x, convert(T, β), y)

mul!(y::StridedVector{T},
     A::Adjoint{T,<:StridedMatrix{T}},
     x::StridedVector{T},
     α::Number,
     β::Number) where {T<:BlasReal} =
    BLAS.gemv!('T', convert(T, α), parent(A), x, convert(T, β), y)

mul!(y::StridedVector{T},
     A::Adjoint{T,<:StridedMatrix{T}},
     x::StridedVector{T},
     α::Number,
     β::Number) where {T<:BlasComplex} =
    BLAS.gemv!('C', convert(T, α), parent(A), x, convert(T, β), y)

function mul!(y::StridedVector, A::StridedMatrix, x::StridedVector, α::Number, β::Number)
    n = length(y)
    for i = 1:n
        y[i] *= β
    end
    for i = 1:n
        for l = 1:size(A,2)
            y[i] += α*A[i,l]*x[l]
        end
    end
    return y
end

function mul!(y::StridedVector, A::Adjoint{<:StridedMatrix}, x::StridedVector, α::Number, β::Number)
    n = length(y)
    for i = 1:n
        y[i] *= β
    end
    for i = 1:n
        for l = 1:size(A,2)
            y[i] += α*A[l,i]'*x[l]
        end
    end
    return y
end
end

function qr!(x::AbstractArray)
    nm = norm(x)
    rmul!(x, inv(nm))
    return x, nm
end
